{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE EmptyDataDecls #-}

module Numeric.Sundials.Arkode ( getDataFromContents
                               , putDataInContents
                               , cV_ADAMS
                               , cV_BDF
                               , vectorToC
                               , cV_SUCCESS
                               , cV_ROOT_RETURN
                               , SunIndexType
                               , SunRealType
                               , SunMatrix(..)
                               , SunVector(..)
                               , sunContentLengthOffset
                               , sunContentDataOffset
                               , hEUN_EULER_2_1_2
                               , bOGACKI_SHAMPINE_4_2_3
                               , aRK324L2SA_ERK_4_2_3
                               , zONNEVELD_5_3_4
                               , aRK436L2SA_ERK_6_3_4
                               , sAYFY_ABURUB_6_3_4
                               , cASH_KARP_6_4_5
                               , fEHLBERG_6_4_5
                               , dORMAND_PRINCE_7_4_5
                               , aRK548L2SA_ERK_8_4_5
                               , vERNER_8_5_6
                               , fEHLBERG_13_7_8
                               , sDIRK_2_1_2
                               , bILLINGTON_3_3_2
                               , tRBDF2_3_3_2
                               , kVAERNO_4_2_3
                               , aRK324L2SA_DIRK_4_2_3
                               , cASH_5_2_4
                               , cASH_5_3_4
                               , sDIRK_5_3_4
                               , kVAERNO_5_3_4
                               , aRK436L2SA_DIRK_6_3_4
                               , kVAERNO_7_4_5
                               , aRK548L2SA_DIRK_8_4_5
                               ) where

import           Foreign
import           Foreign.C.Types

import           Language.C.Types as CT

import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VM

import qualified Language.Haskell.TH as TH
import qualified Data.Map as Map
import           Language.C.Inline.Context

import qualified Data.Vector.Storable as V


#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_matrix.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_butcher_dirk.h>
#include <cvode/cvode.h>


data SunVector = SunVector { sunVecN    :: SunIndexType
                           , sunVecVals :: V.Vector CDouble
                           }

data SunMatrix = SunMatrix { rows :: CInt
                           , cols :: CInt
                           , vals :: V.Vector CDouble
                           }

type SunIndexType = #type sunindextype
type SunRealType = #type realtype

getMatrixDataFromContents :: Ptr SunMatrix -> IO SunMatrix
getMatrixDataFromContents ptr = do
  qtr <- getContentMatrixPtr ptr
  rs  <- getNRows qtr
  cs  <- getNCols qtr
  rtr <- getMatrixData qtr
  vs  <- vectorFromC (fromIntegral $ rs * cs) rtr
  return $ SunMatrix { rows = rs, cols = cs, vals = vs }

putMatrixDataFromContents :: SunMatrix -> Ptr SunMatrix -> IO ()
putMatrixDataFromContents mat ptr = do
  let rs = rows mat
      cs = cols mat
      vs = vals mat
  qtr <- getContentMatrixPtr ptr
  putNRows rs qtr
  putNCols cs qtr
  rtr <- getMatrixData qtr
  vectorToC vs (fromIntegral $ rs * cs) rtr

instance Storable SunVector where
  poke p v    = putDataInContents (sunVecVals v) (fromIntegral $ sunVecN v) p
  peek p      = do (l, v) <- getDataFromContents p
                   return $ SunVector { sunVecN = fromIntegral l
                                      , sunVecVals = v
                                      }
  sizeOf _    = error "sizeOf not supported for SunVector"
  alignment _ = error "alignment not supported for SunVector"

instance Storable SunMatrix where
  poke        = flip putMatrixDataFromContents
  peek        = getMatrixDataFromContents
  sizeOf _    = error "sizeOf not supported for SunMatrix"
  alignment _ = error "alignment not supported for SunMatrix"

vectorFromC :: Storable a => Int -> Ptr a -> IO (VS.Vector a)
vectorFromC len ptr = do
  ptr' <- newForeignPtr_ ptr
  VS.freeze $ VM.unsafeFromForeignPtr0 ptr' len

vectorToC :: Storable a => VS.Vector a -> Int -> Ptr a -> IO ()
vectorToC vec len ptr = do
  ptr' <- newForeignPtr_ ptr
  VS.copy (VM.unsafeFromForeignPtr0 ptr' len) vec

getDataFromContents :: Ptr SunVector -> IO (SunIndexType, VS.Vector CDouble)
getDataFromContents ptr = do
  qtr <- getContentPtr ptr
  rtr <- getData qtr
  len' <- getLength qtr
  v <- vectorFromC (fromIntegral len') rtr
  return (len', v)

putDataInContents :: VS.Vector CDouble -> Int -> Ptr SunVector -> IO ()
putDataInContents vec len ptr = do
  qtr <- getContentPtr ptr
  rtr <- getData qtr
  putLength (fromIntegral len) qtr
  vectorToC vec len rtr

#def typedef struct _generic_N_Vector SunVector;
#def typedef struct _N_VectorContent_Serial SunContent;

#def typedef struct _generic_SUNMatrix SunMatrix;
#def typedef struct _SUNMatrixContent_Dense SunMatrixContent;

sunContentLengthOffset :: Int
sunContentLengthOffset = #offset SunContent, length

sunContentDataOffset :: Int
sunContentDataOffset = #offset SunContent, data

getContentMatrixPtr :: Storable a => Ptr b -> IO a
getContentMatrixPtr ptr = (#peek SunMatrix, content) ptr

getNRows :: Ptr b -> IO CInt
getNRows ptr = (#peek SunMatrixContent, M) ptr
putNRows :: CInt -> Ptr b -> IO ()
putNRows nr ptr = (#poke SunMatrixContent, M) ptr nr

getNCols :: Ptr b -> IO CInt
getNCols ptr = (#peek SunMatrixContent, N) ptr
putNCols :: CInt -> Ptr b -> IO ()
putNCols nc ptr = (#poke SunMatrixContent, N) ptr nc

getMatrixData :: Storable a => Ptr b -> IO a
getMatrixData ptr = (#peek SunMatrixContent, data) ptr

getContentPtr :: Storable a => Ptr b -> IO a
getContentPtr ptr = (#peek SunVector, content) ptr

getData :: Storable a => Ptr b -> IO a
getData ptr = (#peek SunContent, data) ptr

getLength :: Ptr b -> IO SunIndexType
getLength ptr = (#peek SunContent, length) ptr

putLength :: SunIndexType -> Ptr b -> IO ()
putLength l ptr = (#poke SunContent, length) ptr l

cV_SUCCESS :: Int
cV_SUCCESS = #const CV_SUCCESS
cV_ROOT_RETURN :: Int
cV_ROOT_RETURN = #const CV_ROOT_RETURN

cV_ADAMS :: Int
cV_ADAMS = #const CV_ADAMS
cV_BDF :: Int
cV_BDF = #const CV_BDF

mIN_DIRK_NUM, mAX_DIRK_NUM :: Int
mIN_DIRK_NUM = #const MIN_DIRK_NUM
mAX_DIRK_NUM = #const MAX_DIRK_NUM

-- FIXME: We could just use inline-c instead

-- Butcher table accessors -- implicit
sDIRK_2_1_2 :: Int
sDIRK_2_1_2 = #const SDIRK_2_1_2
bILLINGTON_3_3_2 :: Int
bILLINGTON_3_3_2 = #const BILLINGTON_3_3_2
tRBDF2_3_3_2 :: Int
tRBDF2_3_3_2 = #const TRBDF2_3_3_2
kVAERNO_4_2_3 :: Int
kVAERNO_4_2_3 = #const KVAERNO_4_2_3
aRK324L2SA_DIRK_4_2_3 :: Int
aRK324L2SA_DIRK_4_2_3 = #const ARK324L2SA_DIRK_4_2_3
cASH_5_2_4 :: Int
cASH_5_2_4 = #const CASH_5_2_4
cASH_5_3_4 :: Int
cASH_5_3_4 = #const CASH_5_3_4
sDIRK_5_3_4 :: Int
sDIRK_5_3_4 = #const SDIRK_5_3_4
kVAERNO_5_3_4 :: Int
kVAERNO_5_3_4 = #const KVAERNO_5_3_4
aRK436L2SA_DIRK_6_3_4 :: Int
aRK436L2SA_DIRK_6_3_4 = #const ARK436L2SA_DIRK_6_3_4
kVAERNO_7_4_5 :: Int
kVAERNO_7_4_5 = #const KVAERNO_7_4_5
aRK548L2SA_DIRK_8_4_5 :: Int
aRK548L2SA_DIRK_8_4_5 = #const ARK548L2SA_DIRK_8_4_5

-- #define DEFAULT_DIRK_2          SDIRK_2_1_2
-- #define DEFAULT_DIRK_3          ARK324L2SA_DIRK_4_2_3
-- #define DEFAULT_DIRK_4          SDIRK_5_3_4
-- #define DEFAULT_DIRK_5          ARK548L2SA_DIRK_8_4_5

-- Butcher table accessors -- explicit
hEUN_EULER_2_1_2 :: Int
hEUN_EULER_2_1_2 = #const HEUN_EULER_2_1_2
bOGACKI_SHAMPINE_4_2_3 :: Int
bOGACKI_SHAMPINE_4_2_3 = #const BOGACKI_SHAMPINE_4_2_3
aRK324L2SA_ERK_4_2_3 :: Int
aRK324L2SA_ERK_4_2_3 = #const ARK324L2SA_ERK_4_2_3
zONNEVELD_5_3_4 :: Int
zONNEVELD_5_3_4 = #const ZONNEVELD_5_3_4
aRK436L2SA_ERK_6_3_4 :: Int
aRK436L2SA_ERK_6_3_4 = #const ARK436L2SA_ERK_6_3_4
sAYFY_ABURUB_6_3_4 :: Int
sAYFY_ABURUB_6_3_4 = #const SAYFY_ABURUB_6_3_4
cASH_KARP_6_4_5 :: Int
cASH_KARP_6_4_5 = #const CASH_KARP_6_4_5
fEHLBERG_6_4_5 :: Int
fEHLBERG_6_4_5 = #const FEHLBERG_6_4_5
dORMAND_PRINCE_7_4_5 :: Int
dORMAND_PRINCE_7_4_5 = #const DORMAND_PRINCE_7_4_5
aRK548L2SA_ERK_8_4_5 :: Int
aRK548L2SA_ERK_8_4_5 = #const ARK548L2SA_ERK_8_4_5
vERNER_8_5_6 :: Int
vERNER_8_5_6 = #const VERNER_8_5_6
fEHLBERG_13_7_8 :: Int
fEHLBERG_13_7_8 = #const FEHLBERG_13_7_8

-- #define DEFAULT_ERK_2           HEUN_EULER_2_1_2
-- #define DEFAULT_ERK_3           BOGACKI_SHAMPINE_4_2_3
-- #define DEFAULT_ERK_4           ZONNEVELD_5_3_4
-- #define DEFAULT_ERK_5           CASH_KARP_6_4_5
-- #define DEFAULT_ERK_6           VERNER_8_5_6
-- #define DEFAULT_ERK_8           FEHLBERG_13_7_8
