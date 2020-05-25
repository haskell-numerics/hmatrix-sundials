{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE EmptyDataDecls #-}

module Numeric.Sundials.Foreign
  ( getDataFromContents
  , putDataInContents
  , cV_ADAMS
  , cV_BDF
  , vectorToC
  , cV_SUCCESS
  , cV_ROOT_RETURN
  , cV_TOO_CLOSE
  , SunIndexType
  , SunRealType
  , SunMatrix(..)
  , SparsePattern(..)
  , SparseMatrix(..)
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

import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VM
import Data.IORef
import Data.Coerce
import Control.Monad
import Control.Exception
import Text.Printf (printf)

#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_matrix.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <cvode/cvode.h>


data SunVector = SunVector
  { sunVecN    :: SunIndexType
  , sunVecVals :: VS.Vector CDouble
  }

data SunMatrix = SunMatrix
  { rows :: SunIndexType
  , cols :: SunIndexType
  , vals :: VS.Vector CDouble
    -- ^ matrix entries in the column-major order
  }

  -- | A sparse pattern: a column-wise stored matrix that has 0 for zero
  -- entries and 1 for non-zero entries.
newtype SparsePattern = SparsePattern (VS.Vector Int8)
  deriving Show

-- | A sparse matrix
data SparseMatrix = SparseMatrix SparsePattern SunMatrix

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

instance Storable SparseMatrix where
  poke p (SparseMatrix (SparsePattern spat) SunMatrix{..}) = do
    content_ptr <- getContentMatrixPtr p
    sparse_type :: CInt <- #{peek SunMatrixContentSparse, sparsetype} content_ptr
    unless (sparse_type == cSC_MAT) $
      throwIO . ErrorCall $ printf "Sparse SUNMatrix poke: only CSC matrices are supported, got %d" (fromIntegral sparse_type :: Int)
    #{poke SunMatrixContentSparse, M} content_ptr rows
    #{poke SunMatrixContentSparse, N} content_ptr cols
    indexvals :: Ptr SunIndexType <- #{peek SunMatrixContentSparse, indexvals} content_ptr
    indexptrs :: Ptr SunIndexType <- #{peek SunMatrixContentSparse, indexptrs} content_ptr
    data_ :: Ptr SunRealType <- #{peek SunMatrixContentSparse, data} content_ptr
    max_entries :: SunIndexType <- #{peek SunMatrixContentSparse, NNZ} content_ptr

    when (VS.any (`notElem` [0,1]) spat) $
      throwIO $ ErrorCall $ "Illegal sparse pattern: " ++ show spat
    when (max_entries < VS.sum (VS.map fromIntegral spat)) $
      throwIO $ ErrorCall $ "Not enough space in sparse matrix for the sparse pattern"

    out_ix_ref <- newIORef 0
    forM_ [0 .. cols - 1] $ \cur_col -> do
      readIORef out_ix_ref >>= \out_ix ->
        pokeElemOff indexptrs (fromIntegral cur_col) (fromIntegral out_ix)
      forM_ [0 .. rows - 1] $ \cur_row -> do
        let
          in_ix = fromIntegral $ cur_row + cur_col * rows
          non_zero = (spat VS.! in_ix) /= 0
          val = vals VS.! in_ix
        when non_zero $ do
          out_ix <- readIORef out_ix_ref
          writeIORef out_ix_ref $! out_ix+1
          pokeElemOff data_ out_ix (coerce val)
          pokeElemOff indexvals out_ix (fromIntegral cur_row)
    readIORef out_ix_ref >>= \out_ix ->
      pokeElemOff indexptrs (fromIntegral cols) (fromIntegral out_ix)


  peek        = error "peek not supported for a sparse SunMatrix"
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
#def typedef struct _SUNMatrixContent_Sparse SunMatrixContentSparse;

sunContentLengthOffset :: Int
sunContentLengthOffset = #offset SunContent, length

sunContentDataOffset :: Int
sunContentDataOffset = #offset SunContent, data

getContentMatrixPtr :: Storable a => Ptr b -> IO a
getContentMatrixPtr ptr = (#peek SunMatrix, content) ptr

getNRows :: Ptr b -> IO SunIndexType
getNRows ptr = (#peek SunMatrixContent, M) ptr
putNRows :: SunIndexType -> Ptr b -> IO ()
putNRows nr ptr = (#poke SunMatrixContent, M) ptr nr

getNCols :: Ptr b -> IO SunIndexType
getNCols ptr = (#peek SunMatrixContent, N) ptr
putNCols :: SunIndexType -> Ptr b -> IO ()
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

cV_SUCCESS :: CInt
cV_SUCCESS = #const CV_SUCCESS
cV_ROOT_RETURN :: CInt
cV_ROOT_RETURN = #const CV_ROOT_RETURN
cV_TOO_CLOSE :: CInt
cV_TOO_CLOSE = #const CV_TOO_CLOSE

cSC_MAT :: CInt
cSC_MAT = #const CSC_MAT

cV_ADAMS :: CInt
cV_ADAMS = #const CV_ADAMS
cV_BDF :: CInt
cV_BDF = #const CV_BDF

-- Butcher table accessors -- implicit
sDIRK_2_1_2 :: CInt
sDIRK_2_1_2 = #const SDIRK_2_1_2
bILLINGTON_3_3_2 :: CInt
bILLINGTON_3_3_2 = #const BILLINGTON_3_3_2
tRBDF2_3_3_2 :: CInt
tRBDF2_3_3_2 = #const TRBDF2_3_3_2
kVAERNO_4_2_3 :: CInt
kVAERNO_4_2_3 = #const KVAERNO_4_2_3
aRK324L2SA_DIRK_4_2_3 :: CInt
aRK324L2SA_DIRK_4_2_3 = #const ARK324L2SA_DIRK_4_2_3
cASH_5_2_4 :: CInt
cASH_5_2_4 = #const CASH_5_2_4
cASH_5_3_4 :: CInt
cASH_5_3_4 = #const CASH_5_3_4
sDIRK_5_3_4 :: CInt
sDIRK_5_3_4 = #const SDIRK_5_3_4
kVAERNO_5_3_4 :: CInt
kVAERNO_5_3_4 = #const KVAERNO_5_3_4
aRK436L2SA_DIRK_6_3_4 :: CInt
aRK436L2SA_DIRK_6_3_4 = #const ARK436L2SA_DIRK_6_3_4
kVAERNO_7_4_5 :: CInt
kVAERNO_7_4_5 = #const KVAERNO_7_4_5
aRK548L2SA_DIRK_8_4_5 :: CInt
aRK548L2SA_DIRK_8_4_5 = #const ARK548L2SA_DIRK_8_4_5

-- #define DEFAULT_DIRK_2          SDIRK_2_1_2
-- #define DEFAULT_DIRK_3          ARK324L2SA_DIRK_4_2_3
-- #define DEFAULT_DIRK_4          SDIRK_5_3_4
-- #define DEFAULT_DIRK_5          ARK548L2SA_DIRK_8_4_5

-- Butcher table accessors -- explicit
hEUN_EULER_2_1_2 :: CInt
hEUN_EULER_2_1_2 = #const HEUN_EULER_2_1_2
bOGACKI_SHAMPINE_4_2_3 :: CInt
bOGACKI_SHAMPINE_4_2_3 = #const BOGACKI_SHAMPINE_4_2_3
aRK324L2SA_ERK_4_2_3 :: CInt
aRK324L2SA_ERK_4_2_3 = #const ARK324L2SA_ERK_4_2_3
zONNEVELD_5_3_4 :: CInt
zONNEVELD_5_3_4 = #const ZONNEVELD_5_3_4
aRK436L2SA_ERK_6_3_4 :: CInt
aRK436L2SA_ERK_6_3_4 = #const ARK436L2SA_ERK_6_3_4
sAYFY_ABURUB_6_3_4 :: CInt
sAYFY_ABURUB_6_3_4 = #const SAYFY_ABURUB_6_3_4
cASH_KARP_6_4_5 :: CInt
cASH_KARP_6_4_5 = #const CASH_KARP_6_4_5
fEHLBERG_6_4_5 :: CInt
fEHLBERG_6_4_5 = #const FEHLBERG_6_4_5
dORMAND_PRINCE_7_4_5 :: CInt
dORMAND_PRINCE_7_4_5 = #const DORMAND_PRINCE_7_4_5
aRK548L2SA_ERK_8_4_5 :: CInt
aRK548L2SA_ERK_8_4_5 = #const ARK548L2SA_ERK_8_4_5
vERNER_8_5_6 :: CInt
vERNER_8_5_6 = #const VERNER_8_5_6
fEHLBERG_13_7_8 :: CInt
fEHLBERG_13_7_8 = #const FEHLBERG_13_7_8

-- #define DEFAULT_ERK_2           HEUN_EULER_2_1_2
-- #define DEFAULT_ERK_3           BOGACKI_SHAMPINE_4_2_3
-- #define DEFAULT_ERK_4           ZONNEVELD_5_3_4
-- #define DEFAULT_ERK_5           CASH_KARP_6_4_5
-- #define DEFAULT_ERK_6           VERNER_8_5_6
-- #define DEFAULT_ERK_8           FEHLBERG_13_7_8
