{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}

-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Sundials.ARKode.ODE
-- Copyright   :  Dominic Steinitz 2018,
--                Novadiscovery 2018
-- License     :  BSD
-- Maintainer  :  Dominic Steinitz
-- Stability   :  provisional
--
-- Solution of ordinary differential equation (ODE) initial value problems.
-- See <https://computation.llnl.gov/projects/sundials/sundials-software> for more detail.
--
-- A simple example:
--
-- <<diagrams/brusselator.png#diagram=brusselator&height=400&width=500>>
--
-- @
-- import           Numeric.Sundials.ARKode.ODE
-- import           Numeric.LinearAlgebra
--
-- import           Plots as P
-- import qualified Diagrams.Prelude as D
-- import           Diagrams.Backend.Rasterific
--
-- brusselator :: Double -> [Double] -> [Double]
-- brusselator _t x = [ a - (w + 1) * u + v * u * u
--                    , w * u - v * u * u
--                    , (b - w) / eps - w * u
--                    ]
--   where
--     a = 1.0
--     b = 3.5
--     eps = 5.0e-6
--     u = x !! 0
--     v = x !! 1
--     w = x !! 2
--
-- lSaxis :: [[Double]] -> P.Axis B D.V2 Double
-- lSaxis xs = P.r2Axis &~ do
--   let ts = xs!!0
--       us = xs!!1
--       vs = xs!!2
--       ws = xs!!3
--   P.linePlot' $ zip ts us
--   P.linePlot' $ zip ts vs
--   P.linePlot' $ zip ts ws
--
-- main = do
--   let res1 = odeSolve brusselator [1.2, 3.1, 3.0] (fromList [0.0, 0.1 .. 10.0])
--   renderRasterific "diagrams/brusselator.png"
--                    (D.dims2D 500.0 500.0)
--                    (renderAxis $ lSaxis $ [0.0, 0.1 .. 10.0]:(toLists $ tr res1))
-- @
--
-- With Sundials ARKode, it is possible to retrieve the Butcher
-- tableau for the solver. FIXME: Not available just now and hopefully
-- normal service will be resumed soon.
--
-- @
-- import           Numeric.Sundials.ARKode.ODE
-- import           Numeric.LinearAlgebra
--
-- import           Data.List (intercalate)
--
-- import           Text.PrettyPrint.HughesPJClass
--
--
-- butcherTableauTex :: ButcherTable -> String
-- butcherTableauTex (ButcherTable m c b b2) =
--   render $
--   vcat [ text ("\n\\begin{array}{c|" ++ (concat $ replicate n "c") ++ "}")
--        , us
--        , text "\\hline"
--        , text bs <+> text "\\\\"
--        , text b2s <+> text "\\\\"
--        , text "\\end{array}"
--        ]
--   where
--     n = rows m
--     rs = toLists m
--     ss = map (\r -> intercalate " & " $ map show r) rs
--     ts = zipWith (\i r -> show i ++ " & " ++ r) (toList c) ss
--     us = vcat $ map (\r -> text r <+> text "\\\\") ts
--     bs  = " & " ++ (intercalate " & " $ map show $ toList b)
--     b2s = " & " ++ (intercalate " & " $ map show $ toList b2)
--
-- main :: IO ()
-- main = do
--
--   let res = butcherTable (SDIRK_2_1_2 undefined)
--   putStrLn $ show res
--   putStrLn $ butcherTableauTex res
--
--   let resA = butcherTable (KVAERNO_4_2_3 undefined)
--   putStrLn $ show resA
--   putStrLn $ butcherTableauTex resA
--
--   let resB = butcherTable (SDIRK_5_3_4 undefined)
--   putStrLn $ show resB
--   putStrLn $ butcherTableauTex resB
-- @
--
--  Using the code above from the examples gives
--
-- KVAERNO_4_2_3
--
-- \[
-- \begin{array}{c|cccc}
-- 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
-- 0.871733043 & 0.4358665215 & 0.4358665215 & 0.0 & 0.0 \\
-- 1.0 & 0.490563388419108 & 7.3570090080892e-2 & 0.4358665215 & 0.0 \\
-- 1.0 & 0.308809969973036 & 1.490563388254106 & -1.235239879727145 & 0.4358665215 \\
-- \hline
--  & 0.308809969973036 & 1.490563388254106 & -1.235239879727145 & 0.4358665215 \\
--  & 0.490563388419108 & 7.3570090080892e-2 & 0.4358665215 & 0.0 \\
-- \end{array}
-- \]
--
-- SDIRK_2_1_2
--
-- \[
-- \begin{array}{c|cc}
-- 1.0 & 1.0 & 0.0 \\
-- 0.0 & -1.0 & 1.0 \\
-- \hline
--  & 0.5 & 0.5 \\
--  & 1.0 & 0.0 \\
-- \end{array}
-- \]
--
-- SDIRK_5_3_4
--
-- \[
-- \begin{array}{c|ccccc}
-- 0.25 & 0.25 & 0.0 & 0.0 & 0.0 & 0.0 \\
-- 0.75 & 0.5 & 0.25 & 0.0 & 0.0 & 0.0 \\
-- 0.55 & 0.34 & -4.0e-2 & 0.25 & 0.0 & 0.0 \\
-- 0.5 & 0.2727941176470588 & -5.036764705882353e-2 & 2.7573529411764705e-2 & 0.25 & 0.0 \\
-- 1.0 & 1.0416666666666667 & -1.0208333333333333 & 7.8125 & -7.083333333333333 & 0.25 \\
-- \hline
--  & 1.0416666666666667 & -1.0208333333333333 & 7.8125 & -7.083333333333333 & 0.25 \\
--  & 1.2291666666666667 & -0.17708333333333334 & 7.03125 & -7.083333333333333 & 0.0 \\
-- \end{array}
-- \]
-----------------------------------------------------------------------------
module Numeric.Sundials.ARKode.ODE ( odeSolve
                                   , odeSolveV
                                   , odeSolveVWith
                                   , odeSolveVWith'
                                   , ODEMethod(..)
                                   , StepControl(..)
                                   ) where

import qualified Language.C.Inline as C
import qualified Language.C.Inline.Unsafe as CU

import           Data.Monoid ((<>))
import           Data.Maybe (isJust)

import           Foreign.C.Types (CDouble, CInt, CLong)
import           Foreign.Ptr (Ptr)
import           Foreign.Storable (poke, peek)

import qualified Data.Vector.Storable as V

import           Data.Coerce (coerce)
import           System.IO.Unsafe (unsafePerformIO)
import           GHC.Generics (C1, Constructor, (:+:)(..), D1, Rep, Generic, M1(..),
                               from, conName)

import           Numeric.LinearAlgebra.Devel (createVector)

import           Numeric.LinearAlgebra.HMatrix (Vector, Matrix, toList, rows,
                                                cols, toLists, size, reshape)

import           Numeric.Sundials.ODEOpts (ODEOpts(..), Jacobian, SundialsDiagnostics(..))
import qualified Numeric.Sundials.Arkode as T
import           Numeric.Sundials.Arkode (sDIRK_2_1_2,
                                          bILLINGTON_3_3_2,
                                          tRBDF2_3_3_2,
                                          kVAERNO_4_2_3,
                                          aRK324L2SA_DIRK_4_2_3,
                                          cASH_5_2_4,
                                          cASH_5_3_4,
                                          sDIRK_5_3_4,
                                          kVAERNO_5_3_4,
                                          aRK436L2SA_DIRK_6_3_4,
                                          kVAERNO_7_4_5,
                                          aRK548L2SA_DIRK_8_4_5,
                                          hEUN_EULER_2_1_2,
                                          bOGACKI_SHAMPINE_4_2_3,
                                          aRK324L2SA_ERK_4_2_3,
                                          zONNEVELD_5_3_4,
                                          aRK436L2SA_ERK_6_3_4,
                                          sAYFY_ABURUB_6_3_4,
                                          cASH_KARP_6_4_5,
                                          fEHLBERG_6_4_5,
                                          dORMAND_PRINCE_7_4_5,
                                          aRK548L2SA_ERK_8_4_5,
                                          vERNER_8_5_6,
                                          fEHLBERG_13_7_8)


C.context (C.baseCtx <> C.vecCtx <> C.funCtx <> T.sunCtx)

C.include "<stdlib.h>"
C.include "<stdio.h>"
C.include "<math.h>"
C.include "<arkode/arkode_arkstep.h>"                 -- prototypes for ARKODE fcts., consts.
C.include "<arkode/arkode_erkstep.h>"
C.include "<nvector/nvector_serial.h>"        -- serial N_Vector types, fcts., macros
C.include "<sunmatrix/sunmatrix_dense.h>"     -- access to dense SUNMatrix
C.include "<sunlinsol/sunlinsol_dense.h>"     -- access to dense SUNLinearSolver
C.include "<sundials/sundials_types.h>"       -- definition of type realtype
C.include "<sundials/sundials_math.h>"
C.include "../../../helpers.h"
C.include "Numeric/Sundials/Arkode_hsc.h"


-- | Stepping functions
data ODEMethod = SDIRK_2_1_2            Jacobian
               | SDIRK_2_1_2'
               | BILLINGTON_3_3_2       Jacobian
               | BILLINGTON_3_3_2'
               | TRBDF2_3_3_2           Jacobian
               | TRBDF2_3_3_2'
               | KVAERNO_4_2_3          Jacobian
               | KVAERNO_4_2_3'
               | ARK324L2SA_DIRK_4_2_3  Jacobian
               | ARK324L2SA_DIRK_4_2_3'
               | CASH_5_2_4             Jacobian
               | CASH_5_2_4'
               | CASH_5_3_4             Jacobian
               | CASH_5_3_4'
               | SDIRK_5_3_4            Jacobian
               | SDIRK_5_3_4'
               | KVAERNO_5_3_4          Jacobian
               | KVAERNO_5_3_4'
               | ARK436L2SA_DIRK_6_3_4  Jacobian
               | ARK436L2SA_DIRK_6_3_4'
               | KVAERNO_7_4_5          Jacobian
               | KVAERNO_7_4_5'
               | ARK548L2SA_DIRK_8_4_5  Jacobian
               | ARK548L2SA_DIRK_8_4_5'
               | HEUN_EULER_2_1_2         Jacobian
               | HEUN_EULER_2_1_2'
               | BOGACKI_SHAMPINE_4_2_3   Jacobian
               | BOGACKI_SHAMPINE_4_2_3'
               | ARK324L2SA_ERK_4_2_3     Jacobian
               | ARK324L2SA_ERK_4_2_3'
               | ZONNEVELD_5_3_4          Jacobian
               | ZONNEVELD_5_3_4'
               | ARK436L2SA_ERK_6_3_4     Jacobian
               | ARK436L2SA_ERK_6_3_4'
               | SAYFY_ABURUB_6_3_4       Jacobian
               | SAYFY_ABURUB_6_3_4'
               | CASH_KARP_6_4_5          Jacobian
               | CASH_KARP_6_4_5'
               | FEHLBERG_6_4_5         Jacobian
               | FEHLBERG_6_4_5'
               | DORMAND_PRINCE_7_4_5     Jacobian
               | DORMAND_PRINCE_7_4_5'
               | ARK548L2SA_ERK_8_4_5     Jacobian
               | ARK548L2SA_ERK_8_4_5'
               | VERNER_8_5_6            Jacobian
               | VERNER_8_5_6'
               | FEHLBERG_13_7_8         Jacobian
               | FEHLBERG_13_7_8'
  deriving Generic

constrName :: (HasConstructor (Rep a), Generic a)=> a -> String
constrName = genericConstrName . from

class HasConstructor (f :: * -> *) where
  genericConstrName :: f x -> String

instance HasConstructor f => HasConstructor (D1 c f) where
  genericConstrName (M1 x) = genericConstrName x

instance (HasConstructor x, HasConstructor y) => HasConstructor (x :+: y) where
  genericConstrName (L1 l) = genericConstrName l
  genericConstrName (R1 r) = genericConstrName r

instance Constructor c => HasConstructor (C1 c f) where
  genericConstrName x = conName x

instance Show ODEMethod where
  show x = constrName x

-- FIXME: We can probably do better here with generics
getMethod :: ODEMethod -> Int
getMethod (SDIRK_2_1_2 _)            = sDIRK_2_1_2
getMethod (SDIRK_2_1_2')             = sDIRK_2_1_2
getMethod (BILLINGTON_3_3_2 _)       = bILLINGTON_3_3_2
getMethod (BILLINGTON_3_3_2')        = bILLINGTON_3_3_2
getMethod (TRBDF2_3_3_2 _)           = tRBDF2_3_3_2
getMethod (TRBDF2_3_3_2')            = tRBDF2_3_3_2
getMethod (KVAERNO_4_2_3  _)         = kVAERNO_4_2_3
getMethod (KVAERNO_4_2_3')           = kVAERNO_4_2_3
getMethod (ARK324L2SA_DIRK_4_2_3 _)  = aRK324L2SA_DIRK_4_2_3
getMethod (ARK324L2SA_DIRK_4_2_3')   = aRK324L2SA_DIRK_4_2_3
getMethod (CASH_5_2_4 _)             = cASH_5_2_4
getMethod (CASH_5_2_4')              = cASH_5_2_4
getMethod (CASH_5_3_4 _)             = cASH_5_3_4
getMethod (CASH_5_3_4')              = cASH_5_3_4
getMethod (SDIRK_5_3_4 _)            = sDIRK_5_3_4
getMethod (SDIRK_5_3_4')             = sDIRK_5_3_4
getMethod (KVAERNO_5_3_4 _)          = kVAERNO_5_3_4
getMethod (KVAERNO_5_3_4')           = kVAERNO_5_3_4
getMethod (ARK436L2SA_DIRK_6_3_4 _)  = aRK436L2SA_DIRK_6_3_4
getMethod (ARK436L2SA_DIRK_6_3_4')   = aRK436L2SA_DIRK_6_3_4
getMethod (KVAERNO_7_4_5 _)          = kVAERNO_7_4_5
getMethod (KVAERNO_7_4_5')           = kVAERNO_7_4_5
getMethod (ARK548L2SA_DIRK_8_4_5 _)  = aRK548L2SA_DIRK_8_4_5
getMethod (ARK548L2SA_DIRK_8_4_5')   = aRK548L2SA_DIRK_8_4_5
getMethod (HEUN_EULER_2_1_2 _)       = hEUN_EULER_2_1_2
getMethod (HEUN_EULER_2_1_2')        = hEUN_EULER_2_1_2
getMethod (BOGACKI_SHAMPINE_4_2_3 _) = bOGACKI_SHAMPINE_4_2_3
getMethod (BOGACKI_SHAMPINE_4_2_3')  = bOGACKI_SHAMPINE_4_2_3
getMethod (ARK324L2SA_ERK_4_2_3 _)   = aRK324L2SA_ERK_4_2_3
getMethod (ARK324L2SA_ERK_4_2_3')    = aRK324L2SA_ERK_4_2_3
getMethod (ZONNEVELD_5_3_4 _)        = zONNEVELD_5_3_4
getMethod (ZONNEVELD_5_3_4')         = zONNEVELD_5_3_4
getMethod (ARK436L2SA_ERK_6_3_4 _)   = aRK436L2SA_ERK_6_3_4
getMethod (ARK436L2SA_ERK_6_3_4')    = aRK436L2SA_ERK_6_3_4
getMethod (SAYFY_ABURUB_6_3_4 _)     = sAYFY_ABURUB_6_3_4
getMethod (SAYFY_ABURUB_6_3_4')      = sAYFY_ABURUB_6_3_4
getMethod (CASH_KARP_6_4_5 _)        = cASH_KARP_6_4_5
getMethod (CASH_KARP_6_4_5')         = cASH_KARP_6_4_5
getMethod (FEHLBERG_6_4_5 _)         = fEHLBERG_6_4_5
getMethod (FEHLBERG_6_4_5' )         = fEHLBERG_6_4_5
getMethod (DORMAND_PRINCE_7_4_5 _)   = dORMAND_PRINCE_7_4_5
getMethod (DORMAND_PRINCE_7_4_5')    = dORMAND_PRINCE_7_4_5
getMethod (ARK548L2SA_ERK_8_4_5 _)   = aRK548L2SA_ERK_8_4_5
getMethod (ARK548L2SA_ERK_8_4_5')    = aRK548L2SA_ERK_8_4_5
getMethod (VERNER_8_5_6 _)           = vERNER_8_5_6
getMethod (VERNER_8_5_6')            = vERNER_8_5_6
getMethod (FEHLBERG_13_7_8 _)        = fEHLBERG_13_7_8
getMethod (FEHLBERG_13_7_8')         = fEHLBERG_13_7_8

getJacobian :: ODEMethod -> Maybe Jacobian
getJacobian (SDIRK_2_1_2 j)            = Just j
getJacobian (BILLINGTON_3_3_2 j)       = Just j
getJacobian (TRBDF2_3_3_2 j)           = Just j
getJacobian (KVAERNO_4_2_3  j)         = Just j
getJacobian (ARK324L2SA_DIRK_4_2_3 j)  = Just j
getJacobian (CASH_5_2_4 j)             = Just j
getJacobian (CASH_5_3_4 j)             = Just j
getJacobian (SDIRK_5_3_4 j)            = Just j
getJacobian (KVAERNO_5_3_4 j)          = Just j
getJacobian (ARK436L2SA_DIRK_6_3_4 j)  = Just j
getJacobian (KVAERNO_7_4_5 j)          = Just j
getJacobian (ARK548L2SA_DIRK_8_4_5 j)  = Just j
getJacobian (HEUN_EULER_2_1_2 j)       = Just j
getJacobian (BOGACKI_SHAMPINE_4_2_3 j) = Just j
getJacobian (ARK324L2SA_ERK_4_2_3 j)   = Just j
getJacobian (ZONNEVELD_5_3_4 j)        = Just j
getJacobian (ARK436L2SA_ERK_6_3_4 j)   = Just j
getJacobian (SAYFY_ABURUB_6_3_4 j)     = Just j
getJacobian (CASH_KARP_6_4_5 j)        = Just j
getJacobian (FEHLBERG_6_4_5 j)         = Just j
getJacobian (DORMAND_PRINCE_7_4_5 j)   = Just j
getJacobian (ARK548L2SA_ERK_8_4_5 j)   = Just j
getJacobian (VERNER_8_5_6 j)           = Just j
getJacobian (FEHLBERG_13_7_8 j)        = Just j
getJacobian _                          = Nothing

-- | A version of 'odeSolveVWith' with reasonable default step control.
odeSolveV
    :: ODEMethod
    -> Maybe Double      -- ^ initial step size - by default, ARKode
                         -- estimates the initial step size to be the
                         -- solution \(h\) of the equation
                         -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
                         -- \(\ddot{y}\) is an estimated value of the
                         -- second derivative of the solution at \(t_0\)
    -> Double            -- ^ absolute tolerance for the state vector
    -> Double            -- ^ relative tolerance for the state vector
    -> (Double -> Vector Double -> Vector Double) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
    -> Vector Double     -- ^ initial conditions
    -> Vector Double     -- ^ desired solution times
    -> Matrix Double     -- ^ solution
odeSolveV meth hi epsAbs epsRel f y0 ts =
  odeSolveVWith meth (X epsAbs epsRel) hi g y0 ts
    where
      g t x0 = coerce $ f t x0

-- | A version of 'odeSolveV' with reasonable default parameters and
-- system of equations defined using lists. FIXME: we should say
-- something about the fact we could use the Jacobian but don't for
-- compatibility with hmatrix-gsl.
odeSolve :: (Double -> [Double] -> [Double]) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
         -> [Double]                         -- ^ initial conditions
         -> Vector Double                    -- ^ desired solution times
         -> Matrix Double                    -- ^ solution
odeSolve f y0 ts =
  -- FIXME: These tolerances are different from the ones in GSL
  odeSolveVWith SDIRK_5_3_4' (XX' 1.0e-6 1.0e-10 1 1)  Nothing g (V.fromList y0) (V.fromList $ toList ts)
  where
    g t x0 = V.fromList $ f t (V.toList x0)

odeSolveVWith ::
  ODEMethod
  -> StepControl
  -> Maybe Double -- ^ initial step size - by default, ARKode
                  -- estimates the initial step size to be the
                  -- solution \(h\) of the equation
                  -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
                  -- \(\ddot{y}\) is an estimated value of the second
                  -- derivative of the solution at \(t_0\)
  -> (Double -> V.Vector Double -> V.Vector Double) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector Double                     -- ^ Initial conditions
  -> V.Vector Double                     -- ^ Desired solution times
  -> Matrix Double                       -- ^ Error code or solution
odeSolveVWith method control initStepSize f y0 tt =
  case odeSolveVWith' opts method control initStepSize f y0 tt of
    Left  (c, _v) -> error $ show c -- FIXME
    Right (v, _d) -> v
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod   = error "ARKode: unexpected use of ODEOpts.odeMethod"
                   , stepControl = error "ARKode: unexpected use of ODEOpts.stepControl"
                   , initStep    = error "ARKode: unexpected use of ODEOpts.initStep"
                   }

odeSolveVWith' ::
  ODEOpts
  -> ODEMethod
  -> StepControl
  -> Maybe Double -- ^ initial step size - by default, ARKode
                  -- estimates the initial step size to be the
                  -- solution \(h\) of the equation
                  -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
                  -- \(\ddot{y}\) is an estimated value of the second
                  -- derivative of the solution at \(t_0\)
  -> (Double -> V.Vector Double -> V.Vector Double) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector Double                     -- ^ Initial conditions
  -> V.Vector Double                     -- ^ Desired solution times
  -> Either (Matrix Double, Int) (Matrix Double, SundialsDiagnostics) -- ^ Error code or solution
odeSolveVWith' opts method control initStepSize f y0 tt =
  case solveOdeC (fromIntegral $ maxFail opts)
                 (fromIntegral $ maxNumSteps opts) (coerce $ minStep opts)
                 (fromIntegral $ getMethod method) (coerce initStepSize) jacH (scise control)
                 (coerce f) (coerce y0) (coerce tt) of
    Left  (v, c) -> Left  (reshape l (coerce v), fromIntegral c)
    Right (v, d) -> Right (reshape l (coerce v), d)
  where
    l = size y0
    scise (X aTol rTol)                          = coerce (V.replicate l aTol, rTol)
    scise (X' aTol rTol)                         = coerce (V.replicate l aTol, rTol)
    scise (XX' aTol rTol yScale _yDotScale)      = coerce (V.replicate l aTol, yScale * rTol)
    -- FIXME; Should we check that the length of ss is correct?
    scise (ScXX' aTol rTol yScale _yDotScale ss) = coerce (V.map (* aTol) ss, yScale * rTol)
    jacH = fmap (\g t v -> matrixToSunMatrix $ g (coerce t) (coerce v)) $
           getJacobian method
    matrixToSunMatrix m = T.SunMatrix { T.rows = nr, T.cols = nc, T.vals = vs }
      where
        nr = fromIntegral $ rows m
        nc = fromIntegral $ cols m
        -- FIXME: efficiency
        vs = V.fromList $ map coerce $ concat $ toLists m

solveOdeC ::
  CInt ->
  CLong ->
  CDouble ->
  CInt ->
  Maybe CDouble ->
  (Maybe (CDouble -> V.Vector CDouble -> T.SunMatrix)) ->
  (V.Vector CDouble, CDouble) ->
  (CDouble -> V.Vector CDouble -> V.Vector CDouble) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector CDouble -- ^ Initial conditions
  -> V.Vector CDouble -- ^ Desired solution times
  -> Either (V.Vector CDouble, CInt) (V.Vector CDouble, SundialsDiagnostics) -- ^ Partial solution and error code or
                                                                             -- solution and diagnostics
solveOdeC maxErrTestFails maxNumSteps_ minStep_ method initStepSize
          jacH (aTols, rTol) fun f0 ts = unsafePerformIO $ do

  let isInitStepSize :: CInt
      isInitStepSize = fromIntegral $ fromEnum $ isJust initStepSize
      ss :: CDouble
      ss = case initStepSize of
             -- It would be better to put an error message here but
             -- inline-c seems to evaluate this even if it is never
             -- used :(
             Nothing -> 0.0
             Just x  -> x

  let dim = V.length f0
      nEq :: CLong
      nEq = fromIntegral dim
      nTs :: CInt
      nTs = fromIntegral $ V.length ts
  quasiMatrixRes <- createVector ((fromIntegral dim) * (fromIntegral nTs))
  qMatMut <- V.thaw quasiMatrixRes
  diagnostics :: V.Vector CLong <- createVector 10 -- FIXME
  diagMut <- V.thaw diagnostics
  -- We need the types that sundials expects. These are tied together
  -- in 'CLangToHaskellTypes'. FIXME: The Haskell type is currently empty!
  let funIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr () -> IO CInt
      funIO t y f _ptr = do
        sv <- peek y
        poke f $ T.SunVector { T.sunVecN = T.sunVecN sv
                             , T.sunVecVals = fun t (T.sunVecVals sv)
                             }
        [CU.exp| int{ 0 } |]
  let isJac :: CInt
      isJac = fromIntegral $ fromEnum $ isJust jacH
      jacIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunMatrix ->
               Ptr () -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunVector ->
               IO CInt
      jacIO t y _fy jacS _ptr _tmp1 _tmp2 _tmp3 = do
        case jacH of
          Nothing   -> error "Numeric.Sundials.ARKode.ODE: Jacobian not defined"
          Just jacI -> do j <- jacI t <$> (T.sunVecVals <$> peek y)
                          poke jacS j
                          -- FIXME: I don't understand what this comment means
                          -- Unsafe since the function will be called many times.
                          [CU.exp| int{ 0 } |]

  res <- [C.block| int {
                         /* general problem variables */

                         int flag;                  /* reusable error-checking flag                 */
                         int i, j;                  /* reusable loop indices                        */
                         N_Vector y = NULL;         /* empty vector for storing solution            */
                         N_Vector tv = NULL;        /* empty vector for storing absolute tolerances */
                         SUNMatrix A = NULL;        /* empty matrix for linear solver               */
                         SUNLinearSolver LS = NULL; /* empty linear solver object                   */
                         void *arkode_mem = NULL;   /* empty ARKode memory structure                */
                         realtype t;
                         long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;

                         /* general problem parameters */

                         realtype T0 = RCONST(($vec-ptr:(double *ts))[0]); /* initial time              */
                         sunindextype NEQ = $(sunindextype nEq);             /* number of dependent vars. */

                         /* Initialize data structures */

                         y = N_VNew_Serial(NEQ); /* Create serial vector for solution */
                         if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
                         /* Specify initial condition */
                         for (i = 0; i < NEQ; i++) {
                           NV_Ith_S(y,i) = ($vec-ptr:(double *f0))[i];
                         };

                         tv = N_VNew_Serial(NEQ); /* Create serial vector for absolute tolerances */
                         if (check_flag((void *)tv, "N_VNew_Serial", 0)) return 1;
                         /* Specify tolerances */
                         for (i = 0; i < NEQ; i++) {
                           NV_Ith_S(tv,i) = ($vec-ptr:(double *aTols))[i];
                         };

                         /* Call ARKStepCreate to initialize the ARK timestepper module and */
                         /* specify the right-hand side function in y'=f(t,y), the inital time */
                         /* T0, and the initial dependent variable vector y.  Note: since this */
                         /* problem is fully implicit, we set f_E to NULL and f_I to f. */

                         /* Here we use the C types defined in helpers.h which tie up with */
                         /* the Haskell types defined in CLangToHaskellTypes               */
                         if ($(int method) < MIN_DIRK_NUM) {
                           arkode_mem = ARKStepCreate($fun:(int (* funIO) (double t, SunVector y[], SunVector dydt[], void * params)), NULL, T0, y);
                           if (check_flag((void *)arkode_mem, "ARKStepCreate", 0)) return 1;
                             } else {
                           arkode_mem = ARKStepCreate(NULL, $fun:(int (* funIO) (double t, SunVector y[], SunVector dydt[], void * params)), T0, y);
                           if (check_flag(&flag, "ARKStepCreate", 0)) return 1;
                         }

                         flag = ARKStepSetMinStep(arkode_mem, $(double minStep_));
                         if (check_flag(&flag, "ARKStepSetMinStep", 1)) return 1;
                         flag = ARKStepSetMaxNumSteps(arkode_mem, $(long int maxNumSteps_));
                         if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;
                         flag = ARKStepSetMaxErrTestFails(arkode_mem, $(int maxErrTestFails));
                         if (check_flag(&flag, "ARKStepSetMaxErrTestFails", 1)) return 1;

                         /* Set routines */
                         flag = ARKStepSVtolerances(arkode_mem, $(double rTol), tv);
                         if (check_flag(&flag, "ARKStepSVtolerances", 1)) return 1;

                         /* Initialize dense matrix data structure and solver */
                         A = SUNDenseMatrix(NEQ, NEQ);
                         if (check_flag((void *)A, "SUNDenseMatrix", 0)) return 1;
                         LS = SUNDenseLinearSolver(y, A);
                         if (check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return 1;

                         /* Attach matrix and linear solver */
                         flag = ARKStepSetLinearSolver(arkode_mem, LS, A);
                         if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;

                         /* Set the initial step size if there is one */
                         if ($(int isInitStepSize)) {
                           /* FIXME: We could check if the initial step size is 0 */
                           /* or even NaN and then throw an error                 */
                           flag = ARKStepSetInitStep(arkode_mem, $(double ss));
                           if (check_flag(&flag, "ARKStepSetInitStep", 1)) return 1;
                         }

                         /* Set the Jacobian if there is one */
                         if ($(int isJac)) {
                           flag = ARKStepSetJacFn(arkode_mem, $fun:(int (* jacIO) (double t, SunVector y[], SunVector fy[], SunMatrix Jac[], void * params, SunVector tmp1[], SunVector tmp2[], SunVector tmp3[])));
                           if (check_flag(&flag, "ARKStepSetJacFn", 1)) return 1;
                         }

                         /* Store initial conditions */
                         for (j = 0; j < NEQ; j++) {
                           ($vec-ptr:(double *qMatMut))[0 * $(int nTs) + j] = NV_Ith_S(y,j);
                         }

                         /* Explicitly set the method */
                         if ($(int method) >= MIN_DIRK_NUM) {
                           flag = ARKStepSetTableNum(arkode_mem, $(int method), -1);
                           if (check_flag(&flag, "ARKStepSetTableNum", 1)) return 1;
                         } else {
                           flag = ARKStepSetTableNum(arkode_mem, -1, $(int method));
                           if (check_flag(&flag, "ERKStepSetTableNum", 1)) return 1;
                         }

                         /* Main time-stepping loop: calls ARKStep to perform the integration */
                         /* Stops when the final time has been reached                       */
                         for (i = 1; i < $(int nTs); i++) {

                           flag = ARKStepEvolve(arkode_mem, ($vec-ptr:(double *ts))[i], y, &t, ARK_NORMAL); /* call integrator */
                           if (check_flag(&flag, "ARKStep solver failure, stopping integration", 1)) return 1;

                           /* Store the results for Haskell */
                           for (j = 0; j < NEQ; j++) {
                             ($vec-ptr:(double *qMatMut))[i * NEQ + j] = NV_Ith_S(y,j);
                           }
                         }

                         /* Get some final statistics on how the solve progressed */

                         flag = ARKStepGetNumSteps(arkode_mem, &nst);
                         check_flag(&flag, "ARKStepGetNumSteps", 1);
                         ($vec-ptr:(long int *diagMut))[0] = nst;

                         flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
                         check_flag(&flag, "ARKStepGetNumStepAttempts", 1);
                         ($vec-ptr:(long int *diagMut))[1] = nst_a;

                         flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
                         check_flag(&flag, "ARKStepGetNumRhsEvals", 1);
                         ($vec-ptr:(long int *diagMut))[2] = nfe;
                         ($vec-ptr:(long int *diagMut))[3] = nfi;

                         flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
                         check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1);
                         ($vec-ptr:(long int *diagMut))[4] = nsetups;

                         flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
                         check_flag(&flag, "ARKStepGetNumErrTestFails", 1);
                         ($vec-ptr:(long int *diagMut))[5] = netf;

                         flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
                         check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1);
                         ($vec-ptr:(long int *diagMut))[6] = nni;

                         flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
                         check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1);
                         ($vec-ptr:(long int *diagMut))[7] = ncfn;

                         flag = ARKStepGetNumJacEvals(arkode_mem, &nje);
                         check_flag(&flag, "ARKStepGetNumJacEvals", 1);
                         ($vec-ptr:(long int *diagMut))[8] = nje;

                         flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfeLS);
                         check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1);
                         ($vec-ptr:(long int *diagMut))[9] = nfeLS;

                         /* Clean up and return */
                         N_VDestroy(y);            /* Free y vector          */
                         N_VDestroy(tv);           /* Free tv vector         */
                         ARKStepFree(&arkode_mem);  /* Free integrator memory */
                         SUNLinSolFree(LS);        /* Free linear solver     */
                         SUNMatDestroy(A);         /* Free A matrix          */

                         return flag;
                       } |]
  preD <- V.freeze diagMut
  let d = SundialsDiagnostics (fromIntegral $ preD V.!0)
                              (fromIntegral $ preD V.!1)
                              (fromIntegral $ preD V.!2)
                              (fromIntegral $ preD V.!3)
                              (fromIntegral $ preD V.!4)
                              (fromIntegral $ preD V.!5)
                              (fromIntegral $ preD V.!6)
                              (fromIntegral $ preD V.!7)
                              (fromIntegral $ preD V.!8)
                              (fromIntegral $ preD V.!9)
  m <- V.freeze qMatMut
  if res == 0
    then do
      return $ Right (m, d)
    else do
      return $ Left  (m, res)

-- | Adaptive step-size control
-- functions.
--
-- [GSL](https://www.gnu.org/software/gsl/doc/html/ode-initval.html#adaptive-step-size-control)
-- allows the user to control the step size adjustment using
-- \(D_i = \epsilon^{abs}s_i + \epsilon^{rel}(a_{y} |y_i| + a_{dy/dt} h |\dot{y}_i|)\) where
-- \(\epsilon^{abs}\) is the required absolute error, \(\epsilon^{rel}\)
-- is the required relative error, \(s_i\) is a vector of scaling
-- factors, \(a_{y}\) is a scaling factor for the solution \(y\) and
-- \(a_{dydt}\) is a scaling factor for the derivative of the solution \(dy/dt\).
--
-- [ARKode](https://computation.llnl.gov/projects/sundials/arkode)
-- allows the user to control the step size adjustment using
-- \(\eta^{rel}|y_i| + \eta^{abs}_i\). For compatibility with
-- [hmatrix-gsl](https://hackage.haskell.org/package/hmatrix-gsl),
-- tolerances for \(y\) and \(\dot{y}\) can be specified but the latter have no
-- effect.
data StepControl = X     Double Double -- ^ absolute and relative tolerance for \(y\); in GSL terms, \(a_{y} = 1\) and \(a_{dy/dt} = 0\); in ARKode terms, the \(\eta^{abs}_i\) are identical
                 | X'    Double Double -- ^ absolute and relative tolerance for \(\dot{y}\); in GSL terms, \(a_{y} = 0\) and \(a_{dy/dt} = 1\); in ARKode terms, the latter is treated as the relative tolerance for \(y\) so this is the same as specifying 'X' which may be entirely incorrect for the given problem
                 | XX'   Double Double Double Double -- ^ include both via relative tolerance
                                                     -- scaling factors \(a_y\), \(a_{{dy}/{dt}}\); in ARKode terms, the latter is ignored and \(\eta^{rel} = a_{y}\epsilon^{rel}\)
                 | ScXX' Double Double Double Double (Vector Double) -- ^ scale absolute tolerance of \(y_i\); in ARKode terms, \(a_{{dy}/{dt}}\) is ignored, \(\eta^{abs}_i = s_i \epsilon^{abs}\) and \(\eta^{rel} = a_{y}\epsilon^{rel}\)
