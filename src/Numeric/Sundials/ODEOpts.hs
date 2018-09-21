module Numeric.Sundials.ODEOpts where

import           Data.Int (Int32)
import qualified Data.Vector.Storable as VS

import           Numeric.LinearAlgebra.HMatrix (Vector, Matrix)


type Jacobian = Double -> Vector Double -> Matrix Double

data ODEOpts = ODEOpts {
    maxNumSteps :: Int32
  , minStep     :: Double
  , maxFail     :: Int32
  } deriving (Read, Show, Eq, Ord)

data SundialsDiagnostics = SundialsDiagnostics {
    odeGetNumSteps               :: Int
  , odeGetNumStepAttempts        :: Int
  , odeGetNumRhsEvals_fe         :: Int
  , odeGetNumRhsEvals_fi         :: Int
  , odeGetNumLinSolvSetups       :: Int
  , odeGetNumErrTestFails        :: Int
  , odeGetNumNonlinSolvIters     :: Int
  , odeGetNumNonlinSolvConvFails :: Int
  , dlsGetNumJacEvals            :: Int
  , dlsGetNumRhsEvals            :: Int
  } deriving Show

