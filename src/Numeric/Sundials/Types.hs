{-# LANGUAGE DeriveAnyClass, DeriveGeneric #-}
module Numeric.Sundials.Types
  ( OdeRhsCType
  , OdeRhs(..)
  , Jacobian
  , StepControl(..)
  , ODEOpts(..)
  , SundialsDiagnostics(..)
  , emptyDiagnostics
  , SundialsSolution(..)
  , EventInfo(..)
  , CrossingDirection(..)
  , EventSpec(..)
  , SunVector(..)
  , SunIndexType
  , SunRealType
  , sunContentLengthOffset
  , sunContentDataOffset
  )
  where

import           Data.Int (Int32)
import qualified Data.Vector.Storable as VS

import           Numeric.LinearAlgebra.HMatrix (Vector, Matrix)
import           Control.DeepSeq (NFData)
import           GHC.Generics (Generic)
import           Foreign.C.Types
import           Foreign.Ptr
import           Numeric.Sundials.Arkode (SunVector(..),
                                          SunIndexType, SunRealType,
                                          sunContentLengthOffset,
                                          sunContentDataOffset)

-- | The type of the C ODE RHS function.
type OdeRhsCType = CDouble -> Ptr SunVector -> Ptr SunVector -> Ptr () -> IO CInt

-- | The right-hand side of the ODE system.
--
-- Can be either a Haskell function or a pointer to a C function.
data OdeRhs
  = OdeRhsHaskell (CDouble -> VS.Vector CDouble -> VS.Vector CDouble)
  | OdeRhsC (FunPtr OdeRhsCType)

type Jacobian = Double -> Vector Double -> Matrix Double

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
-- [CVode](https://computation.llnl.gov/projects/sundials/cvode)
-- and [ARKode](https://computation.llnl.gov/projects/sundials/arkode)
-- allow the user to control the step size adjustment using
-- \(\eta^{rel}|y_i| + \eta^{abs}_i\). For compatibility with
-- [hmatrix-gsl](https://hackage.haskell.org/package/hmatrix-gsl),
-- tolerances for \(y\) and \(\dot{y}\) can be specified but the latter have no
-- effect.
data StepControl = X     Double Double -- ^ absolute and relative tolerance for \(y\); in GSL terms, \(a_{y} = 1\) and \(a_{dy/dt} = 0\); in ARKode terms, the \(\eta^{abs}_i\) are identical
                 | X'    Double Double -- ^ absolute and relative tolerance for \(\dot{y}\); in GSL terms, \(a_{y} = 0\) and \(a_{dy/dt} = 1\); in ARKode terms, the latter is treated as the relative tolerance for \(y\) so this is the same as specifying 'X' which may be entirely incorrect for the given problem
                 | XX'   Double Double Double Double -- ^ include both via relative tolerance
                                                     -- scaling factors \(a_y\), \(a_{{dy}/{dt}}\); in ARKode terms, the latter is ignored and \(\eta^{rel} = a_{y}\epsilon^{rel}\)
                 | ScXX' Double Double Double Double (Vector Double) -- ^ scale absolute tolerance of \(y_i\); in ARKode terms, \(a_{{dy}/{dt}}\) is ignored, \(\eta^{abs}_i = s_i \epsilon^{abs}\) and \(\eta^{rel} = a_{y}\epsilon^{rel}\)
  deriving (Eq, Ord, Show, Read)

data ODEOpts method = ODEOpts {
    maxNumSteps :: Int32
  , minStep     :: Double
  , maxFail     :: Int32
  , odeMethod   :: method
  , stepControl :: StepControl
  , initStep    :: Maybe Double
    -- ^ initial step size - by default, CVode
    -- estimates the initial step size to be the
    -- solution \(h\) of the equation
    -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
    -- \(\ddot{y}\) is an estimated value of the second
    -- derivative of the solution at \(t_0\)
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

emptyDiagnostics :: SundialsDiagnostics
emptyDiagnostics = SundialsDiagnostics 0 0 0 0 0 0 0 0 0 0

data SundialsSolution =
  SundialsSolution
  { actualTimeGrid :: VS.Vector Double    -- ^ actual time grid returned by the solver (with duplicated event times)
  , solutionMatrix :: Matrix Double       -- ^ matrix of solutions: each column is an unknwown
  , eventInfo      :: [EventInfo]         -- ^ event infos, as many items as triggered events during the simulation
  , diagnostics    :: SundialsDiagnostics -- ^ usual Sundials diagnostics
  }

data EventInfo =
  EventInfo
  { eventTime     :: !Double            -- ^ time at which event was triggered
  , eventIndex    :: !Int               -- ^ which index was triggered
  , rootDirection :: !CrossingDirection -- ^ in which direction ((+)->(-) or (-)->(+)) the root is crossed
  }
  deriving (Generic, Show, NFData)

-- | The direction in which a function should cross the x axis
data CrossingDirection = Upwards | Downwards | AnyDirection
  deriving (Generic, Eq, Show, NFData)

data EventSpec = EventSpec
  { eventCondition :: Double -> VS.Vector Double -> Double
  , eventDirection :: CrossingDirection
  , eventUpdate :: Double -> VS.Vector Double -> VS.Vector Double
  }
