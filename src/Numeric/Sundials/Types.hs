module Numeric.Sundials.Types
  ( OdeProblem(..)
  , EventHandler
  , EventHandlerResult(..)
  , Tolerances(..)
  , OdeRhsCType
  , OdeRhs(..)
  , odeRhsPure
  , UserData
  , Jacobian
  , JacobianRepr(..)
  , ODEOpts(..)
  , SundialsDiagnostics(..)
  , ErrorDiagnostics(..)
  , emptyDiagnostics
  , SundialsSolution(..)
  , CrossingDirection(..)
  , EventSpec(..)
  , TimeEventSpec(..)
  , SunVector(..)
  , SunIndexType
  , SunRealType
  , sunContentLengthOffset
  , sunContentDataOffset
  , sunCtx
  )
  where

import           Data.Int (Int32)
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Map.Strict as Map
import qualified Language.Haskell.TH as TH

import           Numeric.LinearAlgebra.HMatrix (Vector, Matrix)
import           Control.DeepSeq (NFData)
import           Foreign.C.Types
import           Foreign.Ptr
import           Language.C.Types as CT
import           Language.C.Inline.Context
import           Numeric.Sundials.Foreign (SunVector(..), SunMatrix(..),
                                          SunIndexType, SunRealType,
                                          sunContentLengthOffset,
                                          sunContentDataOffset)
import GHC.Generics (Generic)

data EventHandlerResult = EventHandlerResult
  { eventStopSolver :: !Bool
    -- ^ should we stop the solver after handling this event?
  , eventRecord :: !Bool
    -- ^ should we record the state before and after the event in the ODE
    -- solution?
  , eventNewState :: !(VS.Vector Double)
    -- ^ the new state after the event has been applied
  }

type EventHandler
  =  Double -- ^ time
  -> VS.Vector Double -- ^ values of the variables
  -> VS.Vector Int
    -- ^ Vector of triggered event indices.
    -- If the vector is empty, this is a time-based event.
  -> IO EventHandlerResult

data OdeProblem = OdeProblem
  { odeEvents :: V.Vector EventSpec
    -- ^ The events that may occur, including the condition when they occur
    -- and how to update the state of the system when they do.
  , odeMaxEvents :: !Int
    -- ^ The maximal number of events that may occur. This is needed to
    -- allocate enough space to store the events. If more events occur, an
    -- error is returned.
  , odeEventHandler :: EventHandler -- ^ The event handler.
  , odeTimeBasedEvents :: TimeEventSpec
  , odeRhs :: OdeRhs
    -- ^ The right-hand side of the system: either a Haskell function or
    -- a pointer to a compiled function.
  , odeJacobian :: Maybe (Double -> Vector Double -> Matrix Double)
    -- ^ The optional Jacobian (the arguments are the time and the state
    -- vector).
  , odeInitCond :: VS.Vector Double
    -- ^ The initial conditions of the problem.
  , odeSolTimes :: VS.Vector Double
    -- ^ The requested solution times. The actual solution times may be
    -- larger if any events occurred.
  , odeTolerances :: Tolerances
    -- ^ How much error is tolerated in each variable.
  }

data Tolerances = Tolerances
  { relTolerance :: !CDouble
  , absTolerances :: Either CDouble (VS.Vector CDouble)
    -- ^ If 'Left', then the same tolerance is used for all variables.
    --
    -- If 'Right', the vector should contain one tolerance per variable.
  } deriving (Read, Show, Eq, Ord)

-- | The type of the C ODE RHS function.
type OdeRhsCType = CDouble -> Ptr SunVector -> Ptr SunVector -> Ptr UserData -> IO CInt

data UserData

-- | The right-hand side of the ODE system.
--
-- Can be either a Haskell function or a pointer to a C function.
data OdeRhs
  = OdeRhsHaskell (CDouble -> VS.Vector CDouble -> IO (VS.Vector CDouble))
  | OdeRhsC (FunPtr OdeRhsCType) (Ptr UserData)

-- | A version of 'OdeRhsHaskell' that accepts a pure function
odeRhsPure
  :: (CDouble -> VS.Vector CDouble -> VS.Vector CDouble)
  -> OdeRhs
odeRhsPure f = OdeRhsHaskell $ \t y -> return $ f t y

type Jacobian = Double -> Vector Double -> Matrix Double

data JacobianRepr
  = SparseJacobian !Int -- ^ sparse Jacobian with the given number of non-zero elements
  | DenseJacobian
  deriving (Read, Show, Eq)

data ODEOpts method = ODEOpts {
    maxNumSteps :: Int32
  , minStep     :: Double
  , fixedStep   :: Double
      -- ^ If this is greater than 0.0, then a fixed-size step is used.
      --
      -- This is only recommended for testing/debugging, not for production
      -- use.
      --
      -- Also, this only has effect for ARKode; using this with CVode will
      -- trigger an error.
  , maxFail     :: Int32
  , odeMethod   :: method
  , initStep    :: Maybe Double
    -- ^ initial step size - by default, CVode
    -- estimates the initial step size to be the
    -- solution \(h\) of the equation
    -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
    -- \(\ddot{y}\) is an estimated value of the second
    -- derivative of the solution at \(t_0\)
  , jacobianRepr :: JacobianRepr
    -- ^ use a sparse matrix to represent the Jacobian
    -- and a sparse linear solver for Newton iterations
  } deriving (Read, Show, Eq)

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
  , odeMaxEventsReached          :: Bool
  } deriving Show

emptyDiagnostics :: SundialsDiagnostics
emptyDiagnostics = SundialsDiagnostics 0 0 0 0 0 0 0 0 0 0 False

data SundialsSolution =
  SundialsSolution
  { actualTimeGrid :: VS.Vector Double    -- ^ actual time grid returned by the solver (with duplicated event times)
  , solutionMatrix :: Matrix Double       -- ^ matrix of solutions: each column is an unknwown
  , diagnostics    :: SundialsDiagnostics -- ^ usual Sundials diagnostics
  }

data ErrorDiagnostics = ErrorDiagnostics
  { errorCode :: !Int
    -- ^ The numeric error code. Mostly useless at this point, since it is
    -- set to 1 under most error conditions. See 'solveOdeC'.
  , errorEstimates :: !(VS.Vector Double)
    -- ^ The local error estimates as returned by @CVodeGetEstLocalErrors@.
    -- Either an empty vector, or has the same dimensionality as the state
    -- space.
  , varWeights :: !(VS.Vector Double)
    -- ^ The weights with which errors are combined, equal to @1 / (atol_i + y_i * rtol)@.
    -- Either an empty vector, or has the same dimensionality as the state
    -- space.
  , partialResults :: !(Matrix Double)
    -- ^ Partial solution of the ODE system, up until the moment when
    -- solving failed. Contains the time as its first column.
  } deriving Show

-- | The direction in which a function should cross the x axis
data CrossingDirection = Upwards | Downwards | AnyDirection
  deriving (Generic, Eq, Show, NFData)

-- | A time-based event, implemented as an action that returns the time of
-- the next time-based event.
--
-- If there's an additional condition attached to a time-based event, it
-- should be verified in the event handler.
--
-- The action is supposed to be stateful, and the state of the action
-- should be updated by the event handler so that after a given time-based
-- event is handled, the action starts returning the time of the next
-- unhandled time-based event.
--
-- If there is no next time-based event, the action should return +Inf.
newtype TimeEventSpec = TimeEventSpec (IO Double)

data EventSpec = EventSpec
  { eventCondition  :: Double -> VS.Vector Double -> Double
  , eventDirection  :: !CrossingDirection
  }

sunTypesTable :: Map.Map TypeSpecifier TH.TypeQ
sunTypesTable = Map.fromList
  [
    (TypeName "sunindextype", [t| SunIndexType |] )
  , (TypeName "SunVector",    [t| SunVector |] )
  , (TypeName "SunMatrix",    [t| SunMatrix |] )
  , (TypeName "UserData",     [t| UserData |] )
  ]

-- | Allows to map between Haskell and C types
sunCtx :: Context
sunCtx = mempty {ctxTypesTable = sunTypesTable}
