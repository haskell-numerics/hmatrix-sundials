{-# LANGUAGE DeriveAnyClass, DeriveGeneric, TemplateHaskell, OverloadedStrings #-}
module Numeric.Sundials.Types
  ( OdeProblem(..)
  , Tolerances(..)
  , OdeRhsCType
  , OdeRhs(..)
  , UserData
  , Jacobian
  , ODEOpts(..)
  , SundialsDiagnostics(..)
  , ErrorDiagnostics(..)
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
  , sunCtx
  , SundialsErrorContext(..)
  , ReportErrorFn
  , logWithKatip
  )
  where

import           Data.Int (Int32)
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Map.Strict as Map
import qualified Language.Haskell.TH as TH

import           Numeric.LinearAlgebra.HMatrix (Vector, Matrix)
import           Control.DeepSeq (NFData)
import           GHC.Generics (Generic)
import           Foreign.C.Types
import           Foreign.Ptr
import           Language.C.Types as CT
import           Language.C.Inline.Context
import           Numeric.Sundials.Foreign (SunVector(..), SunMatrix(..),
                                          SunIndexType, SunRealType,
                                          sunContentLengthOffset,
                                          sunContentDataOffset)

import qualified Data.Text as T
import Data.Aeson
import Katip
import Foreign.C.String
import qualified Data.Text.Encoding as T
import qualified Data.ByteString as BS
import Control.Monad.Reader

data OdeProblem = OdeProblem
  { odeEvents :: V.Vector EventSpec
    -- ^ The events that may occur, including the condition when they occur
    -- and how to update the state of the system when they do.
  , odeMaxEvents :: !Int
    -- ^ The maximal number of events that may occur. This is needed to
    -- allocate enough space to store the events. If more events occur, an
    -- error is returned.
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
  = OdeRhsHaskell (CDouble -> VS.Vector CDouble -> VS.Vector CDouble)
  | OdeRhsC (FunPtr OdeRhsCType) (Ptr UserData)

type Jacobian = Double -> Vector Double -> Matrix Double

data ODEOpts method = ODEOpts {
    maxNumSteps :: Int32
  , minStep     :: Double
  , maxFail     :: Int32
  , odeMethod   :: method
  , stepControl :: Tolerances
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
  , odeMaxEventsReached          :: Bool
  } deriving Show

emptyDiagnostics :: SundialsDiagnostics
emptyDiagnostics = SundialsDiagnostics 0 0 0 0 0 0 0 0 0 0 False

data SundialsSolution =
  SundialsSolution
  { actualTimeGrid :: VS.Vector Double    -- ^ actual time grid returned by the solver (with duplicated event times)
  , solutionMatrix :: Matrix Double       -- ^ matrix of solutions: each column is an unknwown
  , eventInfo      :: V.Vector EventInfo  -- ^ event infos, as many items as triggered events during the simulation
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
  { eventCondition  :: Double -> VS.Vector Double -> Double
  , eventDirection  :: !CrossingDirection
  , eventUpdate     :: Double -> VS.Vector Double -> VS.Vector Double
  , eventStopSolver :: !Bool
      -- ^ if 'eventUpdate' is 'Nothing', we should just stop the solver
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

-- | The Katip payload for logging Sundials errors
data SundialsErrorContext = SundialsErrorContext
  { sundialsErrorCode :: !Int
  , sundialsErrorModule :: !T.Text
  , sundialsErrorFunction :: !T.Text
  } deriving Generic
instance ToJSON SundialsErrorContext
instance ToObject SundialsErrorContext
instance LogItem SundialsErrorContext where
  payloadKeys _ _ = AllKeys

type ReportErrorFn =
  (  CInt    -- error code
  -> CString -- module name
  -> CString -- function name
  -> CString -- the message
  -> Ptr ()  -- user data (ignored)
  -> IO ()
  )

logWithKatip
  :: Katip m
  => m ReportErrorFn
logWithKatip = do
  log_env <- getLogEnv
  return $
    \err_code c_mod_name c_func_name c_msg _userdata -> do
    let
      toText :: CString -> IO T.Text
      toText = fmap T.decodeUtf8 . BS.packCString
    mod_name <- toText c_mod_name
    func_name <- toText c_func_name
    msg <- toText c_msg
    let
      severity :: Severity
      severity =
        if err_code <= 0
          then ErrorS
          else WarningS
      errCtx :: SundialsErrorContext
      errCtx = SundialsErrorContext
        { sundialsErrorCode = fromIntegral err_code
        , sundialsErrorModule = mod_name
        , sundialsErrorFunction = func_name
        }
    flip runReaderT log_env . unKatipT $ do
      logF errCtx "sundials" severity (logStr msg)
