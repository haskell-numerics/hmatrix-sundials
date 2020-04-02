-- | Common infrastructure for CVode/ARKode
{-# OPTIONS_GHC -Wno-name-shadowing #-}
module Numeric.Sundials.Common where

import Foreign.C.Types
import Foreign.Ptr
import Foreign.Storable (peek, poke)
import Foreign.C.String
import Numeric.Sundials.Types
import qualified Numeric.Sundials.Foreign as T
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import Data.Maybe
import Numeric.LinearAlgebra.HMatrix hiding (Vector)
import GHC.Prim
import Control.Monad.IO.Class
import Control.Monad.Cont
import Control.Exception
import Katip
import Data.Aeson
import qualified Data.Text as T
import qualified Data.Text.Encoding as T
import qualified Data.ByteString as BS
import Control.Monad.Reader
import GHC.Generics (Generic)
import Foreign.ForeignPtr

-- | A collection of variables that we allocate on the Haskell side and
-- pass into the C code to be filled.
data CVars vec = CVars
  { c_diagnostics :: vec SunIndexType
    -- ^ Mutable vector to which we write diagnostic data while
    -- solving. Its size corresponds to the number of fields in
    -- 'SundialsDiagnostics'.
  , c_root_info :: vec CInt
    -- ^ Just a temporary vector (of the size equal to the number of event
    -- specs) that we use to get root info. Isn't used for output.
  , c_event_index :: vec CInt
    -- ^ For each event occurrence, this indicates which of the events
    -- occurred. Size: max_num_events.
  , c_event_time :: vec CDouble
    -- ^ For each event occurrence, this indicates the time of the
    -- occurrence. Size: max_num_events.
  , c_n_events :: vec CInt
    -- ^ Vector of size 1 that gives the total number of events occurred.
  , c_n_rows :: vec CInt
    -- ^ The total number of rows in the output matrix.
  , c_output_mat :: vec CDouble
    -- ^ The output matrix stored in the row-major order.
    -- Dimensions: (1 + dim) * (2 * max_events + nTs).
  , c_actual_event_direction :: vec CInt
    -- ^ Vector of size max_num_events that gives the direction of the
    -- occurred event.
  , c_local_error :: vec CDouble
    -- ^ Vector containing local error estimates. Size: the dimensionality
    -- of the system.
  , c_var_weight :: vec CDouble
    -- ^ Vector containing variable weights (derived from the tolerances).
    -- Size: the dimensionality of the system.
  , c_local_error_set :: vec CInt
    -- The flag (size 1) indicating whether c_local_error is filled with meaningful
    -- values. *Should be initialized with 0.*
  }

allocateCVars :: OdeProblem -> IO (CVars (VS.MVector RealWorld))
allocateCVars OdeProblem{..} = do
  let dim = VS.length odeInitCond
  c_diagnostics <- VSM.new 11
  c_root_info <- VSM.new $ V.length odeEvents
  c_event_index <- VSM.new odeMaxEvents
  c_event_time <- VSM.new odeMaxEvents
  c_actual_event_direction <- VSM.new odeMaxEvents
  c_n_events <- VSM.new 1
  c_n_rows <- VSM.new 1
  c_local_error <- VSM.new dim
  c_var_weight <- VSM.new dim
  c_local_error_set <- VSM.new 1
  VSM.write c_local_error_set 0 0
  c_output_mat <- VSM.new $
    (1 + dim) * (2 * odeMaxEvents + VS.length odeSolTimes)
  return CVars {..}

-- NB: the mutable CVars must not be used after this
freezeCVars :: CVars (VS.MVector RealWorld) -> IO (CVars VS.Vector)
freezeCVars CVars{..} = do
  c_diagnostics <- VS.unsafeFreeze c_diagnostics
  c_root_info <- VS.unsafeFreeze c_root_info
  c_event_index <- VS.unsafeFreeze c_event_index
  c_event_time <- VS.unsafeFreeze c_event_time
  c_actual_event_direction <- VS.unsafeFreeze c_actual_event_direction
  c_n_events <- VS.unsafeFreeze c_n_events
  c_n_rows <- VS.unsafeFreeze c_n_rows
  c_output_mat <- VS.unsafeFreeze c_output_mat
  c_local_error <- VS.unsafeFreeze c_local_error
  c_var_weight <- VS.unsafeFreeze c_var_weight
  c_local_error_set <- VS.unsafeFreeze c_local_error_set
  return CVars {..}

-- | Similar to 'CVars', except these are immutable values that are
-- accessed (read-only) by the C code and specify the system to be solved.
data CConsts = CConsts
  { c_dim :: SunIndexType -- ^ the dimensionality (number of variables/equations)
  , c_method :: CInt -- ^ the ODE method (specific to the solver)
  , c_n_sol_times :: CInt
  , c_sol_time :: VS.Vector CDouble
  , c_init_cond :: VS.Vector CDouble
  , c_rhs :: FunPtr OdeRhsCType
  , c_rhs_userdata :: Ptr UserData
  , c_rtol :: CDouble
  , c_atol :: VS.Vector CDouble
  , c_n_event_specs :: CInt
  , c_event_fn :: CDouble -> Ptr T.SunVector -> Ptr CDouble -> Ptr () -> IO CInt
  , c_apply_event
      :: CInt -- number of triggered events
      -> Ptr CInt -- event indices
      -> CDouble -- time
      -> Ptr T.SunVector -- y
      -> Ptr T.SunVector -- new y
      -> Ptr CInt -- (out) stop the solver?
      -> Ptr CInt -- (out) record the event?
      -> IO CInt
  , c_jac_set :: CInt
  , c_jac :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunMatrix
          -> Ptr () -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunVector
          -> IO CInt
  , c_requested_event_direction :: VS.Vector CInt
  , c_max_events :: CInt
  , c_minstep :: CDouble
  , c_fixedstep :: CDouble
  , c_max_n_steps :: SunIndexType
  , c_max_err_test_fails :: CInt
  , c_init_step_size_set :: CInt
  , c_init_step_size :: CDouble
  }

data Solver = CVode | ARKode
  deriving Show

class Method method where
  methodToInt :: method -> CInt
  methodSolver :: Solver

withCConsts
  :: Method method
  => ODEOpts method
  -> OdeProblem
  -> (CConsts -> IO r)
  -> IO r
withCConsts ODEOpts{..} OdeProblem{..} = runContT $ do
  let
    dim = VS.length c_init_cond
    c_init_cond = coerce odeInitCond
    c_dim = fromIntegral dim
    c_n_sol_times = fromIntegral . VS.length $ odeSolTimes
    c_sol_time = coerce odeSolTimes
    c_rtol = relTolerance odeTolerances
    c_atol = either (VS.replicate dim) id $ absTolerances odeTolerances
    c_minstep = coerce minStep
    c_fixedstep = coerce fixedStep
    c_max_n_steps = fromIntegral maxNumSteps
    c_max_err_test_fails = fromIntegral maxFail
    c_init_step_size_set = fromIntegral . fromEnum $ isJust initStep
    c_init_step_size = coerce . fromMaybe 0 $ initStep
    c_n_event_specs = fromIntegral $ V.length odeEvents
    c_requested_event_direction = V.convert $ V.map (directionToInt . eventDirection) odeEvents
    c_event_fn t y_ptr out_ptr _ptr = do
      y <- sunVecVals <$> peek y_ptr
      let vals = V.convert $ V.map (\ev -> coerce (eventCondition ev) t y) odeEvents
      -- FIXME: We should be able to use poke somehow
      T.vectorToC vals (fromIntegral c_n_event_specs) out_ptr
      return 0
    c_apply_event n_events event_indices_ptr t y_ptr y'_ptr stop_solver_ptr record_event_ptr = do
      event_indices <- vecFromPtr event_indices_ptr (fromIntegral n_events)
      -- Apparently there's no safe version of 
      y_vec <- peek y_ptr
      EventHandlerResult{..} <-
        odeEventHandler
          (coerce t :: Double)
          (coerce $ sunVecVals y_vec :: VS.Vector Double)
          (VS.map fromIntegral event_indices :: VS.Vector Int)
      poke y'_ptr $ SunVector
        { sunVecN = sunVecN y_vec
        , sunVecVals = coerce eventNewState
        }
      poke stop_solver_ptr . fromIntegral $ fromEnum eventStopSolver
      poke record_event_ptr . fromIntegral $ fromEnum eventRecord
      return 0
    c_max_events = fromIntegral odeMaxEvents
    c_jac_set = fromIntegral . fromEnum $ isJust odeJacobian
    c_jac t y _fy jacS _ptr _tmp1 _tmp2 _tmp3 = do
      case odeJacobian of
        Nothing   -> undefined
        Just jacI -> do j <- matrixToSunMatrix . jacI (coerce t) <$> (coerce $ sunVecVals <$> peek y)
                        poke jacS j
                        return 0
    c_method = methodToInt odeMethod

  (c_rhs, c_rhs_userdata) <-
    case odeRhs of
      OdeRhsC ptr u -> return (ptr, u)
      OdeRhsHaskell fun -> do
        let
          funIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr UserData -> IO CInt
          funIO t y f _ptr = do
            sv <- peek y
            r <- fun t (sunVecVals sv)
            poke f $ SunVector { sunVecN = sunVecN sv
                               , sunVecVals = r
                               }
            return 0
        funptr <- ContT $ bracket (mkOdeRhsC funIO) freeHaskellFunPtr
        return (funptr, nullPtr)
  return CConsts{..}

matrixToSunMatrix :: Matrix Double -> T.SunMatrix
matrixToSunMatrix m = T.SunMatrix { T.rows = nr, T.cols = nc, T.vals = vs }
  where
    nr = fromIntegral $ rows m
    nc = fromIntegral $ cols m
    -- FIXME: efficiency
    vs = VS.fromList $ map coerce $ concat $ toLists m

-- Contrary to the documentation, it appears that CVodeGetRootInfo
-- may use both 1 and -1 to indicate a root, depending on the
-- direction of the sign change. See near the end of cvRootfind.
intToDirection :: Integral d => d -> Maybe CrossingDirection
intToDirection d =
  case d of
    1  -> Just Upwards
    -1 -> Just Downwards
    _  -> Nothing

-- | Almost inverse of 'intToDirection'. Map 'Upwards' to 1, 'Downwards' to
-- -1, and 'AnyDirection' to 0.
directionToInt :: Integral d => CrossingDirection -> d
directionToInt d =
  case d of
    Upwards -> 1
    Downwards -> -1
    AnyDirection -> 0

foreign import ccall "wrapper"
  mkOdeRhsC :: OdeRhsCType -> IO (FunPtr OdeRhsCType)

assembleSolverResult
  :: OdeProblem
  -> CInt
  -> CVars VS.Vector
  -> IO (Either ErrorDiagnostics SundialsSolution)
assembleSolverResult OdeProblem{..} ret CVars{..} = do
  let
    dim = VS.length odeInitCond
    n_rows = fromIntegral . VS.head $ c_n_rows
    output_mat = coerce . reshape (dim + 1) . subVector 0 ((dim + 1) * n_rows) $ c_output_mat
    (local_errors, var_weights) =
      if c_local_error_set VS.! 0 == 0
        then (mempty, mempty)
        else coerce (c_local_error, c_var_weight)
    diagnostics = SundialsDiagnostics
      (fromIntegral $ c_diagnostics VS.!0)
      (fromIntegral $ c_diagnostics VS.!1)
      (fromIntegral $ c_diagnostics VS.!2)
      (fromIntegral $ c_diagnostics VS.!3)
      (fromIntegral $ c_diagnostics VS.!4)
      (fromIntegral $ c_diagnostics VS.!5)
      (fromIntegral $ c_diagnostics VS.!6)
      (fromIntegral $ c_diagnostics VS.!7)
      (fromIntegral $ c_diagnostics VS.!8)
      (fromIntegral $ c_diagnostics VS.!9)
      (toEnum . fromIntegral $ c_diagnostics VS.! 10)
  return $
    if ret == T.cV_SUCCESS
      then
        Right $ SundialsSolution
          { actualTimeGrid = extractTimeGrid output_mat
          , solutionMatrix = dropTimeGrid output_mat
          , diagnostics
          }
      else
        Left ErrorDiagnostics
          { partialResults = output_mat
          , errorCode = fromIntegral ret
          , errorEstimates = local_errors
          , varWeights = var_weights
          }
  where
    -- The time grid is the first column of the result matrix
    extractTimeGrid :: Matrix Double -> VS.Vector Double
    extractTimeGrid = head . toColumns
    dropTimeGrid :: Matrix Double -> Matrix Double
    dropTimeGrid = fromColumns . tail . toColumns

-- | The common solving logic between ARKode and CVode
solveCommon
  :: (Method method, Katip m)
  => (CConsts -> CVars (VS.MVector RealWorld) -> ReportErrorFn -> IO CInt)
      -- ^ the CVode/ARKode solving function; mostly inline-C code
  -> ODEOpts method
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveCommon solve_c opts problem@(OdeProblem{..})

  | VS.null odeInitCond = -- 0-dimensional (empty) system

    return . Right $ SundialsSolution
      { actualTimeGrid = odeSolTimes
      , solutionMatrix = (VS.length odeSolTimes >< 0) []
      , diagnostics = emptyDiagnostics
      }

  | otherwise = do

    report_error <- logWithKatip
    liftIO $ do -- the rest is in the IO monad
    vars <- allocateCVars problem
    ret <- withCConsts opts problem $ \consts ->
      solve_c consts vars report_error
    frozenVars <- freezeCVars vars
    assembleSolverResult problem ret frozenVars

-- | An auxiliary function to construct a storable vector from a C pointer
-- and length.
--
-- There doesn't seem to be a safe version of 'VS.unsafeFromForeignPtr0',
-- nor a way to clone an immutable vector, so we emulate it via an
-- intermediate mutable vector.
vecFromPtr
  :: VS.Storable a
  => Ptr a
  -> Int
  -> IO (VS.Vector a)
vecFromPtr ptr n = do
  fptr <- newForeignPtr_ ptr
  let mv = VSM.unsafeFromForeignPtr0 fptr n
  VS.freeze mv -- this does the copying and makes the whole thing safe

----------------------------------------------------------------------
--                           Logging
----------------------------------------------------------------------

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
    \err_code c_mod_name c_func_name c_msg _userdata ->
    -- See Note [CV_TOO_CLOSE]
    if err_code == T.cV_TOO_CLOSE then pure () else do
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
