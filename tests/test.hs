{-# LANGUAGE RecordWildCards, ScopedTypeVariables, OverloadedStrings,
             ViewPatterns, ImplicitParams, OverloadedLists, RankNTypes,
             ExistentialQuantification, LambdaCase, NumDecimals, NamedFieldPuns,
             TypeApplications, DeriveGeneric, StandaloneDeriving, FlexibleInstances #-}
{-# OPTIONS_GHC -Wno-orphans #-}

import Prelude hiding (quot, showList)

import Test.Tasty
import Test.Tasty.Golden
import Test.Tasty.Golden.Advanced
import Test.Tasty.HUnit

import Numeric.Sundials

import Numeric.LinearAlgebra as L hiding ((<.>), (<>))
import qualified Data.Vector.Storable as VS
import qualified Data.Vector as V
import Data.List
import Data.Maybe
import Data.Foldable
import Data.IORef
import Control.Monad
import Control.Exception
import Katip
import System.IO
import System.FilePath
import Text.Printf (printf)
import Data.Aeson
import Data.Aeson.Encode.Pretty
import GHC.Generics

----------------------------------------------------------------------
--                            Helpers
----------------------------------------------------------------------

emptyOdeProblem :: OdeProblem
emptyOdeProblem = OdeProblem
      { odeRhs = error "emptyOdeProblem: no odeRhs provided"
      , odeJacobian = Nothing
      , odeInitCond = error "emptyOdeProblem: no odeInitCond provided"
      , odeEvents = mempty
      , odeTimeBasedEvents = TimeEventSpec $ return $ 1.0 / 0.0
      , odeEventHandler = nilEventHandler
      , odeMaxEvents = 100
      , odeSolTimes = error "emptyOdeProblem: no odeSolTimes provided"
      , odeTolerances = defaultTolerances
      }

data OdeSolver = forall method . (Show method, Method method) => OdeSolver
  String -- name
  [method]

availableSolvers :: [OdeSolver]
availableSolvers =
  [ OdeSolver "CVode"  [BDF, ADAMS]
  , OdeSolver "ARKode" [SDIRK_5_3_4, TRBDF2_3_3_2]
  ]

defaultOpts :: method -> ODEOpts method
defaultOpts method = ODEOpts
  { maxNumSteps = 1e5
  , minStep     = 1.0e-14
  , fixedStep   = 0
  , maxFail     = 10
  , odeMethod   = method
  , initStep    = Nothing
  , jacobianRepr = DenseJacobian
  }

defaultTolerances = Tolerances
  { absTolerances = Left 1.0e-6
  , relTolerance = 1.0e-10
  }

checkDiscrepancy :: HasCallStack => Double -> Double -> Assertion
checkDiscrepancy eps diff = assertBool msg $ diff <= eps
  where
    msg = printf "Difference too large: %.2e > %.2e"
      diff eps

fmod :: RealFrac a => a -> a -> a
fmod arg1 arg2 =
  let
    quot = arg1 / arg2
    quot_floor = realToFrac (floor quot :: Integer)
    quot_floor_times_divisor = quot_floor * arg2
  in
    arg1 - quot_floor_times_divisor

-- | Show list without brackets, so that it's more suitable for file names
showList :: Show a => [a] -> String
showList = intercalate "," . map show

-- derive JSON instances for all sundials types for testing
deriving instance Generic SundialsSolution
instance ToJSON SundialsSolution
instance FromJSON SundialsSolution

deriving instance Generic SundialsDiagnostics
instance ToJSON SundialsDiagnostics
instance FromJSON SundialsDiagnostics

deriving instance Generic ErrorDiagnostics
instance ToJSON ErrorDiagnostics
instance FromJSON ErrorDiagnostics

instance ToJSON CrossingDirection
instance FromJSON CrossingDirection

instance ToJSON (Matrix Double) where
  toJSON = toJSON . toLists
instance FromJSON (Matrix Double) where
  parseJSON = fmap fromLists . parseJSON

compareSolutions
  :: Bool -- ^ same method?
  -> SundialsSolution -- ^ expected
  -> SundialsSolution -- ^ got
  -> Maybe String -- ^ maybe error message
compareSolutions same_method a b = asum @[]
  [ do
      guard . not $ VS.length (actualTimeGrid a) == VS.length (actualTimeGrid b)
      return "Different length of actualTimeGrid"
  , do
      guard . not $ norm_Inf (actualTimeGrid a - actualTimeGrid b) < precision
      return "Different values of actualTimeGrid"
  , do
      guard . not $ size (solutionMatrix a) == size (solutionMatrix b)
      return "Different sizes of the solutionMatrix"
  , do
      guard . not $ norm_Inf (solutionMatrix a - solutionMatrix b) < precision
      return "Different values in the solutionMatrix"
  ]
  where
    precision = if same_method then 1e-10 else 1e-1

odeGoldenTest
  :: forall method . (Method method, Show method)
  => Bool -- ^ compare between methods? (via a canonical golden file)
  -> ODEOpts method -- ^ ode options (affect the golden file name)
  -> String -- ^ name (of both the test and the file)
  -> IO (Either ErrorDiagnostics SundialsSolution)
  -> TestTree
odeGoldenTest do_canonical opts name action =
  let
    method_dir = show (methodSolver @method) </> show (odeMethod opts)
  in
    testGroup name $ do
      (dir, type_, same_method) <-
        [(method_dir, "Method-specific", True)] ++
        [("canonical", "Canonical", False) | do_canonical ]
      let
        golden_path = ("tests/golden" </> dir </> name <.> "json")
        cmp expected got = return $
          case (expected, got) of
            (Right{},Left{}) -> Just "expected success, got error"
            (Left{},Right{}) -> Just "expected error, got success"
            (Left a,Left b) ->
              if errorCode a == errorCode b
                then Nothing
                else Just "got a different error"
            (Right a,Right b) -> compareSolutions same_method a b
        upd val = createDirectoriesAndWriteFile golden_path
          (encodePretty' defConfig {confCompare = compare} val)
      return $ goldenTest
        type_ -- test name
        (fromJust <$> decodeFileStrict golden_path) -- get the golden value
        action -- get the tested value
        cmp
        upd

-- | Make an 'EventHandler' from a vector of per-event handlers and
-- parameters
mkEventHandler
  :: V.Vector (Double -> VS.Vector Double -> VS.Vector Double)
  -> V.Vector Bool -- ^ stop the solver?
  -> V.Vector Bool -- ^ record the event?
  -> EventHandler
mkEventHandler handlers stop_solver_vec record_event_vec t y0 evs
  | VS.null evs = error "mkEventHandler: got a time-based event"
  | otherwise = return EventHandlerResult
      { eventStopSolver = or . map (stop_solver_vec V.!) $ VS.toList evs
      , eventRecord = or . map (record_event_vec V.!) $ VS.toList evs
      , eventNewState = foldl' (\y hndl -> hndl t y) y0 . map (handlers V.!) $ VS.toList evs
      }

mkTimeEvents
  :: [(Double, Double -> VS.Vector Double -> VS.Vector Double, Bool, Bool)] -- ^ time-based events (time, update, stop?, record?)
  -> IO (TimeEventSpec, EventHandler)
mkTimeEvents time_based_events = do
  time_based_events_ref <- newIORef time_based_events
  let handler :: EventHandler
      handler t y0 evs =
        if VS.null evs
          then do
            time_evs <- readIORef time_based_events_ref
            case time_evs of
              [] -> throwIO $ ErrorCall "Unexpected time-based event"
              (t1, update, eventStopSolver, eventRecord) : rest ->
                if t == t1
                  then do
                    writeIORef time_based_events_ref rest
                    return EventHandlerResult
                      { eventStopSolver
                      , eventRecord
                      , eventNewState = update t y0
                      }
                  else throwIO $ ErrorCall "Wrong event time"
          else error "mkTimeEvents: got root-based events"

      event_spec :: TimeEventSpec
      event_spec = TimeEventSpec $ do
        time_evs <- readIORef time_based_events_ref
        return $ case time_evs of
          [] -> 1.0/0.0 -- +Inf
          (t,_,_,_) : _ -> t


  return (event_spec, handler)

combineEventHandlers
  :: EventHandler -- ^ root-based event handler
  -> EventHandler -- ^ time-based event handler
  -> EventHandler -- ^ combined handler
combineEventHandlers rh th t y0 evs =
  if VS.null evs
    then th t y0 evs
    else rh t y0 evs

nilEventHandler :: EventHandler
nilEventHandler _ _ _ = throwIO $ ErrorCall "nilEventHandler"

----------------------------------------------------------------------
--                             The tests
----------------------------------------------------------------------

main = do
  handleScribe <- mkHandleScribe ColorIfTerminal stderr (permitItem InfoS) V2
  log_env <- registerScribe "stderr" handleScribe defaultScribeSettings =<<
    initLogEnv "test" "devel"
  let ?log_env = log_env

  defaultMain $ testGroup "Tests" $
    [
      testGroup solver_name
      [ testGroup (show method) $
          let opts = defaultOpts method in
          [ withVsWithoutJacobian opts
          , stiffishTest opts
          , eventTests opts
          , noErrorTests opts
          , discontinuousRhsTest opts
          , modulusEventTest opts
          , cascadingEventsTest opts
          , simultaneousEventsTest opts
          , timeBasedEventTest opts
          , timeGridTest opts
          , timeAliasingTests opts
          ]
      | method <- methods
      ]
    | OdeSolver solver_name methods <- availableSolvers
    ]

noErrorTests opts = testGroup "Absence of error"
  [ testCase name $ do
      r <- runKatipT ?log_env $ solve opts prob
      case r of
        Right _ -> return ()
        Left e -> assertFailure (show e)
  | (name, prob) <- [ empty ]
  ]

stiffishTest opts = odeGoldenTest True opts "Stiffish" $
  runKatipT ?log_env $ solve opts stiffish

withVsWithoutJacobian opts0 = testGroup "With vs without jacobian"
  [ testCase (name ++ " " ++ (if sparse then "sparse" else "dense")) $ do
      let
        opts = opts0
          { jacobianRepr =
              if sparse
                then SparseJacobian 40 -- this should be enough for our tests
                else DenseJacobian
          }
      Right (solutionMatrix -> solJac)   <- runKatipT ?log_env $ solve opts prob
      Right (solutionMatrix -> solNoJac) <- runKatipT ?log_env $ solve
        opts { jacobianRepr = DenseJacobian }
        prob { odeJacobian = Nothing }
      checkDiscrepancy 1e-2 $ norm_2 (solJac - solNoJac)
  | (name, prob) <- [ brusselator, robertson ]
  , sparse <- [False, True]
  ]

eventTests opts = testGroup "Events"
  [ odeGoldenTest True opts "Exponential" $
      runKatipT ?log_env $ solve opts exponential
  , odeGoldenTest True opts "Robertson" $ do
      let upd _ _ = vector [1.0, 0.0, 0.0]
      runKatipT ?log_env $ solve opts
        (snd robertson)
          { odeEvents =
            [ EventSpec { eventCondition = \_t y -> y ! 0 - 0.0001
                        , eventDirection = AnyDirection
                        }
            , EventSpec { eventCondition = \_t y -> y ! 2 - 0.01
                        , eventDirection = AnyDirection
                        }
            ]
          , odeMaxEvents = 100
          , odeSolTimes = [0,100]
          , odeEventHandler = mkEventHandler
            (V.replicate 2 upd)
            (V.replicate 2 False)
            (V.replicate 2 True)
          }
  , odeGoldenTest True opts "Bounded_sine" $
      runKatipT ?log_env $ solve opts boundedSine
  ]

discontinuousRhsTest opts = odeGoldenTest True opts "Discontinuous_derivative" $
  runKatipT ?log_env $ solve opts discontinuousRHS

modulusEventTest opts0 = localOption (mkTimeout 1e5) $ testGroup "Modulus event"
  [ odeGoldenTest False opts
      (printf "Modulus_recordEvent=%s_initStep=%s"
        (show record_event)
        initStepStr) $
        runKatipT ?log_env $ solve opts (modulusEvent record_event)
        -- the timeout doesn't seem to trigger when the infinite
        -- loop actually happens. I'm not sure why, since we're making safe
        -- calls, which should be interruptible? -fno-omit-yields doesn't
        -- seem to help either. Maybe worth investigating.
  | (initStep, initStepStr::String) <-
      [(Nothing, "Nothing")
      ,(Just 1, "1")
      ,(Just (1 - 2**(-53)), "1-eps")]
  , record_event <- [False, True]
  , let opts = opts0 { initStep }
  ]

cascadingEventsTest opts = odeGoldenTest True opts "Cascading_events" $ do
  runKatipT ?log_env $ solve opts prob
  where
    prob = emptyOdeProblem
      { odeRhs = odeRhsPure $ \_ _ -> [1, 0]
      , odeInitCond = [0, 0]
      , odeEvents = events
      , odeEventHandler = mkEventHandler
          [\_ y -> [7, y ! 1]
          ,\_ y -> [y ! 0, 1]
          ]
          (V.replicate 2 False)
          (V.replicate 2 True)
      , odeMaxEvents = 100
      , odeSolTimes = [0,10]
      }
    events =
      [ EventSpec { eventCondition = \_t y -> y ! 0 - 5
                  , eventDirection = AnyDirection
                  }
      , EventSpec { eventCondition = \_t y -> y ! 0 - 6
                  , eventDirection = AnyDirection
                  }
      ]

simultaneousEventsTest opts = testGroup "Simultaneous events"
  [ odeGoldenTest True opts (printf "Simultaneous_events/maxEvents=%d_eventRecord=%s_stopSolver=%s"
    maxEvents
    (showList record)
    (showList stopSolver)) $ do
      let prob = emptyOdeProblem
            { odeRhs = odeRhsPure $ \_ _ -> [0,0]
            , odeJacobian = Nothing
            , odeInitCond = [0,0]
            , odeEvents = events
            , odeEventHandler = mkEventHandler
                (V.map (\i _ y -> y VS.// [(i,1)]) [0..1])
                (V.fromList stopSolver)
                (V.fromList record)
            , odeMaxEvents = maxEvents
            , odeSolTimes = [0,10]
            , odeTolerances = defaultTolerances
            }
          event =
            EventSpec
              { eventCondition = \t _ -> t - 5
              , eventDirection = AnyDirection
              }
          events = V.replicate 2 event
      runKatipT ?log_env $ solve opts prob
  | maxEvents <- [0..3]
  , record <- replicateM 2 [False,True]
  , stopSolver <- replicateM 2 [False,True]
  ]

-- Four events:
--
-- 1. A time-based event at t = 1.0; recorded
--
-- 2. A time-based event at t = 2.0; not recorded and no update
--
-- 3. A root-based event at t = 3.0
--
-- 4. A time-based event at t = 4.0; stops the solver
timeBasedEventTest opts = odeGoldenTest True opts "Time-based events" $ do
  let
    upd :: Double -> Double -> Vector Double -> Vector Double
    upd x _t y = y VS.// [(1, x)]
  (time_ev_spec, time_ev_handler) <- mkTimeEvents
    [ (1.0, upd 5, False, True)
    , (2.0, const id, False, False)
    , (4.0, upd 8, True, True)
    ]
  let prob = emptyOdeProblem
        { odeRhs = odeRhsPure $ \_ _ -> [1, 0]
        , odeInitCond = [0, 0]
        , odeEvents =
            [ EventSpec
              { eventCondition = \t y -> t/2 + y ! 0 - 4.5
              , eventDirection = Upwards
              }
            ]
        , odeTimeBasedEvents = time_ev_spec
        , odeEventHandler = combineEventHandlers
            (mkEventHandler [upd 13] [False] [True])
            time_ev_handler
        , odeMaxEvents = 100
        , odeSolTimes = [0, 10]
        }
  runKatipT ?log_env $ solve opts prob

-- | A simple test that checks that we do the right thing even when the
-- time grid does not start at 0.
timeGridTest opts = odeGoldenTest True opts "Time grid test" $ do
  runKatipT ?log_env $ solve opts emptyOdeProblem
    { odeRhs = odeRhsPure $ \_ _ -> [1]
    , odeInitCond = [0]
    , odeSolTimes = [-3, -2, 0, 10, 100]
    }

-- | Test time aliasing: when some of the time-based events coincide
-- exactly with time grid points.
--
-- There are two tests: for when events are recorded or not. This
-- corresponds to time-based events with the conditions that either hold or
-- not.
timeAliasingTests opts = testGroup "Time aliasing"
  [ odeGoldenTest True opts ("Time aliasing; record="++show record) $ do
      (time_ev_spec, time_ev_handler) <- mkTimeEvents
        [ (t, if record then \_ y -> [y!0 + 0.2] else const id, False, record)
        | t <- VS.toList ts
        ]
      runKatipT ?log_env $ solve opts emptyOdeProblem
        { odeRhs = odeRhsPure $ \_ _ -> [1]
        , odeInitCond = [0]
        , odeSolTimes = ts
        , odeTimeBasedEvents = time_ev_spec
        , odeEventHandler = time_ev_handler
        }
  | record <- [False, True]
  ]
  where
    ts = [1,2,3]

----------------------------------------------------------------------
--                           ODE problems
----------------------------------------------------------------------

brusselator :: (String, OdeProblem)
brusselator = (,) "brusselator" $ emptyOdeProblem
  { odeRhs = odeRhsPure $ \_t x ->
      let
        u = x VS.! 0
        v = x VS.! 1
        w = x VS.! 2
      in
      [ a - (w + 1) * u + v * u * u
      , w * u - v * u * u
      , (b - w) / eps - w * u
      ]
  , odeJacobian = Just $ \(_t :: Double) x ->
      let
        u = x VS.! 0
        v = x VS.! 1
        w = x VS.! 2
      in (3><3)
      [ (-(w + 1.0)) + 2.0 * u * v, w - 2.0 * u * v, (-w)
      , u * u                     , (-(u * u))     , 0.0
      , (-u)                      , u              , (-1.0) / eps - u
      ]
  , odeEvents = mempty
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = [1.2, 3.1, 3.0]
  , odeSolTimes = [0.0, 0.1 .. 10.0]
  , odeTolerances = defaultTolerances
  }
  where
    a = 1.0
    b = 3.5
    eps :: Fractional a => a
    eps = 5.0e-6

exponential = emptyOdeProblem
  { odeRhs = odeRhsPure $ \_ y -> [y VS.! 0]
  , odeJacobian = Nothing
  , odeInitCond = vector [1]
  , odeEvents = events
  , odeEventHandler = mkEventHandler [\_ _ -> [ 2 ]] [False] [True]
  , odeMaxEvents = 100
  , odeSolTimes = vector [ fromIntegral k / 100 | k <- [0..(22::Int)]]
  , odeTolerances = defaultTolerances
  }
  where
    events =
      [ EventSpec { eventCondition = \_ y -> y ! 0 - 1.1
                  , eventDirection = Upwards
                  }
      ]

robertson = (,) "Robertson" $ emptyOdeProblem
  { odeRhs = odeRhsPure $ \_ (VS.toList -> [y1,y2,y3]) ->
      [ -0.04 * y1 + 1.0e4 * y2 * y3
      , 0.04 * y1 - 1.0e4 * y2 * y3 - 3.0e7 * (y2)^(2 :: Int)
      , 3.0e7 * (y2)^(2 :: Int)
      ]
  , odeJacobian = Just $ \_t (VS.toList -> [_, y2, y3]) -> (3 >< 3)
      [ -0.04, 1.0e4 * y3, 1.0e4 * y2
      , 0.04, -1.0e4*y3 - 3.0e7*2*y2, -1.0e4*y2
      , 0, 3.0e7*2*y2, 0
      ]
  , odeInitCond = [1.0, 0.0, 0.0]
  , odeEvents = []
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeSolTimes = [0,20,100,1000]
  , odeTolerances = defaultTolerances
  }

empty = (,) "Empty system" $ emptyOdeProblem
  { odeRhs = odeRhsPure $ \_ _ -> []
  , odeJacobian = Nothing
  , odeInitCond = []
  , odeEvents = []
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeSolTimes = [0,1]
  , odeTolerances = defaultTolerances
  }

stiffish = emptyOdeProblem
  { odeRhs = odeRhsPure $ \t ((VS.! 0) -> u) -> [ lamda * u + 1.0 / (1.0 + t * t) - lamda * atan t ]
  , odeJacobian = Nothing
  , odeInitCond = [0.0]
  , odeEvents = []
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeSolTimes = [0.0, 0.1 .. 10.0]
  , odeTolerances = defaultTolerances
  }
  where
    lamda = -100.0

-- A sine wave that only changes direction once it reaches ±0.9.
-- Illustrates event-specific reset function
boundedSine = emptyOdeProblem
  { odeRhs = odeRhsPure $ \_t y -> [y VS.! 1, - y VS.! 0]
  , odeJacobian = Nothing
  , odeInitCond = [0,1]
  , odeEvents = events
  , odeEventHandler = mkEventHandler
      [\_ y -> vector [ y ! 0, - abs (y ! 1) ]
      ,\_ y -> vector [ y ! 0, abs (y ! 1) ]
      ]
      (V.replicate 2 False)
      (V.replicate 2 True)
  , odeMaxEvents = 100
  , odeSolTimes = VS.fromList [ 2 * pi * k / 360 | k <- [0..360]]
  , odeTolerances = defaultTolerances
  }
  where
    events =
      [ EventSpec { eventCondition = \_t y -> y ! 0 - 0.9
                  , eventDirection = Upwards
                  }
      , EventSpec { eventCondition = \_t y -> y ! 0 + 0.9
                  , eventDirection = Downwards
                  }
      ]

-- | An example of a system with a discontinuous RHS
discontinuousRHS = emptyOdeProblem
  { odeRhs = odeRhsPure $ \t _ ->
      if t1 <= t && t <= t2
        then [deriv]
        else [0]
  , odeJacobian = Nothing
  , odeInitCond = [0]
  , odeEvents = events
  , odeEventHandler = mkEventHandler
      (V.replicate 2 (\_ y -> y))
      (V.replicate 2 False)
      (V.replicate 2 False)
  , odeMaxEvents = 10
  , odeSolTimes = [0,1]
  , odeTolerances = defaultTolerances
  }
  where
    t1, t2 :: Fractional a => a
    t1 = 0.01
    t2 = 0.02
    deriv = 100
    events =
      [ EventSpec
        { eventCondition = \t _ -> t - t1
        , eventDirection = Upwards
        }
      , EventSpec
        { eventCondition = \t _ -> t - t2
        , eventDirection = Upwards
        }
      ]

modulusEvent record_event = emptyOdeProblem
  { odeRhs = odeRhsPure $ \t _ -> [t `fmod` 1]
  , odeJacobian = Nothing
  , odeInitCond = [0]
  , odeEvents = events
  , odeEventHandler = mkEventHandler
      [\_ y -> y]
      [False]
      [record_event]
  , odeMaxEvents = 11
  , odeSolTimes = [0,10]
  , odeTolerances = defaultTolerances
  }
  where
    events =
      [ EventSpec
        { eventCondition = \t _ ->
            let
              a = t
              b = 1
            in
              abs((a-b/2) `fmod` (2*b) - b) - b/2
        , eventDirection = AnyDirection
        }
      ]
