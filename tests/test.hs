{-# LANGUAGE RecordWildCards, ScopedTypeVariables, OverloadedStrings,
             ViewPatterns, ImplicitParams, OverloadedLists, RankNTypes,
             ExistentialQuantification, LambdaCase, NumDecimals, NamedFieldPuns #-}

import Prelude hiding (quot)

import Test.Tasty
import Test.Tasty.HUnit

import Numeric.Sundials

import Numeric.LinearAlgebra as L
import qualified Data.Vector.Storable as VS
import qualified Data.Vector as V
import Katip
import System.IO
import Text.Printf (printf)

----------------------------------------------------------------------
--                            Helpers
----------------------------------------------------------------------

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
  , maxFail     = 10
  , odeMethod   = method
  , initStep    = Nothing
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

-- | Enumerate all distinct unordered pairs of distinct elements
allPairs :: [a] -> [(a,a)]
allPairs = \case
  [] -> []
  x : xs -> map ((,) x) xs ++ allPairs xs

fmod :: RealFrac a => a -> a -> a
fmod arg1 arg2 =
  let
    quot = arg1 / arg2
    quot_floor = realToFrac (floor quot :: Integer)
    quot_floor_times_divisor = quot_floor * arg2
  in
    arg1 - quot_floor_times_divisor

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
          , eventTests opts
          , noErrorTests opts
          , discontinuousRhsTest opts
          , modulusEventTest opts
          ]
      | method <- methods
      ]
    | OdeSolver solver_name methods <- availableSolvers
    ] ++
    [ testGroup "Method comparison"
      -- FIXME rewrite this to be O(n) instead of O(n^2)
      -- right now, this only compares between different solvers
      [ testGroup (show method1 ++ " vs " ++ show method2) $ compareMethodsTests
          (defaultOpts method1)
          (defaultOpts method2)
      | (OdeSolver _ methods1, OdeSolver _ methods2) <- allPairs availableSolvers
      , method1 <- methods1
      , method2 <- methods2
      ]
    ]

noErrorTests opts = testGroup "Absence of error"
  [ testCase name $ do
      r <- runKatipT ?log_env $ solve opts prob
      case r of
        Right _ -> return ()
        Left e -> assertFailure (show e)
  | (name, prob) <- [ empty ]
  ]

withVsWithoutJacobian opts = testGroup "With vs without jacobian"
  [ testCase name $ do
      Right (solutionMatrix -> solJac)   <- runKatipT ?log_env $ solve opts prob
      Right (solutionMatrix -> solNoJac) <- runKatipT ?log_env $ solve opts prob { odeJacobian = Nothing }
      checkDiscrepancy 1e-2 $ norm_2 (solJac - solNoJac)
  | (name, prob) <- [ brusselator, robertson ]
  ]

compareMethodsTests opts1 opts2 =
  [ testCase name $ do
      Right (solutionMatrix -> sol1) <- runKatipT ?log_env $ solve opts1 prob
      Right (solutionMatrix -> sol2) <- runKatipT ?log_env $ solve opts2 prob
      let diff = maximum $ map abs $
                 zipWith (-) ((toLists $ tr sol1)!!0) ((toLists $ tr sol2)!!0)
      checkDiscrepancy 1e-5 diff
  | (name, prob) <- [ stiffish ]
  ]


eventTests opts = testGroup "Events"
  [ testCase "Exponential" $ do
      Right (eventInfo -> events) <- runKatipT ?log_env $ solve opts exponential
      length events @?= 1
      checkDiscrepancy 1e-4 (abs (eventTime (events V.! 0) - log 1.1))
      rootDirection (events V.! 0) @?= Upwards
      eventIndex (events V.! 0) @?= 0
  , testCase "Robertson" $ do
      let upd _ _ = vector [1.0, 0.0, 0.0]
      Right (eventInfo -> events) <- runKatipT ?log_env $ solve opts
        (snd robertson)
          { odeEvents =
            [ EventSpec { eventCondition = \_t y -> y ! 0 - 0.0001
                        , eventUpdate = upd
                        , eventDirection = AnyDirection
                        , eventStopSolver = False
                        , eventRecord = True
                        }
            , EventSpec { eventCondition = \_t y -> y ! 2 - 0.01
                        , eventUpdate = upd
                        , eventDirection = AnyDirection
                        , eventStopSolver = False
                        , eventRecord = True
                        }
            ]
          , odeMaxEvents = 100
          , odeSolTimes = [0,100]
          }
      length events @?= 100
  , testCase "Bounded sine" $ do
      Right (eventInfo -> events) <- runKatipT ?log_env $ solve opts boundedSine
      length events @?= 3
      V.map rootDirection events @?= [Upwards, Downwards, Upwards]
      V.map eventIndex events @?= [0, 1, 0]
      V.forM_ (V.zip (V.map eventTime events) [1.119766,3.359295,5.598820]) $ \(et_got, et_exp) ->
        checkDiscrepancy 1e-4 (abs (et_exp - et_got))
  ]

discontinuousRhsTest opts = testCaseInfo "Discontinuous derivative" $ do
  Right r <- runKatipT ?log_env $ solve opts discontinuousRHS
  V.length (eventInfo r) @?= 0
  let mx = solutionMatrix r
  rows mx @?= 2 -- because the auxiliary events are not recorded
  let y1 = mx ! 1 ! 0
      diff = abs (y1 - 1)
  checkDiscrepancy 0.1 (abs (y1 - 1))
  return $ printf "%.2e" diff

modulusEventTest opts0 = localOption (mkTimeout 1e5) $ testGroup "Modulus event"
  [ testGroup (if record_event then "Recording events" else "Not recording events")
    [ testCaseInfo ("Init step is " ++ show initStep) $ do
        let opts = opts0 { initStep }
        Right r <- runKatipT ?log_env $ solve opts (modulusEvent record_event)
        V.length (eventInfo r) @?=
          (if record_event then 10 else 0)
        let mx = solutionMatrix r
        rows mx @?=
          (if record_event then 22 else 2) -- because the auxiliary events are not recorded
        let y1 = mx ! 1 ! 0
        -- We don't check the answer here; all we care about is not failing
        -- or entering an infinite loop (hence the timeout).
        -- (Side note: the timeout doesn't seem to trigger when the infinite
        -- loop actually happens. I'm not sure why, since we're making safe
        -- calls, which should be interruptible? -fno-omit-yields doesn't
        -- seem to help either. Maybe worth investigating.)
        -- When step size is big, the result will not be accurate.
        -- However, we display it just FYI:
        return $ printf "Result: %.3f (expected 5.0)" y1
    | initStep <- [Nothing, Just 1, Just (1 - 2**(-53))]
    ]
  | record_event <- [False, True]
  ]

----------------------------------------------------------------------
--                           ODE problems
----------------------------------------------------------------------

brusselator :: (String, OdeProblem)
brusselator = (,) "brusselator" $ OdeProblem
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

exponential = OdeProblem
  { odeRhs = odeRhsPure $ \_ y -> [y VS.! 0]
  , odeJacobian = Nothing
  , odeInitCond = vector [1]
  , odeEvents = events
  , odeMaxEvents = 100
  , odeSolTimes = vector [ fromIntegral k / 100 | k <- [0..(22::Int)]]
  , odeTolerances = defaultTolerances
  }
  where
    events =
      [ EventSpec { eventCondition = \_ y -> y ! 0 - 1.1
                  , eventUpdate = \_ _ -> vector [ 2 ]
                  , eventDirection = Upwards
                  , eventStopSolver = False
                  , eventRecord = True
                  }
      ]

robertson = (,) "Robertson" $ OdeProblem
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
  , odeMaxEvents = 0
  , odeSolTimes = [0,20]
  , odeTolerances = defaultTolerances -- FIXME how to make this integrate indefinitely, as in the sundials example?
  }

empty = (,) "Empty system" $ OdeProblem
  { odeRhs = odeRhsPure $ \_ _ -> []
  , odeJacobian = Nothing
  , odeInitCond = []
  , odeEvents = []
  , odeMaxEvents = 0
  , odeSolTimes = [0,1]
  , odeTolerances = defaultTolerances
  }

stiffish = (,) "Stiffish" $ OdeProblem
  { odeRhs = odeRhsPure $ \t ((VS.! 0) -> u) -> [ lamda * u + 1.0 / (1.0 + t * t) - lamda * atan t ]
  , odeJacobian = Nothing
  , odeInitCond = [0.0]
  , odeEvents = []
  , odeMaxEvents = 0
  , odeSolTimes = [0.0, 0.1 .. 10.0]
  , odeTolerances = defaultTolerances
  }
  where
    lamda = -100.0

-- A sine wave that only changes direction once it reaches ±0.9.
-- Illustrates event-specific reset function
boundedSine = OdeProblem
  { odeRhs = odeRhsPure $ \_t y -> [y VS.! 1, - y VS.! 0]
  , odeJacobian = Nothing
  , odeInitCond = [0,1]
  , odeEvents = events
  , odeMaxEvents = 100
  , odeSolTimes = VS.fromList [ 2 * pi * k / 360 | k <- [0..360]]
  , odeTolerances = defaultTolerances
  }
  where
    events =
      [ EventSpec { eventCondition = \_t y -> y ! 0 - 0.9
                     , eventUpdate = \_ y -> vector [ y ! 0, - abs (y ! 1) ]
                     , eventDirection = Upwards
                     , eventStopSolver = False
                     , eventRecord = True
                     }
      , EventSpec { eventCondition = \_t y -> y ! 0 + 0.9
                     , eventUpdate = \_ y -> vector [ y ! 0, abs (y ! 1) ]
                     , eventDirection = Downwards
                     , eventStopSolver = False
                     , eventRecord = True
                     }
      ]

-- | An example of a system with a discontinuous RHS
discontinuousRHS = OdeProblem
  { odeRhs = odeRhsPure $ \t _ ->
      if t1 <= t && t <= t2
        then [deriv]
        else [0]
  , odeJacobian = Nothing
  , odeInitCond = [0]
  , odeEvents = events
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
        , eventUpdate = \_ y -> y
        , eventDirection = Upwards
        , eventStopSolver = False
        , eventRecord = False
        }
      , EventSpec
        { eventCondition = \t _ -> t - t2
        , eventUpdate = \_ y -> y
        , eventDirection = Upwards
        , eventStopSolver = False
        , eventRecord = False
        }
      ]

modulusEvent record_event = OdeProblem
  { odeRhs = odeRhsPure $ \t _ -> [t `fmod` 1]
  , odeJacobian = Nothing
  , odeInitCond = [0]
  , odeEvents = events
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
        , eventUpdate = \_ y -> y
        , eventDirection = AnyDirection
        , eventStopSolver = False
        , eventRecord = record_event
        }
      ]
