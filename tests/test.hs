{-# LANGUAGE RecordWildCards, ScopedTypeVariables, OverloadedStrings,
             ViewPatterns, ImplicitParams, OverloadedLists, RankNTypes,
             ExistentialQuantification, LambdaCase, NumDecimals #-}
import Test.Tasty
import Test.Tasty.HUnit

import Numeric.Sundials

import Numeric.LinearAlgebra as L
import qualified Data.Vector.Storable as VS
import qualified Data.Vector as V
import Katip
import Foreign.C.Types
import System.IO
import Text.Printf (printf)
import GHC.Stack
import Control.Monad
import Data.Bifunctor
import Data.Coerce

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
                        }
            , EventSpec { eventCondition = \_t y -> y ! 2 - 0.01
                        , eventUpdate = upd
                        , eventDirection = AnyDirection
                        , eventStopSolver = False
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
  Right (solutionMatrix -> mx) <- runKatipT ?log_env $ solve opts discontinuousRHS
  let y1 = VS.last $ flatten mx
      diff = abs (y1 - 1)
  checkDiscrepancy 0.1 (abs (y1 - 1))
  return $ printf "%.2e" diff

----------------------------------------------------------------------
--                           ODE problems
----------------------------------------------------------------------

brusselator :: (String, OdeProblem)
brusselator = (,) "brusselator" $ OdeProblem
  { odeRhs = OdeRhsHaskell $ \_t x ->
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
  { odeRhs = OdeRhsHaskell $ \_ y -> [y VS.! 0]
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
                  }
      ]

robertson = (,) "Robertson" $ OdeProblem
  { odeRhs = OdeRhsHaskell $ \_ (VS.toList -> [y1,y2,y3]) ->
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
  { odeRhs = OdeRhsHaskell $ \_ _ -> []
  , odeJacobian = Nothing
  , odeInitCond = []
  , odeEvents = []
  , odeMaxEvents = 0
  , odeSolTimes = [0,1]
  , odeTolerances = defaultTolerances
  }

stiffish = (,) "Stiffish" $ OdeProblem
  { odeRhs = OdeRhsHaskell $ \t ((VS.! 0) -> u) -> [ lamda * u + 1.0 / (1.0 + t * t) - lamda * atan t ]
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
  { odeRhs = OdeRhsHaskell $ \_t y -> [y VS.! 1, - y VS.! 0]
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
                     }
      , EventSpec { eventCondition = \_t y -> y ! 0 + 0.9
                     , eventUpdate = \_ y -> vector [ y ! 0, abs (y ! 1) ]
                     , eventDirection = Downwards
                     , eventStopSolver = False
                     }
      ]

-- | An example of a system with a discontinuous RHS
discontinuousRHS = OdeProblem
  { odeRhs = OdeRhsHaskell $ \t _ ->
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
        , eventUpdate = \t y -> y
        , eventDirection = Upwards
        , eventStopSolver = False
        }
      , EventSpec
        { eventCondition = \t _ -> t - t2
        , eventUpdate = \t y -> y
        , eventDirection = Upwards
        , eventStopSolver = False
        }
      ]

largeTs :: V.Vector Double
largeTs = V.fromList $ 0.0 : take 12 (iterate (*10) 0.04)
