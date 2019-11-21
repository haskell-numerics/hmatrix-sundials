{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -Wall #-}

import qualified Numeric.Sundials.ARKode.ODE as ARK
import qualified Numeric.Sundials.CVode.ODE  as CV
import           Numeric.LinearAlgebra as L
import           Numeric.Sundials.Types

import           Plots as P
import qualified Diagrams.Prelude as D
import           Diagrams.Backend.Rasterific
import qualified Data.Vector.Storable as V

import           Control.Lens
import           Control.Monad
import           Data.Coerce
import           Foreign.C.Types

import           Test.Hspec

import Control.Monad.IO.Class
import System.IO
import Katip


lorenz :: Double -> [Double] -> [Double]
lorenz _t u = [ sigma * (y - x)
              , x * (rho - z) - y
              , x * y - beta * z
              ]
  where
    rho = 28.0
    sigma = 10.0
    beta = 8.0 / 3.0
    x = u !! 0
    y = u !! 1
    z = u !! 2

_lorenzJac :: Double -> Vector Double -> Matrix Double
_lorenzJac _t u = (3><3) [ (-sigma), rho - z, y
                        , sigma   , -1.0   , x
                        , 0.0     , (-x)   , (-beta)
                        ]
  where
    rho = 28.0
    sigma = 10.0
    beta = 8.0 / 3.0
    x = u ! 0
    y = u ! 1
    z = u ! 2

brusselator :: Double -> [Double] -> [Double]
brusselator _t x = [ a - (w + 1) * u + v * u * u
                   , w * u - v * u * u
                   , (b - w) / eps - w * u
                   ]
  where
    a = 1.0
    b = 3.5
    eps = 5.0e-6
    u = x !! 0
    v = x !! 1
    w = x !! 2

brussJac :: Double -> Vector Double -> Matrix Double
brussJac _t x = tr $
  (3><3) [ (-(w + 1.0)) + 2.0 * u * v, w - 2.0 * u * v, (-w)
         , u * u                     , (-(u * u))     , 0.0
         , (-u)                      , u              , (-1.0) / eps - u
         ]
  where
    y = toList x
    u = y !! 0
    v = y !! 1
    w = y !! 2
    eps = 5.0e-6

brusselatorWithJacobian :: (MonadIO m, Katip m) => Vector Double -> Bool -> m CV.SolverResult
brusselatorWithJacobian ts usejac = CV.odeSolveRootVWith' opts
                      (OdeRhsHaskell . coerce $ \t v -> vector $ brusselator t (toList v))
                      (if usejac then Just brussJac else Nothing)
                      (vector [1.2, 3.1, 3.0])
                      [] 0
                      ts
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod = CV.BDF
                   , stepControl = CV.XX' 1.0e-6 1.0e-10 1 1
                   , initStep = Nothing
                   }

brussRoot :: (MonadIO m, Katip m) => m CV.SolverResult
brussRoot = CV.odeSolveRootVWith' opts
                      (OdeRhsHaskell . coerce $ \t v -> vector $ brusselator t (toList v))
                      Nothing
                      (vector [1.2, 3.1, 3.0])
                      events 100
                      (vector [0.0, 0.1 .. 10.0])
  where
    events =
      [ EventSpec { eventCondition = brussRootFn
                     , eventUpdate =
                         \_ev x -> let y = toList x in vector [(y!!0) + 0.5 , (y!!1), (y!!2)]
                     , eventDirection = AnyDirection
                     , eventStopSolver = False
                     }
      ]
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod = CV.BDF
                   , stepControl = CV.XX' 1.0e-6 1.0e-10 1 1
                   , initStep = Nothing
                   }

brussRootFn :: Double -> Vector Double -> Double
brussRootFn _ v = case xs of
                    [y1, _y2, y3] -> y1 - y3
                    _            -> error "brusselator root function RHS not defined"
  where
    xs = toList v

exponential :: (MonadIO m, Katip m) => m CV.SolverResult
exponential = CV.odeSolveRootVWith' opts
                      (OdeRhsHaskell . coerce $ \(t :: Double) y -> vector [y ! 0])
                      Nothing
                      (vector [1])
                      events 100
                      (vector [ fromIntegral k / 100 | k <- [0..(22::Int)]])
  where
    events =
      [ EventSpec { eventCondition = \t y -> y ! 0 - 1.1
                     , eventUpdate = \ev y -> vector [ 2 ]
                     , eventDirection = Upwards
                     , eventStopSolver = False
                     }
      ]
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod = CV.BDF
                   , stepControl = CV.XX' 1.0e-6 1.0e-10 1 1
                   , initStep = Nothing
                   }

-- A sine wave that only changes direction once it reaches ±0.9.
-- Illustrates event-specific reset function
boundedSine :: (MonadIO m, Katip m) => m CV.SolverResult
boundedSine = CV.odeSolveRootVWith'
  opts
  (OdeRhsHaskell . coerce $ \(_t :: Double) y -> vector [y ! 1, - y ! 0]) -- ODE RHS
  Nothing
  (vector [0, 1]) -- initial conditions
  events
  100 -- maximum number of events
  (vector [ 2 * pi * k / 360 | k <- [0..360]]) -- solution times
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod = CV.ADAMS
                   , stepControl = CV.XX' 1.0e-6 1.0e-10 1 1
                   , initStep = Nothing
                   }
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

stiffish :: Double -> [Double] -> [Double]
stiffish t v = [ lamda * u + 1.0 / (1.0 + t * t) - lamda * atan t ]
  where
    lamda = -100.0
    u = v !! 0

stiffishV :: Double -> Vector Double -> Vector Double
stiffishV t v = fromList [ lamda * u + 1.0 / (1.0 + t * t) - lamda * atan t ]
  where
    lamda = -100.0
    u = v ! 0

_stiffJac :: Double -> Vector Double -> Matrix Double
_stiffJac _t _v = (1><1) [ lamda ]
  where
    lamda = -100.0

predatorPrey :: Double -> [Double] -> [Double]
predatorPrey _t v = [ x * a - b * x * y
                    , d * x * y - c * y - e * y * z
                    , (-f) * z + g * y * z
                    ]
  where
    x = v!!0
    y = v!!1
    z = v!!2
    a = 1.0
    b = 1.0
    c = 1.0
    d = 1.0
    e = 1.0
    f = 1.0
    g = 1.0

roberts :: OdeRhs
roberts = OdeRhsHaskell . coerce $ \(t :: Double) v -> vector $ robertsAux t (toList v)
  where
    robertsAux _ [y1, y2, y3] =
      [ -0.04 * y1 + 1.0e4 * y2 * y3
      , 0.04 * y1 - 1.0e4 * y2 * y3 - 3.0e7 * (y2)^(2 :: Int)
      , 3.0e7 * (y2)^(2 :: Int)
      ]
    robertsAux _ _ = error "roberts RHS not defined"

robertsJac :: Double -> Vector Double -> Matrix Double
robertsJac _t (toList -> [y1, y2, y3]) = (3 >< 3)
  [ -0.04, 1.0e4 * y3, 1.0e4 * y2
  , 0.04, -1.0e4*y3 - 3.0e7*2*y2, -1.0e4*y2
  , 0, 3.0e7*2*y2, 0
  ]

ts :: [Double]
ts = take 12 $ map (* 10.0) (0.04 : ts)

solve :: (MonadIO m, Katip m) => m CV.SolverResult
solve = CV.odeSolveRootVWith' opts
                      roberts Nothing (vector [1.0, 0.0, 0.0])
                      events 100
                      (vector (0.0 : ts))
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod = CV.BDF
                   , stepControl = CV.ScXX' 1.0 1.0e-4 1.0 1.0 (vector [1.0e-8, 1.0e-14, 1.0e-6])
                   , initStep = Nothing
                   }
    events =
      [ EventSpec { eventCondition = \_t y -> y ! 0 - 0.0001
                     , eventUpdate = const id
                     , eventDirection = AnyDirection
                     , eventStopSolver = False
                     }
      , EventSpec { eventCondition = \_t y -> y ! 2 - 0.01
                     , eventUpdate = const id
                     , eventDirection = AnyDirection
                     , eventStopSolver = False
                     }
      ]

solve2 :: (MonadIO m, Katip m) => m CV.SolverResult
solve2 = CV.odeSolveRootVWith' opts
                      roberts Nothing (vector [1.0, 0.0, 0.0])
                      events 100
                      (vector (0.0 : ts))
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod = CV.BDF
                   , stepControl = CV.ScXX' 1.0 1.0e-4 1.0 1.0 (vector [1.0e-8, 1.0e-14, 1.0e-6])
                   , initStep = Nothing
                   }
    events =
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
    upd _ _ = vector [1.0, 0.0, 0.0]

solve1 :: (MonadIO m, Katip m) => m CV.SolverResult
solve1 = CV.odeSolveRootVWith' opts
                      roberts Nothing (vector [1.0, 0.0, 0.0])
                      events 100
                      (vector (0.0 : ts))
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod = CV.BDF
                   , stepControl = CV.ScXX' 1.0 1.0e-4 1.0 1.0 (vector [1.0e-8, 1.0e-14, 1.0e-6])
                   , initStep = Nothing
                   }
    events =
      [ EventSpec { eventCondition = \t _y -> t - 1.0
                     , eventUpdate = \t y -> vector [2.0, y!1, y!2]
                     , eventDirection = AnyDirection
                     , eventStopSolver = False
                     }
      ]

robertsonWithJacobian :: (MonadIO m, Katip m) => Vector Double -> Bool -> m CV.SolverResult
robertsonWithJacobian ts usejac = CV.odeSolveRootVWith' opts
                      roberts (if usejac then Just robertsJac else Nothing) (vector [1.0, 0.0, 0.0])
                      [] 0
                      ts
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod = CV.BDF
                   , stepControl = CV.ScXX' 1.0 1.0e-4 1.0 1.0 (vector [1.0e-8, 1.0e-14, 1.0e-6])
                   , initStep = Nothing
                   }

lSaxis :: [[Double]] -> P.Axis B D.V2 Double
lSaxis xs = P.r2Axis &~ do
  let zs = xs!!0
      us = xs!!1
      vs = xs!!2
      ws = xs!!3
  P.linePlot' $ zip zs us
  P.linePlot' $ zip zs vs
  P.linePlot' $ zip zs ws

lSaxis2 :: [[Double]] -> P.Axis B D.V2 Double
lSaxis2 xs = P.r2Axis &~ do
  let zs = xs!!0
      us = xs!!1
      vs = xs!!2
  P.linePlot' $ zip zs us
  P.linePlot' $ zip zs vs

kSaxis :: [(Double, Double)] -> P.Axis B D.V2 Double
kSaxis xs = P.r2Axis &~ do
  P.linePlot' xs


main :: IO ()
main = do
  handleScribe <- mkHandleScribe ColorIfTerminal stderr (permitItem InfoS) V2
  log_env <- registerScribe "stderr" handleScribe defaultScribeSettings =<< initLogEnv "test" "devel"

  runKatipT log_env $ do

  res1 <- ARK.odeSolve brusselator [1.2, 3.1, 3.0] (fromList [0.0, 0.1 .. 10.0])
  liftIO $ renderRasterific "diagrams/brusselator.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ lSaxis $ [0.0, 0.1 .. 10.0]:(toLists $ tr res1))

  res1a <- ARK.odeSolve brusselator [1.2, 3.1, 3.0] (fromList [0.0, 0.1 .. 10.0])
  liftIO $ renderRasterific "diagrams/brusselatorA.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ lSaxis $ [0.0, 0.1 .. 10.0]:(toLists $ tr res1a))

  res2 <- ARK.odeSolve stiffish [0.0] (fromList [0.0, 0.1 .. 10.0])
  liftIO $ renderRasterific "diagrams/stiffish.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip [0.0, 0.1 .. 10.0] (concat $ toLists res2))

  res2a <- ARK.odeSolveV (ARK.SDIRK_5_3_4') Nothing 1e-6 1e-10 stiffishV (fromList [0.0]) (fromList [0.0, 0.1 .. 10.0])

  res2b <- ARK.odeSolveV (ARK.TRBDF2_3_3_2') Nothing 1e-6 1e-10 stiffishV (fromList [0.0]) (fromList [0.0, 0.1 .. 10.0])

  let maxDiffA = maximum $ map abs $
                 zipWith (-) ((toLists $ tr res2a)!!0) ((toLists $ tr res2b)!!0)

  res2c <- CV.odeSolveV (CV.BDF) Nothing 1e-6 1e-10 stiffishV (fromList [0.0]) (fromList [0.0, 0.1 .. 10.0])

  let maxDiffB = maximum $ map abs $
                 zipWith (-) ((toLists $ tr res2a)!!0) ((toLists $ tr res2c)!!0)

  let maxDiffC = maximum $ map abs $
                 zipWith (-) ((toLists $ tr res2b)!!0) ((toLists $ tr res2c)!!0)

  res3 <- ARK.odeSolve lorenz [-5.0, -5.0, 1.0] (fromList [0.0, 0.01 .. 20.0])

  liftIO $ renderRasterific "diagrams/lorenz.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res3)!!0) ((toLists $ tr res3)!!1))

  liftIO $ renderRasterific "diagrams/lorenz1.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res3)!!0) ((toLists $ tr res3)!!2))

  liftIO $ renderRasterific "diagrams/lorenz2.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res3)!!1) ((toLists $ tr res3)!!2))

  res4 <- CV.odeSolve predatorPrey [0.5, 1.0, 2.0] (fromList [0.0, 0.01 .. 10.0])

  liftIO $ renderRasterific "diagrams/predatorPrey.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res4)!!0) ((toLists $ tr res4)!!1))

  liftIO $ renderRasterific "diagrams/predatorPrey1.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res4)!!0) ((toLists $ tr res4)!!2))

  liftIO $ renderRasterific "diagrams/predatorPrey2.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res4)!!1) ((toLists $ tr res4)!!2))

  res4a <- ARK.odeSolve predatorPrey [0.5, 1.0, 2.0] (fromList [0.0, 0.01 .. 10.0])

  let maxDiffPpA = maximum $ map abs $
                   zipWith (-) ((toLists $ tr res4)!!0) ((toLists $ tr res4a)!!0)

  let cond5 =
        runKatipT log_env solve >>= \case
          CV.SolverSuccess events _ _ -> do
            length events `shouldBe` 2
            (abs (eventTime (events!!0) - 0.2640208751331032) / 0.2640208751331032 < 1.0e-8) `shouldBe` True
            (abs (eventTime (events!!1) - 2.0786731062254436e7) / 2.0786731062254436e7 < 1.0e-8) `shouldBe` True
          CV.SolverError e ->
            error $ "Root finding error!\n" ++ show e

  let cond6 =
        runKatipT log_env solve1 >>= \case
          CV.SolverSuccess events _ _ -> do
            length events `shouldBe` 1
            (abs (eventTime (events!!0) - 1.0) / 1.0 < 1.0e-10) `shouldBe` True
          CV.SolverError _ ->
            error "Root finding error!"

  let cond7 =
        runKatipT log_env solve2 >>= \case
          CV.SolverSuccess events _ _ -> length events `shouldBe` 100
          CV.SolverError _ -> error "solver failed"
            True

  brussRoot >>= \case
    CV.SolverSuccess events m _diagn -> do
      liftIO $ renderRasterific
        "diagrams/brussRoot.png"
        (D.dims2D 500.0 500.0)
        (renderAxis $ lSaxis $ toLists $ tr m)
    CV.SolverError e ->
      liftIO $ expectationFailure $ show $ errorCode e

  let boundedSineSpec = do
        runKatipT log_env boundedSine >>= \case
          CV.SolverSuccess events m _ -> do
            liftIO $ renderRasterific
              "diagrams/boundedSine.png"
              (D.dims2D 500.0 500.0)
              (renderAxis $ lSaxis2 $ toLists $ tr m)
            length events `shouldBe` 3
            map rootDirection events `shouldBe` [Upwards, Downwards, Upwards]
            map eventIndex events `shouldBe` [0, 1, 0]
            forM_ (zip (map eventTime events) [1.1197660081724263,3.3592952656818404,5.5988203973243]) $ \(et_got, et_exp) ->
              et_got `shouldSatisfy` ((< 1e-8) . abs . subtract et_exp)
          CV.SolverError {} ->
            expectationFailure "Solver error"
  let exponentialSpec = do
        runKatipT log_env exponential >>= \case
          CV.SolverSuccess events _m _diagn -> do
            length events `shouldBe` 1
            (abs (eventTime (events!!0) - log 1.1) < 1e-4) `shouldBe` True
            rootDirection (events!!0) `shouldBe` Upwards
            eventIndex (events!!0) `shouldBe` 0
          CV.SolverError e ->
            expectationFailure $ show $ errorCode e

      robertsonJac = do
        let ts = vector [0, 1 .. 10]
        CV.SolverSuccess _ m1 _ <- runKatipT log_env $ robertsonWithJacobian ts True
        CV.SolverSuccess _ m2 _ <- runKatipT log_env $ robertsonWithJacobian ts False
        norm_2 (m1-m2) `shouldSatisfy` (< 1e-4)

      brusselatorJac = do
        let ts = [0.0, 0.1 .. 10.0]
        CV.SolverSuccess _ m1 _ <- runKatipT log_env $ brusselatorWithJacobian (vector ts) True
        CV.SolverSuccess _ m2 _ <- runKatipT log_env $ brusselatorWithJacobian (vector ts) False
        norm_2 (m1-m2) `shouldSatisfy` (< 1e-3)

  liftIO $ hspec $ do
    describe "Compare results" $ do
      it "Robertson should stop early" cond7
      it "Robertson time only" $ cond6
      it "Robertson from SUNDIALS manual" $ cond5
      it "Robertson with explicit Jacobian up to t=10" robertsonJac
      it "Brusselator with explicit Jacobian" brusselatorJac
      it "for SDIRK_5_3_4' and TRBDF2_3_3_2'" $ maxDiffA < 1.0e-6
      it "for SDIRK_5_3_4' and BDF" $ maxDiffB < 1.0e-6
      it "for TRBDF2_3_3_2' and BDF" $ maxDiffC < 1.0e-6
      it "for CV and ARK for the Predator Prey model" $ maxDiffPpA < 1.0e-3
    describe "Handling empty systems" $
      forM_ [("CVOde",CV.odeSolve),("ARKOde",ARK.odeSolve)] $ \(name, solveFn) ->
        it name $ do
          r <- runKatipT log_env $ solveFn (\_ _ -> []) [] (V.enumFromTo 0 10)
          r `shouldSatisfy` \sol -> L.size sol == (11,0)
    describe "Events" $ do
      it "Bounded sine events" $ boundedSineSpec
      it "Exponential events" $ exponentialSpec
      describe "Discontinuous zero crossings" $ do
        let
          eq :: OdeRhs
          eq = OdeRhsHaskell $ \_ _ -> V.singleton 1

          cond
            :: (Double -> Double -> Bool)
            -> (Double -> Vector Double -> Double)
          cond op _t y =
            if V.head y `op` 0
              then 1
              else -1

          solve op = CV.odeSolveWithEvents
            opts
            [EventSpec
              { eventCondition = cond op
              , eventDirection = AnyDirection
              , eventUpdate = \_t y -> V.map (+1) y
              , eventStopSolver = False
              }
            ]
            5 -- max # of events
            eq
            Nothing
            (V.singleton (-1))
            (V.fromList [0, 2])

          ops :: [(String, Double -> Double -> Bool)]
          ops =
            [ (">=", (>=))
            , (">",  (>))
            , ("<=", (<=))
            , ("<",  (<))
            ]
        forM_ ops $ \(op_name, op) -> it ("Event condition expressed as " ++ op_name) $ do
          Right soln <- runKatipT log_env $ solve op
          let
            evs = eventInfo soln
            [ev] = evs
          length evs `shouldBe` 1
          eventTime ev `shouldSatisfy` (\t -> abs (t-1) < 1e-3)
          -- row 0 is time 0, rows 1 and 2 are right before and right after
          -- the event (time 1), row 3 is the end point (time 2)
          solutionMatrix soln ! 1 ! 0 `shouldSatisfy` (\y -> abs y < 1e-3)
          solutionMatrix soln ! 2 ! 0 `shouldSatisfy` (\y -> abs (y-1) < 1e-3)
          solutionMatrix soln ! 3 ! 0 `shouldSatisfy` (\y -> abs (y-2) < 1e-3)
      it "Indicates when max_events is exceeded" $ do
        let solve max_events = CV.odeSolveWithEvents
              opts
              [EventSpec
                { eventCondition = \_ y -> y V.! 0 - 0.99999
                , eventDirection = Upwards
                , eventUpdate = \_t _y -> V.singleton 0
                , eventStopSolver = False
                }
              ]
              max_events
              (OdeRhsHaskell $ \_t _y -> V.singleton 1) -- the rhs
              Nothing
              (V.singleton 0) -- initial condition
              (V.fromList [0, 10]) -- solution times
        -- the event will be triggered 10 times on [0,10]
        Right soln1 <- runKatipT log_env $ solve 10
        length (eventInfo soln1) `shouldBe` 10
        odeMaxEventsReached (diagnostics soln1) `shouldBe` True
        Right soln2 <- runKatipT log_env $ solve 11
        length (eventInfo soln1) `shouldBe` 10
        odeMaxEventsReached (diagnostics soln2) `shouldBe` False
      it "Stops the solver when requested" $ do
        let solve max_events = CV.odeSolveWithEvents
              opts
              [EventSpec
                { eventCondition = \_ y -> y V.! 0 - 1
                , eventDirection = Upwards
                , eventUpdate = \_ y -> V.singleton 42
                , eventStopSolver = True
                }
              ]
              1 -- max_events
              (OdeRhsHaskell $ \_t _y -> V.singleton 1) -- the rhs
              Nothing
              (V.singleton 0) -- initial condition
              (V.fromList [0, 10]) -- solution times
        Right soln1 <- runKatipT log_env $ solve 10
        length (eventInfo soln1) `shouldBe` 1
        odeMaxEventsReached (diagnostics soln1) `shouldBe` False
        let ev = head $ eventInfo soln1
        eventTime ev `shouldSatisfy` (\x -> abs (x-1) < 1e-7)
        eventIndex ev `shouldBe` 0
        rootDirection ev `shouldBe` Upwards
        V.length (actualTimeGrid soln1) `shouldBe` 3
        actualTimeGrid soln1 V.! 1 `shouldSatisfy` (\x -> abs (x-1) < 1e-7)
        actualTimeGrid soln1 V.! 2 `shouldSatisfy` (\x -> abs (x-1) < 1e-7)
        L.size (solutionMatrix soln1) `shouldBe` (3,1)
        solutionMatrix soln1 ! 2 ! 0 `shouldBe` 42

  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod = CV.BDF
                   , stepControl = CV.ScXX' 1.0 1.0e-4 1.0 1.0 (vector [1.0e-8, 1.0e-14, 1.0e-6])
                   , initStep = Nothing
                   }
