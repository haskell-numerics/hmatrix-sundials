{-# OPTIONS_GHC -Wall #-}

import qualified Numeric.Sundials.ARKode.ODE as ARK
import qualified Numeric.Sundials.CVode.ODE  as CV
import           Numeric.LinearAlgebra
import           Numeric.Sundials.ODEOpts (ODEOpts(..))

import           Plots as P
import qualified Diagrams.Prelude as D
import           Diagrams.Backend.Rasterific

import           Control.Lens

import           Data.Functor.Compose

import           Test.Hspec


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

_brussJac :: Double -> Vector Double -> Matrix Double
_brussJac _t x = (3><3) [ (-(w + 1.0)) + 2.0 * u * v, w - 2.0 * u * v, (-w)
                       , u * u                     , (-(u * u))     , 0.0
                       , (-u)                      , u              , (-1.0) / eps - u
                       ]
  where
    y = toList x
    u = y !! 0
    v = y !! 1
    w = y !! 2
    eps = 5.0e-6

brussRoot :: CV.SolverResult
brussRoot = CV.odeSolveRootVWith' opts
                      (\t v -> vector $ brusselator t (toList v))
                      Nothing
                      (vector [1.2, 3.1, 3.0])
                      events 100
                      (vector [0.0, 0.1 .. 10.0])
  where
    events =
      [ CV.EventSpec { CV.eventCondition = brussRootFn
                     , CV.eventUpdate =
                         \_ev x -> let y = toList x in vector [(y!!0) + 0.5 , (y!!1), (y!!2)]
                     , CV.eventDirection = CV.AnyDirection
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

exponential :: CV.SolverResult
exponential = CV.odeSolveRootVWith' opts
                      (\t y -> vector [y ! 0])
                      Nothing
                      (vector [1])
                      events 100
                      (vector [ fromIntegral k / 100 | k <- [0..(22::Int)]])
  where
    events =
      [ CV.EventSpec { CV.eventCondition = \t y -> y ! 0 - 1.1
                     , CV.eventUpdate = \ev y -> vector [ 2 ]
                     , CV.eventDirection = CV.Upwards
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
boundedSine :: CV.SolverResult
boundedSine = CV.odeSolveRootVWith'
  opts
  (\_t y -> vector [y ! 1, - y ! 0]) -- ODE RHS
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
      [ CV.EventSpec { CV.eventCondition = \_t y -> y ! 0 - 0.9
                     , CV.eventUpdate = \_ y -> vector [ y ! 0, - abs (y ! 1) ]
                     , CV.eventDirection = CV.Upwards
                     }
      , CV.EventSpec { CV.eventCondition = \_t y -> y ! 0 + 0.9
                     , CV.eventUpdate = \_ y -> vector [ y ! 0, abs (y ! 1) ]
                     , CV.eventDirection = CV.Downwards
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

roberts :: Double -> Vector Double -> Vector Double
roberts t v = vector $ robertsAux t (toList v)
  where
    robertsAux _ [y1, y2, y3] =
      [ -0.04 * y1 + 1.0e4 * y2 * y3
      , 0.04 * y1 - 1.0e4 * y2 * y3 - 3.0e7 * (y2)^(2 :: Int)
      , 3.0e7 * (y2)^(2 :: Int)
      ]
    robertsAux _ _ = error "roberts RHS not defined"

ts :: [Double]
ts = take 12 $ map (* 10.0) (0.04 : ts)

solve :: CV.SolverResult
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
      [ CV.EventSpec { CV.eventCondition = \_t y -> y ! 0 - 0.0001
                     , CV.eventUpdate = const id
                     , CV.eventDirection = CV.AnyDirection
                     }
      , CV.EventSpec { CV.eventCondition = \_t y -> y ! 2 - 0.01
                     , CV.eventUpdate = const id
                     , CV.eventDirection = CV.AnyDirection
                     }
      ]

solve2 :: CV.SolverResult
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
      [ CV.EventSpec { CV.eventCondition = \_t y -> y ! 0 - 0.0001
                     , CV.eventUpdate = upd
                     , CV.eventDirection = CV.AnyDirection
                     }
      , CV.EventSpec { CV.eventCondition = \_t y -> y ! 2 - 0.01
                     , CV.eventUpdate = upd
                     , CV.eventDirection = CV.AnyDirection
                     }
      ]
    upd _ _ = vector [1.0, 0.0, 0.0]

solve1 :: CV.SolverResult
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
      [ CV.EventSpec { CV.eventCondition = \t _y -> t - 1.0
                     , CV.eventUpdate = \t y -> vector [2.0, y!1, y!2]
                     , CV.eventDirection = CV.AnyDirection
                     }
      ]

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
  let res1 = ARK.odeSolve brusselator [1.2, 3.1, 3.0] (fromList [0.0, 0.1 .. 10.0])
  renderRasterific "diagrams/brusselator.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ lSaxis $ [0.0, 0.1 .. 10.0]:(toLists $ tr res1))

  let res1a = ARK.odeSolve brusselator [1.2, 3.1, 3.0] (fromList [0.0, 0.1 .. 10.0])
  renderRasterific "diagrams/brusselatorA.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ lSaxis $ [0.0, 0.1 .. 10.0]:(toLists $ tr res1a))

  let res2 = ARK.odeSolve stiffish [0.0] (fromList [0.0, 0.1 .. 10.0])
  renderRasterific "diagrams/stiffish.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip [0.0, 0.1 .. 10.0] (concat $ toLists res2))

  let res2a = ARK.odeSolveV (ARK.SDIRK_5_3_4') Nothing 1e-6 1e-10 stiffishV (fromList [0.0]) (fromList [0.0, 0.1 .. 10.0])

  let res2b = ARK.odeSolveV (ARK.TRBDF2_3_3_2') Nothing 1e-6 1e-10 stiffishV (fromList [0.0]) (fromList [0.0, 0.1 .. 10.0])

  let maxDiffA = maximum $ map abs $
                 zipWith (-) ((toLists $ tr res2a)!!0) ((toLists $ tr res2b)!!0)

  let res2c = CV.odeSolveV (CV.BDF) Nothing 1e-6 1e-10 stiffishV (fromList [0.0]) (fromList [0.0, 0.1 .. 10.0])

  let maxDiffB = maximum $ map abs $
                 zipWith (-) ((toLists $ tr res2a)!!0) ((toLists $ tr res2c)!!0)

  let maxDiffC = maximum $ map abs $
                 zipWith (-) ((toLists $ tr res2b)!!0) ((toLists $ tr res2c)!!0)

  let res3 = ARK.odeSolve lorenz [-5.0, -5.0, 1.0] (fromList [0.0, 0.01 .. 20.0])

  renderRasterific "diagrams/lorenz.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res3)!!0) ((toLists $ tr res3)!!1))

  renderRasterific "diagrams/lorenz1.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res3)!!0) ((toLists $ tr res3)!!2))

  renderRasterific "diagrams/lorenz2.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res3)!!1) ((toLists $ tr res3)!!2))

  let res4 = CV.odeSolve predatorPrey [0.5, 1.0, 2.0] (fromList [0.0, 0.01 .. 10.0])

  renderRasterific "diagrams/predatorPrey.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res4)!!0) ((toLists $ tr res4)!!1))

  renderRasterific "diagrams/predatorPrey1.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res4)!!0) ((toLists $ tr res4)!!2))

  renderRasterific "diagrams/predatorPrey2.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip ((toLists $ tr res4)!!1) ((toLists $ tr res4)!!2))

  let res4a = ARK.odeSolve predatorPrey [0.5, 1.0, 2.0] (fromList [0.0, 0.01 .. 10.0])

  let maxDiffPpA = maximum $ map abs $
                   zipWith (-) ((toLists $ tr res4)!!0) ((toLists $ tr res4a)!!0)

  let cond5 =
        case solve of
          CV.SolverSuccess events _ _ -> do
            length events `shouldBe` 2
            (abs (CV.eventTime (events!!0) - 0.2640208751331032) / 0.2640208751331032 < 1.0e-8) `shouldBe` True
            (abs (CV.eventTime (events!!1) - 2.0786731062254436e7) / 2.0786731062254436e7 < 1.0e-8) `shouldBe` True
          CV.SolverError _ _ ->
            error "Root finding error!"

  let cond6 =
        case solve1 of
          CV.SolverSuccess events _ _ -> do
            length events `shouldBe` 1
            (abs (CV.eventTime (events!!0) - 1.0) / 1.0 < 1.0e-10) `shouldBe` True
          CV.SolverError _ _ ->
            error "Root finding error!"

  let cond7 =
        case solve2 of
          CV.SolverSuccess {} ->
            error "Solver returned Success"
          CV.SolverError _ _ ->
            True

  case brussRoot of
    CV.SolverSuccess events m _diagn -> do
      renderRasterific
        "diagrams/brussRoot.png"
        (D.dims2D 500.0 500.0)
        (renderAxis $ lSaxis $ toLists $ tr m)
    CV.SolverError m n ->
      expectationFailure $ show n

  let boundedSineSpec = do
        case boundedSine of
          CV.SolverSuccess events m _ -> do
            renderRasterific
              "diagrams/boundedSine.png"
              (D.dims2D 500.0 500.0)
              (renderAxis $ lSaxis2 $ toLists $ tr m)
            length events `shouldBe` 3
            map CV.rootDirection events `shouldBe` [CV.Upwards, CV.Downwards, CV.Upwards]
            map CV.eventIndex events `shouldBe` [0, 1, 0]
            all ((< 1e-8) . abs) (zipWith (-)
              (map CV.eventTime events)
              [1.1197660081724263,3.3592952656818404,5.5988203973243])
                `shouldBe` True
          CV.SolverError m n ->
            expectationFailure "Solver error"
  let exponentialSpec = do
        case exponential of
          CV.SolverSuccess events _m _diagn -> do
            length events `shouldBe` 1
            (abs (CV.eventTime (events!!0) - log 1.1) < 1e-4) `shouldBe` True
            CV.rootDirection (events!!0) `shouldBe` CV.Upwards
            CV.eventIndex (events!!0) `shouldBe` 0
          CV.SolverError m n ->
            expectationFailure $ show n

  hspec $ describe "Compare results" $ do
    it "Robertson should fail" $ cond7
    it "Robertson time only" $ cond6
    it "Robertson from SUNDIALS manual" $ cond5
    it "for SDIRK_5_3_4' and TRBDF2_3_3_2'" $ maxDiffA < 1.0e-6
    it "for SDIRK_5_3_4' and BDF" $ maxDiffB < 1.0e-6
    it "for TRBDF2_3_3_2' and BDF" $ maxDiffC < 1.0e-6
    it "for CV and ARK for the Predator Prey model" $ maxDiffPpA < 1.0e-3
    it "Bounded sine events" $ boundedSineSpec
    it "Exponential events" $ exponentialSpec
