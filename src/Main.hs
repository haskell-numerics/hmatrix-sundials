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

brussRoot :: CV.SolverResult Matrix Vector (Compose [] []) Int Double
brussRoot = CV.odeSolveRootVWith' opts CV.BDF
                      (CV.XX' 1.0e-6 1.0e-10 1 1)
                      Nothing (\t v -> vector $ brusselator t (toList v))
                      (vector [1.2, 3.1, 3.0])
                      1 brussRootFn 100
                      (\_ev x -> let y = toList x in vector [(y!!0) + 0.5 , (y!!1), (y!!2)])
                      (vector [0.0, 0.1 .. 10.0])
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   }

brussRootFn :: Double -> Vector Double -> Vector Double
brussRootFn _ v = case xs of
                    [y1, _y2, y3] -> vector [ y1 - y3
                                            ]
                    _            -> error "brusselator root function RHS not defined"
  where
    xs = toList v

-- A sine wave that only changes direction once it reaches ±0.9.
-- Illustrates event-specific reset function
boundedSine :: CV.SolverResult Matrix Vector (Compose [] []) Int Double
boundedSine = CV.odeSolveRootVWith'
  opts
  CV.ADAMS -- Adams-Moulton multistep method
  (CV.XX' 1.0e-6 1.0e-10 1 1) -- adaptive step control
  Nothing -- initial step size: use the auto-calculated one
  (\_t y -> vector [y ! 1, - y ! 0]) -- ODE RHS
  (vector [0, 1]) -- initial conditions
  2 -- number of event equations
  (\_t y -> vector [ y ! 0 - 0.9, y ! 0 + 0.9 ]) -- event equations
  100 -- maximum number of events
  (\ev y -> vector [y ! 0, (if ev == 0 then -1 else 1) * abs (y ! 1)]) -- event handler
  (vector [ 2 * pi * k / 360 | k <- [0..360]]) -- solution times
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   }

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

rootFn :: Double -> Vector Double -> Vector Double
rootFn _ v = case xs of
               [y1, _y2, y3] -> vector [ y1 - 0.0001
                                       , y3 - 0.01
                                       ]
               _             -> error "roberts root function RHS not defined"
  where
    xs = toList v

rootFn1 :: Double -> Vector Double -> Vector Double
rootFn1 t _ = fromList [t - 1.0]

ts :: [Double]
ts = take 12 $ map (* 10.0) (0.04 : ts)

solve :: CV.SolverResult Matrix Vector (Compose [] []) Int Double
solve = CV.odeSolveRootVWith' opts CV.BDF
                      (CV.ScXX' 1.0 1.0e-4 1.0 1.0 (vector [1.0e-8, 1.0e-14, 1.0e-6]))
                      Nothing roberts (vector [1.0, 0.0, 0.0])
                      2 rootFn 100
                      (const id)
                      (vector (0.0 : ts))
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   }

solve2 :: CV.SolverResult Matrix Vector (Compose [] []) Int Double
solve2 = CV.odeSolveRootVWith' opts CV.BDF
                      (CV.ScXX' 1.0 1.0e-4 1.0 1.0 (vector [1.0e-8, 1.0e-14, 1.0e-6]))
                      Nothing roberts (vector [1.0, 0.0, 0.0])
                      2 rootFn 100
                      (const . const $ vector [1.0, 0.0, 0.0])
                      (vector (0.0 : ts))
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   }

solve1 :: CV.SolverResult Matrix Vector (Compose [] []) Int Double
solve1 = CV.odeSolveRootVWith' opts CV.BDF
                      (CV.ScXX' 1.0 1.0e-4 1.0 1.0 (vector [1.0e-8, 1.0e-14, 1.0e-6]))
                      Nothing roberts (vector [1.0, 0.0, 0.0])
                      1 rootFn1 100
                      (\_ev x -> let y = toList x in vector [2.0, y!!1, y!!2])
                      (vector (0.0 : ts))
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
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
          CV.SolverRoot rootTimes _ _ _ ->
            abs (rootTimes!0 - 0.2640208751331032) / 0.2640208751331032 < 1.0e-10 &&
            abs (rootTimes!1 - 2.0786731062254436e7) / 2.0786731062254436e7 < 1.0e-10
          CV.SolverSuccess _ _ ->
            error "No roots found!"
          CV.SolverError _ _ ->
            error "Root finding error!"

  let cond6 =
        case solve1 of
          CV.SolverRoot rootTimes _ _ _ ->
            abs (rootTimes!0 - 1.0) / 1.0 < 1.0e-10
          CV.SolverSuccess _ _ ->
            error "No roots found!"
          CV.SolverError _ _ ->
            error "Root finding error!"

  let cond7 =
        case solve2 of
          CV.SolverRoot _ _ _ _ ->
            error "Roots found!"
          CV.SolverSuccess _ _ ->
            error "No roots found!"
          CV.SolverError _ _ ->
            True

  case brussRoot of
    CV.SolverRoot _a _b m _c -> do
      renderRasterific
        "diagrams/brussRoot.png"
        (D.dims2D 500.0 500.0)
        (renderAxis $ lSaxis $ toLists $ tr m)
    CV.SolverSuccess m _ ->
      error $ "No roots found!\n" ++ show m
    CV.SolverError m n ->
      error $ "No roots found!\n" ++ show m ++ "\n" ++ show n

  case boundedSine of
    CV.SolverRoot _a _b m _c -> do
      renderRasterific
        "diagrams/boundedSine.png"
        (D.dims2D 500.0 500.0)
        (renderAxis $ lSaxis2 $ toLists $ tr m)
    CV.SolverSuccess m _ ->
      error $ "No roots found!\n" ++ show m
    CV.SolverError m n ->
      error $ "No roots found!\n" ++ show m ++ "\n" ++ show n

  hspec $ describe "Compare results" $ do
    it "Robertson should fail" $ cond7
    it "Robertson time only" $ cond6
    it "Robertson from SUNDIALS manual" $ cond5
    it "for SDIRK_5_3_4' and TRBDF2_3_3_2'" $ maxDiffA < 1.0e-6
    it "for SDIRK_5_3_4' and BDF" $ maxDiffB < 1.0e-6
    it "for TRBDF2_3_3_2' and BDF" $ maxDiffC < 1.0e-6
    it "for CV and ARK for the Predator Prey model" $ maxDiffPpA < 1.0e-3

