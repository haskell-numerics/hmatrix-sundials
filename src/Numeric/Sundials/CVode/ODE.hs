{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}

-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Sundials.CVode.ODE
-- Copyright   :  Dominic Steinitz 2018,
--                Novadiscovery 2018
-- License     :  BSD
-- Maintainer  :  Dominic Steinitz
-- Stability   :  provisional
--
-- Solution of ordinary differential equation (ODE) initial value problems.
--
-- <https://computation.llnl.gov/projects/sundials/sundials-software>
--
-- A simple example:
--
-- <<diagrams/brusselator.png#diagram=brusselator&height=400&width=500>>
--
-- @
-- import           Numeric.Sundials.CVode.ODE
-- import           Numeric.LinearAlgebra
--
-- import           Plots as P
-- import qualified Diagrams.Prelude as D
-- import           Diagrams.Backend.Rasterific
--
-- brusselator :: Double -> [Double] -> [Double]
-- brusselator _t x = [ a - (w + 1) * u + v * u * u
--                    , w * u - v * u * u
--                    , (b - w) / eps - w * u
--                    ]
--   where
--     a = 1.0
--     b = 3.5
--     eps = 5.0e-6
--     u = x !! 0
--     v = x !! 1
--     w = x !! 2
--
-- lSaxis :: [[Double]] -> P.Axis B D.V2 Double
-- lSaxis xs = P.r2Axis &~ do
--   let ts = xs!!0
--       us = xs!!1
--       vs = xs!!2
--       ws = xs!!3
--   P.linePlot' $ zip ts us
--   P.linePlot' $ zip ts vs
--   P.linePlot' $ zip ts ws
--
-- main = do
--   let res1 = odeSolve brusselator [1.2, 3.1, 3.0] (fromList [0.0, 0.1 .. 10.0])
--   renderRasterific "diagrams/brusselator.png"
--                    (D.dims2D 500.0 500.0)
--                    (renderAxis $ lSaxis $ [0.0, 0.1 .. 10.0]:(toLists $ tr res1))
-- @
--
-----------------------------------------------------------------------------
module Numeric.Sundials.CVode.ODE ( odeSolve
                                   , odeSolveV
                                   , odeSolveVWith
                                   , odeSolveVWith'
                                   , odeSolveRootVWith'
                                   , ODEMethod(..)
                                   , StepControl(..)
                                   , SolverResult(..)
                                   ) where

import qualified Language.C.Inline as C
import qualified Language.C.Inline.Unsafe as CU

import           Data.Monoid ((<>))
import           Data.Maybe (isJust)
import           Data.List.Split (chunksOf)
import           Data.Functor.Compose

import           Foreign.C.Types (CDouble, CInt, CLong)
import           Foreign.Ptr (Ptr)
import           Foreign.Storable (poke)

import qualified Data.Vector.Storable as V

import           Data.Coerce (coerce)
import           System.IO.Unsafe (unsafePerformIO)

import           Numeric.LinearAlgebra.Devel (createVector)

import           Numeric.LinearAlgebra.HMatrix (Vector, Matrix, toList, rows,
                                                cols, toLists, size, reshape,
                                                subVector, subMatrix)

import           Numeric.Sundials.Arkode (cV_ADAMS, cV_BDF,
                                          getDataFromContents, putDataInContents,
                                          vectorToC, cV_SUCCESS, cV_ROOT_RETURN)
import qualified Numeric.Sundials.Arkode as T
import           Numeric.Sundials.ODEOpts (ODEOpts(..), Jacobian, SundialsDiagnostics(..))


C.context (C.baseCtx <> C.vecCtx <> C.funCtx <> T.sunCtx)

C.include "<stdlib.h>"
C.include "<stdio.h>"
C.include "<math.h>"
C.include "<cvode/cvode.h>"               -- prototypes for CVODE fcts., consts.
C.include "<nvector/nvector_serial.h>"    -- serial N_Vector types, fcts., macros
C.include "<sunmatrix/sunmatrix_dense.h>" -- access to dense SUNMatrix
C.include "<sunlinsol/sunlinsol_dense.h>" -- access to dense SUNLinearSolver
C.include "<cvode/cvode_direct.h>"        -- access to CVDls interface
C.include "<sundials/sundials_types.h>"   -- definition of type realtype
C.include "<sundials/sundials_math.h>"
C.include "../../../helpers.h"
C.include "Numeric/Sundials/Arkode_hsc.h"


-- | Stepping functions
data ODEMethod = ADAMS
               | BDF

getMethod :: ODEMethod -> Int
getMethod (ADAMS) = cV_ADAMS
getMethod (BDF)   = cV_BDF

getJacobian :: ODEMethod -> Maybe Jacobian
getJacobian _ = Nothing

-- | A version of 'odeSolveVWith' with reasonable default step control.
odeSolveV
    :: ODEMethod
    -> Maybe Double      -- ^ initial step size - by default, CVode
                         -- estimates the initial step size to be the
                         -- solution \(h\) of the equation
                         -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
                         -- \(\ddot{y}\) is an estimated value of the
                         -- second derivative of the solution at \(t_0\)
    -> Double            -- ^ absolute tolerance for the state vector
    -> Double            -- ^ relative tolerance for the state vector
    -> (Double -> Vector Double -> Vector Double) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
    -> Vector Double     -- ^ initial conditions
    -> Vector Double     -- ^ desired solution times
    -> Matrix Double     -- ^ solution
odeSolveV meth hi epsAbs epsRel f y0 ts =
  odeSolveVWith meth (X epsAbs epsRel) hi g y0 ts
  where
    g t x0 = coerce $ f t x0

-- | A version of 'odeSolveV' with reasonable default parameters and
-- system of equations defined using lists. FIXME: we should say
-- something about the fact we could use the Jacobian but don't for
-- compatibility with hmatrix-gsl.
odeSolve :: (Double -> [Double] -> [Double]) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
         -> [Double]                         -- ^ initial conditions
         -> Vector Double                    -- ^ desired solution times
         -> Matrix Double                    -- ^ solution
odeSolve f y0 ts =
  -- FIXME: These tolerances are different from the ones in GSL
  odeSolveVWith BDF (XX' 1.0e-6 1.0e-10 1 1)  Nothing g (V.fromList y0) (V.fromList $ toList ts)
  where
    g t x0 = V.fromList $ f t (V.toList x0)

-- | A version of 'odeSolveVWith'' with reasonable default solver
-- options.
odeSolveVWith ::
  ODEMethod
  -> StepControl
  -> Maybe Double -- ^ initial step size - by default, CVode
                  -- estimates the initial step size to be the
                  -- solution \(h\) of the equation
                  -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
                  -- \(\ddot{y}\) is an estimated value of the second
                  -- derivative of the solution at \(t_0\)
  -> (Double -> V.Vector Double -> V.Vector Double) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector Double                     -- ^ Initial conditions
  -> V.Vector Double                     -- ^ Desired solution times
  -> Matrix Double                       -- ^ Error code or solution
odeSolveVWith method control initStepSize f y0 tt =
  case odeSolveVWith' opts method control initStepSize f y0 tt of
    Left  (c, _v) -> error $ show c -- FIXME
    Right (v, _d) -> v
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   }

odeSolveVWith' ::
  ODEOpts
  -> ODEMethod
  -> StepControl
  -> Maybe Double -- ^ initial step size - by default, CVode
                  -- estimates the initial step size to be the
                  -- solution \(h\) of the equation
                  -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
                  -- \(\ddot{y}\) is an estimated value of the second
                  -- derivative of the solution at \(t_0\)
  -> (Double -> V.Vector Double -> V.Vector Double) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector Double                     -- ^ Initial conditions
  -> V.Vector Double                     -- ^ Desired solution times
  -> Either (Matrix Double, Int) (Matrix Double, SundialsDiagnostics) -- ^ Error code or solution
odeSolveVWith' opts method control initStepSize f y0 tt =
  case solveOdeC (fromIntegral $ maxFail opts)
                  (fromIntegral $ maxNumSteps opts) (coerce $ minStep opts)
                  (fromIntegral $ getMethod method) (coerce initStepSize) jacH (scise control)
                  (coerce f) (coerce y0)
                  0 (\_ x -> x) 0 id (coerce tt) of
    -- Remove the time column for backwards compatibility
    SolverError v c         -> Left
                               ( subMatrix (0, 1) (V.length tt, l) (reshape (l + 1) (coerce v))
                               , fromIntegral c
                               )
    SolverSuccess v d       -> Right
                               ( subMatrix (0, 1) (V.length tt, l) (reshape (l + 1) (coerce v))
                               , d
                               )
    SolverRoot _t _rs _v _d -> error "Roots found with no root equations!"
  where
    l = size y0
    scise (X aTol rTol)                          = coerce (V.replicate l aTol, rTol)
    scise (X' aTol rTol)                         = coerce (V.replicate l aTol, rTol)
    scise (XX' aTol rTol yScale _yDotScale)      = coerce (V.replicate l aTol, yScale * rTol)
    -- FIXME; Should we check that the length of ss is correct?
    scise (ScXX' aTol rTol yScale _yDotScale ss) = coerce (V.map (* aTol) ss, yScale * rTol)
    jacH = fmap (\g t v -> matrixToSunMatrix $ g (coerce t) (coerce v)) $
           getJacobian method

matrixToSunMatrix :: Matrix Double -> T.SunMatrix
matrixToSunMatrix m = T.SunMatrix { T.rows = nr, T.cols = nc, T.vals = vs }
  where
    nr = fromIntegral $ rows m
    nc = fromIntegral $ cols m
    -- FIXME: efficiency
    vs = V.fromList $ map coerce $ concat $ toLists m

solveOdeC ::
  CInt ->
  CLong ->
  CDouble ->
  CInt ->
  Maybe CDouble ->
  (Maybe (CDouble -> V.Vector CDouble -> T.SunMatrix)) ->
  (V.Vector CDouble, CDouble) ->
  (CDouble -> V.Vector CDouble -> V.Vector CDouble) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector CDouble -- ^ Initial conditions
  -> CInt -- ^ Number of event equations
  -> (CDouble -> V.Vector CDouble -> V.Vector CDouble) -- ^ The event equations themselves
  -> CInt -- ^ Maximum number of events
  -> (V.Vector CDouble -> V.Vector CDouble)
  -> V.Vector CDouble -- ^ Desired solution times
  -> SolverResult V.Vector V.Vector (Compose [] []) CInt CDouble
solveOdeC maxErrTestFails maxNumSteps_ minStep_ method initStepSize
          jacH (aTols, rTol) fun f0 nr g nRootEvs resetFun ts =
  unsafePerformIO $ do

  let isInitStepSize :: CInt
      isInitStepSize = fromIntegral $ fromEnum $ isJust initStepSize
      ss :: CDouble
      ss = case initStepSize of
             -- It would be better to put an error message here but
             -- inline-c seems to evaluate this even if it is never
             -- used :(
             Nothing -> 0.0
             Just x  -> x

  let dim = V.length f0
      nEq :: CLong
      nEq = fromIntegral dim
      nTs :: CInt
      nTs = fromIntegral $ V.length ts
  quasiMatrixRes <- createVector ((1 + fromIntegral dim) * (fromIntegral (2 * nRootEvs) + fromIntegral nTs))
  qMatMut <- V.thaw quasiMatrixRes
  diagnostics :: V.Vector CLong <- createVector 10 -- FIXME
  diagMut <- V.thaw diagnostics
  -- We need the types that sundials expects.
  -- FIXME: The Haskell type is currently empty!
  let funIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr () -> IO CInt
      funIO t y f _ptr = do
        -- Convert the pointer we get from C (y) to a vector, and then
        -- apply the user-supplied function.
        fImm <- fun t <$> getDataFromContents dim y
        -- Fill in the provided pointer with the resulting vector.
        putDataInContents fImm dim f
        -- FIXME: I don't understand what this comment means
        -- Unsafe since the function will be called many times.
        [CU.exp| int{ 0 } |]

  let nrPre = fromIntegral nr
  gResults :: V.Vector CInt <- createVector nrPre
  gResultss :: V.Vector CInt <- createVector $ nrPre * fromIntegral nRootEvs
  -- FIXME: Do we need to do this here? Maybe as it will get GC'd and
  -- we'd have to do a malloc in C otherwise :(
  gResMut <- V.thaw gResults
  gRessMut <- V.thaw gResultss
  tRoot :: V.Vector CDouble <- createVector $ fromIntegral nRootEvs
  tRootMut <- V.thaw tRoot
  nCrossings :: V.Vector CInt <- createVector 1
  nCrossingsMut <- V.thaw nCrossings

  let gIO :: CDouble -> Ptr T.SunVector -> Ptr CDouble -> Ptr () -> IO CInt
      gIO x y f _ptr = do
        -- Convert the pointer we get from C (y) to a vector, and then
        -- apply the user-supplied function.
        gImm <- g x <$> getDataFromContents dim y
        -- Fill in the provided pointer with the resulting vector.
        vectorToC gImm nrPre f
        -- FIXME: I don't understand what this comment means
        -- Unsafe since the function will be called many times.
        [CU.exp| int{ 0 } |]

  let rIO :: Ptr T.SunVector -> Ptr T.SunVector -> IO CInt
      rIO y f = do
        -- Convert the pointer we get from C (y) to a vector, and then
        -- apply the user-supplied function.
        rImm <- resetFun <$> getDataFromContents dim y
        -- Fill in the provided pointer with the resulting vector.
        putDataInContents rImm dim f
        -- FIXME: I don't understand what this comment means
        -- Unsafe since the function will be called many times.
        [CU.exp| int{ 0 } |]

  let isJac :: CInt
      isJac = fromIntegral $ fromEnum $ isJust jacH
      jacIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunMatrix ->
               Ptr () -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunVector ->
               IO CInt
      jacIO t y _fy jacS _ptr _tmp1 _tmp2 _tmp3 = do
        case jacH of
          Nothing   -> error "Numeric.Sundials.CVode.ODE: Jacobian not defined"
          Just jacI -> do j <- jacI t <$> getDataFromContents dim y
                          poke jacS j
                          -- FIXME: I don't understand what this comment means
                          -- Unsafe since the function will be called many times.
                          [CU.exp| int{ 0 } |]

  res <- [C.block| int {
                         /* general problem variables */

                         int flag;                  /* reusable error-checking flag                 */
                         int flagr;                 /* root finding flag                            */

                         int i, j, k, l;            /* reusable loop indices                        */
                         N_Vector y = NULL;         /* empty vector for storing solution            */
                         N_Vector tv = NULL;        /* empty vector for storing absolute tolerances */

                         SUNMatrix A = NULL;        /* empty matrix for linear solver               */
                         SUNLinearSolver LS = NULL; /* empty linear solver object                   */
                         void *cvode_mem = NULL;    /* empty CVODE memory structure                 */
                         realtype t;
                         long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;

                         realtype tout;

                         /* general problem parameters */

                         realtype T0 = RCONST(($vec-ptr:(double *ts))[0]); /* initial time              */
                         sunindextype NEQ = $(sunindextype nEq);           /* number of dependent vars. */

                         /* Initialize data structures */

                         y = N_VNew_Serial(NEQ); /* Create serial vector for solution */
                         if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
                         /* Specify initial condition */
                         for (i = 0; i < NEQ; i++) {
                           NV_Ith_S(y,i) = ($vec-ptr:(double *f0))[i];
                         };

                         cvode_mem = CVodeCreate($(int method), CV_NEWTON);
                         if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

                         /* Call CVodeInit to initialize the integrator memory and specify the
                          * user's right hand side function in y'=f(t,y), the inital time T0, and
                          * the initial dependent variable vector y. */
                         flag = CVodeInit(cvode_mem, $fun:(int (* funIO) (double t, SunVector y[], SunVector dydt[], void * params)), T0, y);
                         if (check_flag(&flag, "CVodeInit", 1)) return(1);

                         tv = N_VNew_Serial(NEQ); /* Create serial vector for absolute tolerances */
                         if (check_flag((void *)tv, "N_VNew_Serial", 0)) return 1;
                         /* Specify tolerances */
                         for (i = 0; i < NEQ; i++) {
                           NV_Ith_S(tv,i) = ($vec-ptr:(double *aTols))[i];
                         };

                         flag = CVodeSetMinStep(cvode_mem, $(double minStep_));
                         if (check_flag(&flag, "CVodeSetMinStep", 1)) return 1;
                         flag = CVodeSetMaxNumSteps(cvode_mem, $(long int maxNumSteps_));
                         if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return 1;
                         flag = CVodeSetMaxErrTestFails(cvode_mem, $(int maxErrTestFails));
                         if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1)) return 1;

                         /* Call CVodeSVtolerances to specify the scalar relative tolerance
                          * and vector absolute tolerances */
                         flag = CVodeSVtolerances(cvode_mem, $(double rTol), tv);
                         if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

                         /* Call CVodeRootInit to specify the root function g with nr components */
                         flag = CVodeRootInit(cvode_mem, $(int nr), $fun:(int (* gIO) (double t, SunVector y[], double gout[], void * params)));

                         if (check_flag(&flag, "CVodeRootInit", 1)) return(1);

                         /* Initialize dense matrix data structure and solver */
                         A = SUNDenseMatrix(NEQ, NEQ);
                         if (check_flag((void *)A, "SUNDenseMatrix", 0)) return 1;
                         LS = SUNDenseLinearSolver(y, A);
                         if (check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return 1;

                         /* Attach matrix and linear solver */
                         flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
                         if (check_flag(&flag, "CVDlsSetLinearSolver", 1)) return 1;

                         /* Set the initial step size if there is one */
                         if ($(int isInitStepSize)) {
                           /* FIXME: We could check if the initial step size is 0 */
                           /* or even NaN and then throw an error                 */
                           flag = CVodeSetInitStep(cvode_mem, $(double ss));
                           if (check_flag(&flag, "CVodeSetInitStep", 1)) return 1;
                         }

                         /* Set the Jacobian if there is one */
                         if ($(int isJac)) {
                           flag = CVDlsSetJacFn(cvode_mem, $fun:(int (* jacIO) (double t, SunVector y[], SunVector fy[], SunMatrix Jac[], void * params, SunVector tmp1[], SunVector tmp2[], SunVector tmp3[])));
                           if (check_flag(&flag, "CVDlsSetJacFn", 1)) return 1;
                         }

                         /* FIXME: These only work by accident */
                         /* Store initial conditions */
			 ($vec-ptr:(double *qMatMut))[0 * NEQ + 0] = ($vec-ptr:(double *ts))[0];
                         for (j = 0; j < NEQ; j++) {
                           ($vec-ptr:(double *qMatMut))[0 * NEQ + (j + 1)] = NV_Ith_S(y,j);
                         }

                         /* FIXME: This comment is no longer correct */
                         /* Main time-stepping loop: calls CVode to perform the integration */
                         /* Stops when the final time has been reached                      */
                         i = 1; k = 0;
                         while (1) {
                           flag = CVode(cvode_mem, ($vec-ptr:(double *ts))[i], y, &t, CV_NORMAL); /* call integrator */
                           if (check_flag(&flag, "CVode solver failure, stopping integration", 1)) return 1;

                           /* Store the results for Haskell */
			   ($vec-ptr:(double *qMatMut))[(i + k) * (NEQ + 1) + 0] = t;
                           for (j = 0; j < NEQ; j++) {
                             ($vec-ptr:(double *qMatMut))[(i + k) * (NEQ + 1) + (j + 1)] = NV_Ith_S(y,j);
                           }

                           if (flag == CV_ROOT_RETURN) {
			     ($vec-ptr:(double *tRootMut))[k / 2] = t;
			     flagr = CVodeGetRootInfo(cvode_mem, ($vec-ptr:(int *gResMut)));
			     for (l = 0; l < $(int nr); l++) {
			       ($vec-ptr:(int *gRessMut))[l + k * $(int nr) / 2] = ($vec-ptr:(int *gResMut))[l];
			     }
			     if (check_flag(&flagr, "CVodeGetRootInfo", 1)) return(1);
			     flagr = flag;

                             /* Update the state with the supplied function */
                             $fun:(int (* rIO) (SunVector y[], SunVector z[]))(y, y);

			     ($vec-ptr:(double *qMatMut))[(i + k  + 1) * (NEQ + 1) + 0] = t;
			     for (j = 0; j < NEQ; j++) {
			       ($vec-ptr:(double *qMatMut))[(i + k + 1) * (NEQ + 1) + (j + 1)] = NV_Ith_S(y,j);
			     }

                             flag = CVodeReInit(cvode_mem, t, y);
			     if (check_flag(&flag, "CVodeReInit", 1)) return(1);
                             if (k > 2 * $(int nRootEvs)) break; else k += 2;
                           }
			   else {
			     i++;
			     if (i >= $(int nTs)) break;
			   }
                         }

			 /* The number of actual roots we found */
			 ($vec-ptr:(int *nCrossingsMut))[0] = k / 2;

                         /* Get some final statistics on how the solve progressed */
                         flag = CVodeGetNumSteps(cvode_mem, &nst);
                         check_flag(&flag, "CVodeGetNumSteps", 1);
                         ($vec-ptr:(long int *diagMut))[0] = nst;

                         /* FIXME */
                         ($vec-ptr:(long int *diagMut))[1] = 0;

                         flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
                         check_flag(&flag, "CVodeGetNumRhsEvals", 1);
                         ($vec-ptr:(long int *diagMut))[2] = nfe;
                         /* FIXME */
                         ($vec-ptr:(long int *diagMut))[3] = 0;

                         flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
                         check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
                         ($vec-ptr:(long int *diagMut))[4] = nsetups;

                         flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
                         check_flag(&flag, "CVodeGetNumErrTestFails", 1);
                         ($vec-ptr:(long int *diagMut))[5] = netf;

                         flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
                         check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
                         ($vec-ptr:(long int *diagMut))[6] = nni;

                         flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
                         check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
                         ($vec-ptr:(long int *diagMut))[7] = ncfn;

                         flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
                         check_flag(&flag, "CVDlsGetNumJacEvals", 1);
                         ($vec-ptr:(long int *diagMut))[8] = ncfn;

                         flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
                         check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
                         ($vec-ptr:(long int *diagMut))[9] = ncfn;

                         /* Clean up and return */

                         N_VDestroy(y);          /* Free y vector          */
                         N_VDestroy(tv);         /* Free tv vector         */
                         CVodeFree(&cvode_mem);  /* Free integrator memory */
                         SUNLinSolFree(LS);      /* Free linear solver     */
                         SUNMatDestroy(A);       /* Free A matrix          */

                         if (flag == CV_SUCCESS && flagr == CV_ROOT_RETURN) {
                           return CV_ROOT_RETURN;
                         }
                         else {
                           return flag;
                         }
                       } |]
  preD <- V.freeze diagMut
  let d = SundialsDiagnostics (fromIntegral $ preD V.!0)
                              (fromIntegral $ preD V.!1)
                              (fromIntegral $ preD V.!2)
                              (fromIntegral $ preD V.!3)
                              (fromIntegral $ preD V.!4)
                              (fromIntegral $ preD V.!5)
                              (fromIntegral $ preD V.!6)
                              (fromIntegral $ preD V.!7)
                              (fromIntegral $ preD V.!8)
                              (fromIntegral $ preD V.!9)
  m  <- V.freeze qMatMut
  t  <- V.freeze tRootMut
  rs <- V.freeze gRessMut
  n <-  V.freeze nCrossingsMut
  let f r | r == cV_SUCCESS     = SolverSuccess (subVector 0 (fromIntegral (fromIntegral (dim + 1) * nTs)) m) d
          | r == cV_ROOT_RETURN =
              if (n V.! 0) <= nRootEvs
                then
                  SolverRoot (V.take (fromIntegral (n V.! 0)) t)
                             (Compose $ take (fromIntegral (n V.! 0)) $
                              chunksOf (fromIntegral nr) (V.toList rs))
                             (subVector 0 (fromIntegral (fromIntegral (dim + 1) * (nTs + (n V.! 0)))) m)
                             d
                else SolverError m (fromIntegral r)
          | otherwise           = SolverError m res
  return $ f $ fromIntegral res

data SolverResult f g h a b =
    SolverError (f b) a                            -- ^ Partial results and error code
  | SolverSuccess (f b) SundialsDiagnostics        -- ^ Results and diagnostics
  | SolverRoot (g b) (h a) (f b) SundialsDiagnostics   -- ^ Times at which the root was found, information about which root and the
                                                   -- results and diagnostics.
    deriving Show

odeSolveRootVWith' ::
  ODEOpts
  -> ODEMethod
  -> StepControl
  -> Maybe Double -- ^ initial step size - by default, CVode
                  -- estimates the initial step size to be the
                  -- solution \(h\) of the equation
                  -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
                  -- \(\ddot{y}\) is an estimated value of the second
                  -- derivative of the solution at \(t_0\)
  -> (Double -> V.Vector Double -> V.Vector Double) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector Double                     -- ^ Initial conditions
  -> Int                                 -- ^ Dimension of the range of the roots function
  -> (Double -> V.Vector Double -> V.Vector Double) -- ^ Roots function
  -> Int                                  -- ^ Maximum number of events
  -> (V.Vector Double -> V.Vector Double) -- ^ Function to reset the state
  -> V.Vector Double                      -- ^ Desired solution times
  -> SolverResult Matrix Vector (Compose [] []) Int Double
odeSolveRootVWith' opts method control initStepSize f y0 is gg nRootEvs hh tt =
  case solveOdeC (fromIntegral $ maxFail opts)
                 (fromIntegral $ maxNumSteps opts) (coerce $ minStep opts)
                 (fromIntegral $ getMethod method) (coerce initStepSize) jacH (scise control)
                 (coerce f) (coerce y0) (fromIntegral is) (coerce gg) (fromIntegral nRootEvs) (coerce hh)
                 (coerce tt) of
    SolverError v c     -> SolverError
                           (reshape l1 (coerce v)) (fromIntegral c)
    SolverSuccess v d   -> SolverSuccess
                           (reshape l1 (coerce v)) d
    SolverRoot t rs v d -> SolverRoot (coerce t) (Compose $ map (map fromIntegral) $ getCompose rs)
                           (reshape l1 (coerce v)) d
  where
    l = size y0
    l1 = l + 1 -- one more for the time column
    scise (X aTol rTol)                          = coerce (V.replicate l aTol, rTol)
    scise (X' aTol rTol)                         = coerce (V.replicate l aTol, rTol)
    scise (XX' aTol rTol yScale _yDotScale)      = coerce (V.replicate l aTol, yScale * rTol)
    -- FIXME; Should we check that the length of ss is correct?
    scise (ScXX' aTol rTol yScale _yDotScale ss) = coerce (V.map (* aTol) ss, yScale * rTol)
    jacH = fmap (\g t v -> matrixToSunMatrix $ g (coerce t) (coerce v)) $
           getJacobian method

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
-- allows the user to control the step size adjustment using
-- \(\eta^{rel}|y_i| + \eta^{abs}_i\). For compatibility with
-- [hmatrix-gsl](https://hackage.haskell.org/package/hmatrix-gsl),
-- tolerances for \(y\) and \(\dot{y}\) can be specified but the latter have no
-- effect.
data StepControl = X     Double Double -- ^ absolute and relative tolerance for \(y\); in GSL terms, \(a_{y} = 1\) and \(a_{dy/dt} = 0\); in ARKode terms, the \(\eta^{abs}_i\) are identical
                 | X'    Double Double -- ^ absolute and relative tolerance for \(\dot{y}\); in GSL terms, \(a_{y} = 0\) and \(a_{dy/dt} = 1\); in ARKode terms, the latter is treated as the relative tolerance for \(y\) so this is the same as specifying 'X' which may be entirely incorrect for the given problem
                 | XX'   Double Double Double Double -- ^ include both via relative tolerance
                                                     -- scaling factors \(a_y\), \(a_{{dy}/{dt}}\); in ARKode terms, the latter is ignored and \(\eta^{rel} = a_{y}\epsilon^{rel}\)
                 | ScXX' Double Double Double Double (Vector Double) -- ^ scale absolute tolerance of \(y_i\); in ARKode terms, \(a_{{dy}/{dt}}\) is ignored, \(\eta^{abs}_i = s_i \epsilon^{abs}\) and \(\eta^{rel} = a_{y}\epsilon^{rel}\)
