{-# OPTIONS_GHC -Wall -Wno-partial-type-signatures #-}

{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE PartialTypeSignatures #-}

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
                                   , odeSolveWithEvents
                                   , ODEMethod(..)
                                   , StepControl(..)
                                   , SolverResult(..)
                                   ) where

import qualified Language.C.Inline as C
import qualified Language.C.Inline.Unsafe as CU

import           Data.Monoid ((<>))
import           Data.Maybe (isJust, fromJust)
import           Data.List (genericLength)

import           Foreign.C.Types (CDouble, CInt, CLong)
import           Foreign.Ptr
import           Foreign.Storable (peek, poke)

import qualified Data.Vector.Storable as V

import           Data.Coerce (coerce)
import           System.IO.Unsafe (unsafePerformIO)

import           Numeric.LinearAlgebra.Devel (createVector)

import           Numeric.LinearAlgebra.HMatrix (Vector, Matrix, toList, rows,
                                                cols, toLists, size, reshape,
                                                subVector, subMatrix, toColumns, fromColumns, asColumn)

import           Numeric.Sundials.Arkode (cV_ADAMS, cV_BDF,
                                          vectorToC, cV_SUCCESS,
                                          SunVector(..))
import qualified Numeric.Sundials.Arkode as T
import           Numeric.Sundials.ODEOpts


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
  deriving (Eq, Ord, Show, Read)

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
  case odeSolveVWith' opts f y0 tt of
    Left  (c, _v) -> error $ show c -- FIXME
    Right (v, _d) -> v
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod   = method
                   , stepControl = control
                   , initStep    = initStepSize
                   }

odeSolveVWith' ::
  ODEOpts ODEMethod
  -> (Double -> V.Vector Double -> V.Vector Double) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector Double                     -- ^ Initial conditions
  -> V.Vector Double                     -- ^ Desired solution times
  -> Either (Matrix Double, Int) (Matrix Double, SundialsDiagnostics) -- ^ Error code or solution
odeSolveVWith' opts f y0 tt =
  case solveOdeC (fromIntegral $ maxFail opts)
                  (fromIntegral $ maxNumSteps opts) (coerce $ minStep opts)
                  (fromIntegral . getMethod . odeMethod $ opts) (coerce $ initStep opts) jacH (scise $ stepControl opts)
                  (OdeRhsHaskell $ coerce f) (coerce y0)
                  0 (\_ x -> x) [] 0 (\_ _ y -> y) (coerce tt) of
    -- Remove the time column for backwards compatibility
    SolverError m c         -> Left
                               ( subMatrix (0, 1) (V.length tt, l) m
                               , fromIntegral c
                               )
    SolverSuccess _ m d     -> Right
                               ( subMatrix (0, 1) (V.length tt, l) m
                               , d
                               )
  where
    l = size y0
    scise (X aTol rTol)                          = coerce (V.replicate l aTol, rTol)
    scise (X' aTol rTol)                         = coerce (V.replicate l aTol, rTol)
    scise (XX' aTol rTol yScale _yDotScale)      = coerce (V.replicate l aTol, yScale * rTol)
    -- FIXME; Should we check that the length of ss is correct?
    scise (ScXX' aTol rTol yScale _yDotScale ss) = coerce (V.map (* aTol) ss, yScale * rTol)
    jacH = fmap (\g t v -> matrixToSunMatrix $ g (coerce t) (coerce v)) $
           getJacobian $ odeMethod opts

matrixToSunMatrix :: Matrix Double -> T.SunMatrix
matrixToSunMatrix m = T.SunMatrix { T.rows = nr, T.cols = nc, T.vals = vs }
  where
    nr = fromIntegral $ rows m
    nc = fromIntegral $ cols m
    -- FIXME: efficiency
    vs = V.fromList $ map coerce $ concat $ toLists m

foreign import ccall "wrapper"
  mkOdeRhsC :: OdeRhsCType -> IO (FunPtr OdeRhsCType)

solveOdeC ::
  CInt ->
  CLong ->
  CDouble ->
  CInt ->
  Maybe CDouble ->
  (Maybe (CDouble -> V.Vector CDouble -> T.SunMatrix)) ->
  (V.Vector CDouble, CDouble) ->
  OdeRhs -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector CDouble -- ^ Initial conditions
  -> CInt -- ^ Number of event equations
  -> (CDouble -> V.Vector CDouble -> V.Vector CDouble) -- ^ The event equations themselves
  -> [CrossingDirection] -- ^ The required crossing direction for each event
  -> CInt -- ^ Maximum number of events
  -> (Int -> CDouble -> V.Vector CDouble -> V.Vector CDouble)
      -- ^ Function to reset/update the state when an event occurs. The
      -- 'Int' argument is the 0-based number of the event that has
      -- occurred. If multiple events have occurred at the same time, they
      -- are handled in the increasing order of the event index. The other
      -- arguments are the time and the point in the state space. Return
      -- the updated point in the state space.
  -> V.Vector CDouble -- ^ Desired solution times
  -> SolverResult
solveOdeC maxErrTestFails maxNumSteps_ minStep_ method initStepSize
          jacH (aTols, rTol) rhs f0 nr event_fn directions max_events apply_event ts
  | V.null f0 = -- 0-dimensional (empty) system
    SolverSuccess [] (asColumn (coerce ts)) emptyDiagnostics
  | otherwise =
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
  output_mat_mut :: V.MVector _ CDouble <- V.thaw =<< createVector ((1 + fromIntegral dim) * (fromIntegral (2 * max_events) + fromIntegral nTs))
  diagMut :: V.MVector _ CLong <- V.thaw =<< createVector 10 -- FIXME
  rhs_funptr :: FunPtr OdeRhsCType <-
    case rhs of
      OdeRhsC ptr -> return ptr
      OdeRhsHaskell fun -> do
        let
          funIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr () -> IO CInt
          funIO t y f _ptr = do
            sv <- peek y
            poke f $ SunVector { sunVecN = sunVecN sv
                               , sunVecVals = fun t (sunVecVals sv)
                               }
            return 0
        mkOdeRhsC funIO

  let nrPre = fromIntegral nr
  gResults :: V.Vector CInt <- createVector nrPre
  -- FIXME: Do we need to do this here? Maybe as it will get GC'd and
  -- we'd have to do a malloc in C otherwise :(
  gResMut <- V.thaw gResults
  event_index_mut :: V.MVector _ CInt <- V.thaw =<< createVector (fromIntegral max_events)
  event_time_mut :: V.MVector _ CDouble <- V.thaw =<< createVector (fromIntegral max_events)
  -- Total number of events. This is *not* directly re
  n_events_mut :: V.MVector _ CInt <- V.thaw =<< createVector 1
  -- Total number of rows in the output_mat_mut matrix. It *cannot* be
  -- inferred from n_events_mut because when an event occurs k times, it
  -- contributes k to n_events_mut but only 2 to n_rows_mut.
  n_rows_mut :: V.MVector _ CInt <- V.thaw =<< createVector 1
  actual_event_direction_mut :: V.MVector _ CInt <- V.thaw =<< createVector (fromIntegral max_events)

  let event_fn_c :: CDouble -> Ptr T.SunVector -> Ptr CDouble -> Ptr () -> IO CInt
      event_fn_c x y f _ptr = do
        vals <- event_fn x <$> (sunVecVals <$> peek y)
        -- FIXME: We should be able to use poke somehow
        vectorToC vals nrPre f
        return 0

      apply_event_c :: CInt -> CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> IO CInt
      apply_event_c event_index t y y' = do
        sv <- peek y
        poke y' $ SunVector
          { sunVecN = sunVecN sv
          , sunVecVals = apply_event (fromIntegral event_index) t (sunVecVals sv)
          }
        return 0

      requested_event_directions :: V.Vector CInt
      requested_event_directions = V.fromList $ map directionToInt directions

  let isJac :: CInt
      isJac = fromIntegral $ fromEnum $ isJust jacH
      jacIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunMatrix ->
               Ptr () -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunVector ->
               IO CInt
      jacIO t y _fy jacS _ptr _tmp1 _tmp2 _tmp3 = do
        case jacH of
          Nothing   -> error "Numeric.Sundials.CVode.ODE: Jacobian not defined"
          Just jacI -> do j <- jacI t <$> (sunVecVals <$> peek y)
                          poke jacS j
                          -- FIXME: I don't understand what this comment means
                          -- Unsafe since the function will be called many times.
                          [CU.exp| int{ 0 } |]

  res :: Int <- fromIntegral <$> [C.block| int {
                         /* general problem variables */

                         int flag;                  /* reusable error-checking flag                 */

                         int i, j;                  /* reusable loop indices                        */
                         N_Vector y = NULL;         /* empty vector for storing solution            */
                         N_Vector tv = NULL;        /* empty vector for storing absolute tolerances */

                         SUNMatrix A = NULL;        /* empty matrix for linear solver               */
                         SUNLinearSolver LS = NULL; /* empty linear solver object                   */
                         void *cvode_mem = NULL;    /* empty CVODE memory structure                 */
                         realtype t;
                         long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;

                         realtype tout;

                         /* input_ind tracks the current index into the ts array */
                         int input_ind = 1;
                         /* output_ind tracks the current row into the output_mat_mut matrix.
                            If differs from input_ind because of the extra rows corresponding to events. */
                         int output_ind = 1;
                         /* We need to update n_rows_mut every time we update output_ind because
                            of the possibility of early return (in which case we still need to assemble
                            the partial results matrix). We could even work with n_rows_mut only and ditch
                            output_ind, but the inline-c expression is quite verbose, and output_ind is
                            more convenient to use in index calculations.
                         */
                         ($vec-ptr:(int *n_rows_mut))[0] = output_ind;
                         /* event_ind tracks the current event number */
                         int event_ind = 0;

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
                         flag = CVodeInit(cvode_mem, $(int (* rhs_funptr) (double t, SunVector y[], SunVector dydt[], void * params)), T0, y);
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

                         /* Call CVodeRootInit to specify the root function event_fn_c with nr components */
                         flag = CVodeRootInit(cvode_mem, $(int nr), $fun:(int (* event_fn_c) (double t, SunVector y[], double gout[], void * params)));

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

                         /* Store initial conditions */
                         ($vec-ptr:(double *output_mat_mut))[0 * (NEQ + 1) + 0] = ($vec-ptr:(double *ts))[0];
                         for (j = 0; j < NEQ; j++) {
                           ($vec-ptr:(double *output_mat_mut))[0 * (NEQ + 1) + (j + 1)] = NV_Ith_S(y,j);
                         }

                         while (1) {
                           flag = CVode(cvode_mem, ($vec-ptr:(double *ts))[input_ind], y, &t, CV_NORMAL); /* call integrator */
                           if (check_flag(&flag, "CVode solver failure, stopping integration", 1)) return 1;

                           /* Store the results for Haskell */
                           ($vec-ptr:(double *output_mat_mut))[output_ind * (NEQ + 1) + 0] = t;
                           for (j = 0; j < NEQ; j++) {
                             ($vec-ptr:(double *output_mat_mut))[output_ind * (NEQ + 1) + (j + 1)] = NV_Ith_S(y,j);
                           }
                           output_ind++;
                           ($vec-ptr:(int *n_rows_mut))[0] = output_ind;

                           if (flag == CV_ROOT_RETURN) {
                             if (event_ind >= $(int max_events)) {
                               /* We reached the maximum number of events.
                                  Either the maximum number of events is set to 0,
                                  or there's a bug in our code below. In any case return an error.
                               */
                               return 1;
                             }

                             /* Are we interested in this event?
                                If not, continue without any observable side-effects.
                             */
                             int good_event = 0;
                             flag = CVodeGetRootInfo(cvode_mem, ($vec-ptr:(int *gResMut)));
                             if (check_flag(&flag, "CVodeGetRootInfo", 1)) return 1;
                             for (i = 0; i < $(int nr); i++) {
                               int ev = ($vec-ptr:(int *gResMut))[i];
                               int req_dir = ($vec-ptr:(const int *requested_event_directions))[i];
                               if (ev != 0 && ev * req_dir >= 0) {
                                 good_event = 1;

                                 ($vec-ptr:(int *actual_event_direction_mut))[event_ind] = ev;
                                 ($vec-ptr:(int *event_index_mut))[event_ind] = i;
                                 ($vec-ptr:(double *event_time_mut))[event_ind] = t;
                                 event_ind++;

                                 /* Update the state with the supplied function */
                                 $fun:(int (* apply_event_c) (int, double, SunVector y[], SunVector z[]))(i, t, y, y);
                               }
                             }

                             if (good_event) {
                               ($vec-ptr:(double *output_mat_mut))[output_ind * (NEQ + 1) + 0] = t;
                               for (j = 0; j < NEQ; j++) {
                                 ($vec-ptr:(double *output_mat_mut))[output_ind * (NEQ + 1) + (j + 1)] = NV_Ith_S(y,j);
                               }
                               output_ind++;
                               ($vec-ptr:(int *n_rows_mut))[0] = output_ind;

                               if (event_ind >= $(int max_events)) {
                                 /* We collected the requested number of events. Stop the solver. */
                                 break;
                               }
                               flag = CVodeReInit(cvode_mem, t, y);
                               if (check_flag(&flag, "CVodeReInit", 1)) return(1);
                             } else {
                               /* Since this is not a wanted event, it shouldn't get a row */
                               output_ind--;
                               ($vec-ptr:(int *n_rows_mut))[0] = output_ind;
                             }
                           }
                           else {
                             if (++input_ind >= $(int nTs))
                               break;
                           }
                         }

                         /* The number of actual roots we found */
                         ($vec-ptr:(int *n_events_mut))[0] = event_ind;

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

                         return CV_SUCCESS;
                       } |]

  -- Free the allocated FunPtr. Ideally this should be done within
  -- a bracket...
  case rhs of
    OdeRhsHaskell {} -> freeHaskellFunPtr rhs_funptr
    OdeRhsC {} -> return () -- we didn't allocate this

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
  n_rows <- fromIntegral . V.head <$> V.freeze n_rows_mut
  output_mat <- coerce . reshape (dim + 1) . subVector 0 ((dim + 1) * n_rows) <$>
    V.freeze output_mat_mut
  n_events <- fromIntegral . V.head <$> V.freeze n_events_mut
  event_time             :: V.Vector Double
    <- coerce . V.take n_events <$> V.freeze event_time_mut
  event_index            :: V.Vector Int
    <- V.map fromIntegral . V.take n_events <$> V.freeze event_index_mut
  actual_event_direction :: V.Vector CInt
    <- V.take n_events <$> V.freeze actual_event_direction_mut
  let
    events :: [EventInfo]
    events = zipWith3 EventInfo
      (V.toList event_time)
      (V.toList event_index)
      (map (fromJust . intToDirection) $ V.toList actual_event_direction)
  return $
    if res == cV_SUCCESS
      then
        SolverSuccess events output_mat d
      else
        SolverError output_mat res

data SolverResult
  = SolverError !(Matrix Double) !Int
      -- ^ Partial results and error code
  | SolverSuccess
      [EventInfo]
      !(Matrix Double)
      !SundialsDiagnostics
      -- ^ Times at which the event was triggered, information about which root and the
                                                   -- results and diagnostics.
    deriving Show

odeSolveRootVWith' ::
  ODEOpts ODEMethod
  -> OdeRhs
      -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> Maybe (Double -> Vector Double -> Matrix Double)
      -- ^ The Jacobian (optional)
  -> V.Vector Double                      -- ^ Initial conditions
  -> [EventSpec]                          -- ^ Event specifications
  -> Int                                  -- ^ Maximum number of events
  -> V.Vector Double                      -- ^ Desired solution times
  -> SolverResult
odeSolveRootVWith' opts rhs mb_jacobian y0 event_specs nRootEvs tt =
  solveOdeC (fromIntegral $ maxFail opts)
                 (fromIntegral $ maxNumSteps opts) (coerce $ minStep opts)
                 (fromIntegral . getMethod . odeMethod $ opts) (coerce $ initStep opts) jacH (scise $ stepControl opts)
                 rhs (coerce y0)
                 (genericLength event_specs) event_equations event_directions
                 (fromIntegral nRootEvs) reset_state
                 (coerce tt)
  where
    l = size y0
    scise (X aTol rTol)                          = coerce (V.replicate l aTol, rTol)
    scise (X' aTol rTol)                         = coerce (V.replicate l aTol, rTol)
    scise (XX' aTol rTol yScale _yDotScale)      = coerce (V.replicate l aTol, yScale * rTol)
    -- FIXME; Should we check that the length of ss is correct?
    scise (ScXX' aTol rTol yScale _yDotScale ss) = coerce (V.map (* aTol) ss, yScale * rTol)
    jacH = fmap (\g t v -> matrixToSunMatrix $ g (coerce t) (coerce v)) $ mb_jacobian
    event_equations :: CDouble -> Vector CDouble -> Vector CDouble
    event_equations t y = V.fromList $
      map (\ev -> coerce (eventCondition ev) t y) event_specs
    event_directions :: [CrossingDirection]
    event_directions = map eventDirection event_specs
    reset_state :: Int -> CDouble -> Vector CDouble -> Vector CDouble
    reset_state n_event = coerce $ eventUpdate (event_specs !! n_event)

odeSolveWithEvents
  :: ODEOpts ODEMethod
  -> [EventSpec]
    -- ^ Event specifications
  -> Int
    -- ^ Maximum number of events
  -> OdeRhs
    -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> Maybe (Double -> Vector Double -> Matrix Double)
    -- ^ The Jacobian (optional)
  -> V.Vector Double
    -- ^ Initial conditions
  -> V.Vector Double
    -- ^ Desired solution times
  -> Either Int SundialsSolution
    -- ^ Either an error code or a solution
odeSolveWithEvents opts event_specs max_events rhs mb_jacobian initial sol_times =
  let
    result :: SolverResult
    result =
      odeSolveRootVWith' opts rhs mb_jacobian initial event_specs
        max_events sol_times
  in
    case result of
      SolverError _ code -> Left code
      SolverSuccess events mx diagn ->
        Right $ SundialsSolution
            { actualTimeGrid = extractTimeGrid mx
            , solutionMatrix = dropTimeGrid mx
            , eventInfo = events
            , diagnostics = diagn
            }
  where
    -- The time grid is the first column of the result matrix
    extractTimeGrid :: Matrix Double -> Vector Double
    extractTimeGrid = head . toColumns
    dropTimeGrid :: Matrix Double -> Matrix Double
    dropTimeGrid = fromColumns . tail . toColumns
