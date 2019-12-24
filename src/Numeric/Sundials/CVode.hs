-- | Solution of ordinary differential equation (ODE) initial value problems.
--
-- <https://computation.llnl.gov/projects/sundials/sundials-software>
module Numeric.Sundials.CVode
  ( CVMethod(..)
  , solveC
  ) where

import qualified Language.C.Inline as C
import qualified Data.Vector.Storable as VS
import Foreign.C.Types
import GHC.Prim

import Numeric.Sundials.Foreign
import Numeric.Sundials.Types
import Numeric.Sundials.Common

C.context (C.baseCtx <> C.vecCtx <> C.funCtx <> sunCtx)

C.include "<stdlib.h>"
C.include "<stdio.h>"
C.include "<string.h>"
C.include "<math.h>"
C.include "<arkode/arkode.h>"
C.include "<cvode/cvode.h>"               -- prototypes for CVODE fcts., consts.
C.include "<nvector/nvector_serial.h>"    -- serial N_Vector types, fcts., macros
C.include "<sunmatrix/sunmatrix_dense.h>" -- access to dense SUNMatrix
C.include "<sunlinsol/sunlinsol_dense.h>" -- access to dense SUNLinearSolver
C.include "<cvode/cvode_direct.h>"        -- access to CVDls interface
C.include "<sundials/sundials_types.h>"   -- definition of type realtype
C.include "<sundials/sundials_math.h>"
C.include "../../helpers.h"
C.include "Numeric/Sundials/Foreign_hsc.h"

-- | Available methods for CVode
data CVMethod = ADAMS
              | BDF
  deriving (Eq, Ord, Show, Read)

instance Method CVMethod where
  methodToInt ADAMS = cV_ADAMS
  methodToInt BDF   = cV_BDF
  methodSolver = CVode

solveC :: CConsts -> CVars (VS.MVector RealWorld) -> ReportErrorFn -> IO CInt
solveC CConsts{..} CVars{..} report_error =
  [C.block| int {
  /* general problem variables */

  int flag;                  /* reusable error-checking flag                 */

  int i, j;                  /* reusable loop indices                        */
  N_Vector y = NULL;         /* empty vector for storing solution            */
  N_Vector tv = NULL;        /* empty vector for storing absolute tolerances */

  SUNMatrix A = NULL;        /* empty matrix for linear solver               */
  SUNLinearSolver LS = NULL; /* empty linear solver object                   */
  void *cvode_mem = NULL;    /* empty CVODE memory structure                 */
  realtype t;
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf;

  /* input_ind tracks the current index into the c_sol_time array */
  int input_ind = 1;
  /* output_ind tracks the current row into the c_output_mat matrix.
     If differs from input_ind because of the extra rows corresponding to events. */
  int output_ind = 1;
  /* We need to update c_n_rows every time we update output_ind because
     of the possibility of early return (in which case we still need to assemble
     the partial results matrix). We could even work with c_n_rows only and ditch
     output_ind, but the inline-c expression is quite verbose, and output_ind is
     more convenient to use in index calculations.
  */
  ($vec-ptr:(int *c_n_rows))[0] = output_ind;
  /* event_ind tracks the current event number */
  int event_ind = 0;

  /* general problem parameters */

  realtype T0 = RCONST(($vec-ptr:(double *c_sol_time))[0]); /* initial time              */
  sunindextype c_dim = $(sunindextype c_dim);           /* number of dependent vars. */

  /* Initialize data structures */

  ARKErrHandlerFn report_error = $fun:(void (*report_error)(int,const char*, const char*, char*, void*));

  /* Initialize odeMaxEventsReached to False */
  ($vec-ptr:(sunindextype *c_diagnostics))[10] = 0;

  y = N_VNew_Serial(c_dim); /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0, report_error)) return 1;
  /* Specify initial condition */
  for (i = 0; i < c_dim; i++) {
    NV_Ith_S(y,i) = ($vec-ptr:(double *c_init_cond))[i];
  };

  // NB: Uses the Newton solver by default
  cvode_mem = CVodeCreate($(int c_method));
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0, report_error)) return(1);
  flag = CVodeInit(cvode_mem, $(int (* c_rhs) (double t, SunVector y[], SunVector dydt[], UserData* params)), T0, y);
  if (check_flag(&flag, "CVodeInit", 1, report_error)) return(1);

  /* Set the error handler */
  flag = CVodeSetErrHandlerFn(cvode_mem, report_error, NULL);
  if (check_flag(&flag, "CVodeSetErrHandlerFn", 1, report_error)) return 1;

  /* Set the user data */
  flag = CVodeSetUserData(cvode_mem, $(UserData* c_rhs_userdata));
  if (check_flag(&flag, "CVodeSetUserData", 1, report_error)) return(1);

  tv = N_VNew_Serial(c_dim); /* Create serial vector for absolute tolerances */
  if (check_flag((void *)tv, "N_VNew_Serial", 0, report_error)) return 1;
  /* Specify tolerances */
  for (i = 0; i < c_dim; i++) {
    NV_Ith_S(tv,i) = ($vec-ptr:(double *c_atol))[i];
  };

  flag = CVodeSetMinStep(cvode_mem, $(double c_minstep));
  if (check_flag(&flag, "CVodeSetMinStep", 1, report_error)) return 1;
  flag = CVodeSetMaxNumSteps(cvode_mem, $(sunindextype c_max_n_steps));
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1, report_error)) return 1;
  flag = CVodeSetMaxErrTestFails(cvode_mem, $(int c_max_err_test_fails));
  if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1, report_error)) return 1;

  /* Specify the scalar relative tolerance and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, $(double c_rtol), tv);
  if (check_flag(&flag, "CVodeSVtolerances", 1, report_error)) return(1);

  /* Specify the root function */
  flag = CVodeRootInit(cvode_mem, $(int c_n_event_specs), $fun:(int (* c_event_fn) (double t, SunVector y[], double gout[], void * params)));
  if (check_flag(&flag, "CVodeRootInit", 1, report_error)) return(1);

  /* Initialize dense matrix data structure and solver */
  A = SUNDenseMatrix(c_dim, c_dim);
  if (check_flag((void *)A, "SUNDenseMatrix", 0, report_error)) return 1;
  LS = SUNDenseLinearSolver(y, A);
  if (check_flag((void *)LS, "SUNDenseLinearSolver", 0, report_error)) return 1;

  /* Attach matrix and linear solver */
  flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
  if (check_flag(&flag, "CVDlsSetLinearSolver", 1, report_error)) return 1;

  /* Set the initial step size if there is one */
  if ($(int c_init_step_size_set)) {
    /* FIXME: We could check if the initial step size is 0 */
    /* or even NaN and then throw an error                 */
    flag = CVodeSetInitStep(cvode_mem, $(double c_init_step_size));
    if (check_flag(&flag, "CVodeSetInitStep", 1, report_error)) return 1;
  }

  /* Set the Jacobian if there is one */
  if ($(int c_jac_set)) {
    flag = CVDlsSetJacFn(cvode_mem, $fun:(int (* c_jac) (double t, SunVector y[], SunVector fy[], SunMatrix Jac[], void * params, SunVector tmp1[], SunVector tmp2[], SunVector tmp3[])));
    if (check_flag(&flag, "CVDlsSetJacFn", 1, report_error)) return 1;
  }

  /* Store initial conditions */
  ($vec-ptr:(double *c_output_mat))[0 * (c_dim + 1) + 0] = ($vec-ptr:(double *c_sol_time))[0];
  for (j = 0; j < c_dim; j++) {
    ($vec-ptr:(double *c_output_mat))[0 * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
  }

  while (1) {
    double ti = ($vec-ptr:(double *c_sol_time))[input_ind];
    flag = CVode(cvode_mem, ti, y, &t, CV_NORMAL); /* call integrator */
    if (check_flag(&flag, "CVode", 1, report_error)) {
      N_Vector ele = N_VNew_Serial(c_dim);
      N_Vector weights = N_VNew_Serial(c_dim);
      flag = CVodeGetEstLocalErrors(cvode_mem, ele);
      // CV_SUCCESS is defined is 0, so we OR the flags
      flag = flag || CVodeGetErrWeights(cvode_mem, weights);
      if (flag == CV_SUCCESS) {
        double *arr_ptr = N_VGetArrayPointer(ele);
        memcpy(($vec-ptr:(double *c_local_error)), arr_ptr, c_dim * sizeof(double));

        arr_ptr = N_VGetArrayPointer(weights);
        memcpy(($vec-ptr:(double *c_var_weight)), arr_ptr, c_dim * sizeof(double));

        ($vec-ptr:(int *c_local_error_set))[0] = 1;
      }
      N_VDestroy(ele);
      N_VDestroy(weights);
      return 1;
    }

    /* Store the results for Haskell */
    ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + 0] = t;
    for (j = 0; j < c_dim; j++) {
      ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
    }
    output_ind++;
    ($vec-ptr:(int *c_n_rows))[0] = output_ind;

    if (flag == CV_ROOT_RETURN) {
      if (event_ind >= $(int c_max_events)) {
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
      int stop_solver = 0;
      flag = CVodeGetRootInfo(cvode_mem, ($vec-ptr:(int *c_root_info)));
      if (check_flag(&flag, "CVodeGetRootInfo", 1, report_error)) return 1;
      for (i = 0; i < $(int c_n_event_specs); i++) {
        int ev = ($vec-ptr:(int *c_root_info))[i];
        int req_dir = ($vec-ptr:(const int *c_requested_event_direction))[i];
        if (ev != 0 && ev * req_dir >= 0) {
          good_event = 1;

          ($vec-ptr:(int *c_actual_event_direction))[event_ind] = ev;
          ($vec-ptr:(int *c_event_index))[event_ind] = i;
          ($vec-ptr:(double *c_event_time))[event_ind] = t;
          event_ind++;
          stop_solver = ($vec-ptr:(int *c_event_stops_solver))[i];

          /* Update the state with the supplied function */
          $fun:(int (* c_apply_event) (int, double, SunVector y[], SunVector z[]))(i, t, y, y);
        }
      }

      if (good_event) {
        ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + 0] = t;
        for (j = 0; j < c_dim; j++) {
          ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
        }
        output_ind++;
        ($vec-ptr:(int *c_n_rows))[0] = output_ind;

        if (stop_solver) {
          break;
        }
        if (event_ind >= $(int c_max_events)) {
          /* We collected the requested number of events. Stop the solver. */
          ($vec-ptr:(sunindextype *c_diagnostics))[10] = 1;
          break;
        }

        flag = CVodeReInit(cvode_mem, t, y);
        if (check_flag(&flag, "CVodeReInit", 1, report_error)) return(1);
      } else {
        /* Since this is not a wanted event, it shouldn't get a row */
        output_ind--;
        ($vec-ptr:(int *c_n_rows))[0] = output_ind;
      }
    }
    else {
      if (++input_ind >= $(int c_n_sol_times))
        break;
    }
  }

  /* The number of actual roots we found */
  ($vec-ptr:(int *c_n_events))[0] = event_ind;

  /* Get some final statistics on how the solve progressed */
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[0] = nst;

  /* FIXME */
  ($vec-ptr:(sunindextype *c_diagnostics))[1] = 0;

  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[2] = nfe;
  /* FIXME */
  ($vec-ptr:(sunindextype *c_diagnostics))[3] = 0;

  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[4] = nsetups;

  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[5] = netf;

  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[6] = nni;

  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[7] = ncfn;

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[8] = ncfn;

  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[9] = ncfn;

  /* Clean up and return */
  N_VDestroy(y);          /* Free y vector          */
  N_VDestroy(tv);         /* Free tv vector         */
  CVodeFree(&cvode_mem);  /* Free integrator memory */
  SUNLinSolFree(LS);      /* Free linear solver     */
  SUNMatDestroy(A);       /* Free A matrix          */

  return CV_SUCCESS;
 } |]
