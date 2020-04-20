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
  int retval = CV_SUCCESS;

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

  /* t_start tracks the starting point of the integration in order to detect
     empty integration interval and avoid a potential infinite loop;
     see Note [CV_TOO_CLOSE]. Unlike T0, t_start is updated every time we
     restart the solving after handling (or not) an event.
     Why not just look for the last recorded time in c_output_mat? Because
     an event may have eventRecord = False and not be present there.
  */
  double t_start = T0;

  /* Initialize data structures */

  ARKErrHandlerFn report_error = $fun:(void (*report_error)(int,const char*, const char*, char*, void*));

  if ($(double c_fixedstep) > 0.0) {
    report_error(0, "hmatrix-sundials", "solveC", "fixedStep cannot be used with CVode", NULL);
    return 6426;
  }

  /* Initialize odeMaxEventsReached to False */
  ($vec-ptr:(sunindextype *c_diagnostics))[10] = 0;

  y = N_VNew_Serial(c_dim); /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0, report_error)) return 6896;
  /* Specify initial condition */
  for (i = 0; i < c_dim; i++) {
    NV_Ith_S(y,i) = ($vec-ptr:(double *c_init_cond))[i];
  };

  // NB: Uses the Newton solver by default
  cvode_mem = CVodeCreate($(int c_method));
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0, report_error)) return(8396);
  flag = CVodeInit(cvode_mem, $(int (* c_rhs) (double t, SunVector y[], SunVector dydt[], UserData* params)), T0, y);
  if (check_flag(&flag, "CVodeInit", 1, report_error)) return(1960);

  /* Set the error handler */
  flag = CVodeSetErrHandlerFn(cvode_mem, report_error, NULL);
  if (check_flag(&flag, "CVodeSetErrHandlerFn", 1, report_error)) return 1093;

  /* Set the user data */
  flag = CVodeSetUserData(cvode_mem, $(UserData* c_rhs_userdata));
  if (check_flag(&flag, "CVodeSetUserData", 1, report_error)) return(1949);

  tv = N_VNew_Serial(c_dim); /* Create serial vector for absolute tolerances */
  if (check_flag((void *)tv, "N_VNew_Serial", 0, report_error)) return 6471;
  /* Specify tolerances */
  for (i = 0; i < c_dim; i++) {
    NV_Ith_S(tv,i) = ($vec-ptr:(double *c_atol))[i];
  };

  flag = CVodeSetMinStep(cvode_mem, $(double c_minstep));
  if (check_flag(&flag, "CVodeSetMinStep", 1, report_error)) return 6433;
  flag = CVodeSetMaxNumSteps(cvode_mem, $(sunindextype c_max_n_steps));
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1, report_error)) return 9904;
  flag = CVodeSetMaxErrTestFails(cvode_mem, $(int c_max_err_test_fails));
  if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1, report_error)) return 2512;

  /* Specify the scalar relative tolerance and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, $(double c_rtol), tv);
  if (check_flag(&flag, "CVodeSVtolerances", 1, report_error)) return(6212);

  /* Specify the root function */
  flag = CVodeRootInit(cvode_mem, $(int c_n_event_specs), $fun:(int (* c_event_fn) (double t, SunVector y[], double gout[], void * params)));
  if (check_flag(&flag, "CVodeRootInit", 1, report_error)) return(6290);

  /* Initialize dense matrix data structure and solver */
  A = SUNDenseMatrix(c_dim, c_dim);
  if (check_flag((void *)A, "SUNDenseMatrix", 0, report_error)) return 9061;
  LS = SUNDenseLinearSolver(y, A);
  if (check_flag((void *)LS, "SUNDenseLinearSolver", 0, report_error)) return 9316;

  /* Attach matrix and linear solver */
  flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
  if (check_flag(&flag, "CVDlsSetLinearSolver", 1, report_error)) return 2625;

  /* Set the initial step size if there is one */
  if ($(int c_init_step_size_set)) {
    /* FIXME: We could check if the initial step size is 0 */
    /* or even NaN and then throw an error                 */
    flag = CVodeSetInitStep(cvode_mem, $(double c_init_step_size));
    if (check_flag(&flag, "CVodeSetInitStep", 1, report_error)) return 4010;
  }

  /* Set the Jacobian if there is one */
  if ($(int c_jac_set)) {
    flag = CVDlsSetJacFn(cvode_mem, $fun:(int (* c_jac) (double t, SunVector y[], SunVector fy[], SunMatrix Jac[], void * params, SunVector tmp1[], SunVector tmp2[], SunVector tmp3[])));
    if (check_flag(&flag, "CVDlsSetJacFn", 1, report_error)) return 3124;
  }

  /* Store initial conditions */
  ($vec-ptr:(double *c_output_mat))[0 * (c_dim + 1) + 0] = ($vec-ptr:(double *c_sol_time))[0];
  for (j = 0; j < c_dim; j++) {
    ($vec-ptr:(double *c_output_mat))[0 * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
  }

  while (1) {
    double ti = ($vec-ptr:(double *c_sol_time))[input_ind];
    double next_time_event = ($fun:(double (*c_next_time_event)()))();
    if (next_time_event <= t_start) {
      report_error(0, "hmatrix-sundials", "solveC", "time-based event is in the past", NULL);
      retval = 1;
      goto finish;
    }
    double next_stop_time = fmin(ti, next_time_event);
    flag = CVode(cvode_mem, next_stop_time, y, &t, CV_NORMAL); /* call integrator */
    if (flag == CV_TOO_CLOSE) {
      /* See Note [CV_TOO_CLOSE]
         No solving was required; just set the time t manually and continue
         as if solving succeeded. */
      t = next_stop_time;
    }
    else
    if (t == next_stop_time && t == t_start && flag == CV_ROOT_RETURN) {
      /* See Note [CV_TOO_CLOSE]
         Probably the initial step size was set, and that's why we didn't
         get CV_TOO_CLOSE.
         Pretend that the root didn't happen, lest we keep handling it
         forever. */
      flag = CV_SUCCESS;
    }
    else
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
      return 45;
    }
    int root_based_event = flag == CV_ROOT_RETURN;
    int time_based_event = t == next_time_event;

    /* Store the results for Haskell */
    ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + 0] = t;
    for (j = 0; j < c_dim; j++) {
      ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
    }
    output_ind++;
    ($vec-ptr:(int *c_n_rows))[0] = output_ind;

    if (root_based_event || time_based_event) {
      if (event_ind >= $(int c_max_events)) {
        /* We reached the maximum number of events.
           Either the maximum number of events is set to 0,
           or there's a bug in our code below. In any case return an error.
        */
        return 8630;
      }

      /* How many events triggered? */
      int n_events_triggered = 0;
      int *c_root_info = ($vec-ptr:(int *c_root_info));
      if (root_based_event) {
        flag = CVodeGetRootInfo(cvode_mem, c_root_info);
        if (check_flag(&flag, "CVodeGetRootInfo", 1, report_error)) return 2829;
        for (i = 0; i < $(int c_n_event_specs); i++) {
          int ev = c_root_info[i];
          int req_dir = ($vec-ptr:(const int *c_requested_event_direction))[i];
          if (ev != 0 && ev * req_dir >= 0) {
            /* After the above call to CVodeGetRootInfo, c_root_info has an
            entry per EventSpec. Here we reuse the same array but convert it
            into one that contains indices of triggered events. */
            c_root_info[n_events_triggered++] = i;
          }
        }
      }

      /* Should we stop the solver? */
      int stop_solver = 0;
      /* Should we record the state before/after the event in the output matrix? */
      int record_events = 0;
      if (n_events_triggered > 0 || time_based_event) {
        /* Update the state with the supplied function */
        $fun:(int (* c_apply_event) (int, int*, double, SunVector y[], SunVector z[], int*, int*))(n_events_triggered, c_root_info, t, y, y, &stop_solver, &record_events);
      }

      if (record_events) {
        ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + 0] = t;
        for (j = 0; j < c_dim; j++) {
          ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
        }
        event_ind++;
        output_ind++;
        ($vec-ptr:(int *c_n_rows))[0] = output_ind;
      } else {
        /* Remove the saved row */
        output_ind--;
        ($vec-ptr:(int *c_n_rows))[0] = output_ind;
      }
      if (event_ind >= $(int c_max_events)) {
        ($vec-ptr:(sunindextype *c_diagnostics))[10] = 1;
        stop_solver = 1;
      }
      if (stop_solver) {
        goto finish;
      }

      t_start = t;
      if (n_events_triggered > 0 || time_based_event) {
        flag = CVodeReInit(cvode_mem, t, y);
        if (check_flag(&flag, "CVodeReInit", 1, report_error)) return(1576);
      }
    }
    else {
      if (++input_ind >= $(int c_n_sol_times))
        goto finish;
    }
  }

  finish:

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

{- Note [CV_TOO_CLOSE]
   ~~~~~~~~~~~~~~~~~~~
   One edge condition that may occur is that an event time may exactly
   coincide with a solving time (e.g. they are both exactly equal to an
   integer). Then the following will happen:

   * Sundials will indicate a root at t1.
   * We will handle the event and re-initialize the system at t1.
   * We restart Sundials with the tout being equal to the next solving time,
     which also happens to be equal t1.
   * Sundials sees that the start and end solving times are equal, and
     returns the CV_TOO_CLOSE error.

   Calculating on our side when the start and end times are "too close" by
   Sundials standards is a bit complicated (see the code at the beginning
   of the cvHin function). It's much easier just to call Sundials and
   handle the error.

   For that, however, we need to make sure we ignore CV_TOO_CLOSE in our
   error handler so as not to confuse the end users with mysterious error
   messages in the logs.

   That said, we can't always rely on CV_TOO_CLOSE. When the initial step
   size is set, cvHin is not called, and CV_TOO_CLOSE is not triggered.
   Therefore we also add an explicit check to avoid an infinite loop of
   integrating over an empty interval.
-}
