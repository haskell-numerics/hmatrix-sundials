-- |
-- Solution of ordinary differential equation (ODE) initial value problems.
-- See <https://computation.llnl.gov/projects/sundials/sundials-software> for more detail.
module Numeric.Sundials.ARKode
  ( ARKMethod(..)
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
C.include "<arkode/arkode_ls.h>"
C.include "<arkode/arkode_arkstep.h>"
C.include "<nvector/nvector_serial.h>"
C.include "<sunmatrix/sunmatrix_dense.h>"
C.include "<sunlinsol/sunlinsol_dense.h>"
C.include "<sunlinsol/sunlinsol_klu.h>"
C.include "<sundials/sundials_types.h>"
C.include "<sundials/sundials_math.h>"
C.include "../../helpers.h"
C.include "Numeric/Sundials/Foreign_hsc.h"

-- | Available methods for ARKode
data ARKMethod = SDIRK_2_1_2
               | BILLINGTON_3_3_2
               | TRBDF2_3_3_2
               | KVAERNO_4_2_3
               | ARK324L2SA_DIRK_4_2_3
               | CASH_5_2_4
               | CASH_5_3_4
               | SDIRK_5_3_4
               | KVAERNO_5_3_4
               | ARK436L2SA_DIRK_6_3_4
               | KVAERNO_7_4_5
               | ARK548L2SA_DIRK_8_4_5
               | HEUN_EULER_2_1_2
               | BOGACKI_SHAMPINE_4_2_3
               | ARK324L2SA_ERK_4_2_3
               | ZONNEVELD_5_3_4
               | ARK436L2SA_ERK_6_3_4
               | SAYFY_ABURUB_6_3_4
               | CASH_KARP_6_4_5
               | FEHLBERG_6_4_5
               | DORMAND_PRINCE_7_4_5
               | ARK548L2SA_ERK_8_4_5
               | VERNER_8_5_6
               | FEHLBERG_13_7_8
  deriving Show

instance Method ARKMethod where
  methodToInt SDIRK_2_1_2             = sDIRK_2_1_2
  methodToInt BILLINGTON_3_3_2        = bILLINGTON_3_3_2
  methodToInt TRBDF2_3_3_2            = tRBDF2_3_3_2
  methodToInt KVAERNO_4_2_3           = kVAERNO_4_2_3
  methodToInt ARK324L2SA_DIRK_4_2_3   = aRK324L2SA_DIRK_4_2_3
  methodToInt CASH_5_2_4              = cASH_5_2_4
  methodToInt CASH_5_3_4              = cASH_5_3_4
  methodToInt SDIRK_5_3_4             = sDIRK_5_3_4
  methodToInt KVAERNO_5_3_4           = kVAERNO_5_3_4
  methodToInt ARK436L2SA_DIRK_6_3_4   = aRK436L2SA_DIRK_6_3_4
  methodToInt KVAERNO_7_4_5           = kVAERNO_7_4_5
  methodToInt ARK548L2SA_DIRK_8_4_5   = aRK548L2SA_DIRK_8_4_5
  methodToInt HEUN_EULER_2_1_2        = hEUN_EULER_2_1_2
  methodToInt BOGACKI_SHAMPINE_4_2_3  = bOGACKI_SHAMPINE_4_2_3
  methodToInt ARK324L2SA_ERK_4_2_3    = aRK324L2SA_ERK_4_2_3
  methodToInt ZONNEVELD_5_3_4         = zONNEVELD_5_3_4
  methodToInt ARK436L2SA_ERK_6_3_4    = aRK436L2SA_ERK_6_3_4
  methodToInt SAYFY_ABURUB_6_3_4      = sAYFY_ABURUB_6_3_4
  methodToInt CASH_KARP_6_4_5         = cASH_KARP_6_4_5
  methodToInt FEHLBERG_6_4_5          = fEHLBERG_6_4_5
  methodToInt DORMAND_PRINCE_7_4_5    = dORMAND_PRINCE_7_4_5
  methodToInt ARK548L2SA_ERK_8_4_5    = aRK548L2SA_ERK_8_4_5
  methodToInt VERNER_8_5_6            = vERNER_8_5_6
  methodToInt FEHLBERG_13_7_8         = fEHLBERG_13_7_8

  methodSolver = ARKode

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
  void *arkode_mem = NULL;   /* empty ARKode memory structure                */
  realtype t;
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;

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

  /* Initialize odeMaxEventsReached to False */
  ($vec-ptr:(sunindextype *c_diagnostics))[10] = 0;

  y = N_VNew_Serial(c_dim); /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0, report_error)) return 6896;
  /* Specify initial condition */
  for (i = 0; i < c_dim; i++) {
    NV_Ith_S(y,i) = ($vec-ptr:(double *c_init_cond))[i];
  };

  ARKRhsFn c_rhs = $(int (*c_rhs)(double, SunVector*, SunVector*, UserData*));
  if ($(int c_method) < MIN_DIRK_NUM) {
    arkode_mem = ARKStepCreate(c_rhs, NULL, T0, y);
  } else {
    arkode_mem = ARKStepCreate(NULL, c_rhs, T0, y);
  }
  if (check_flag(arkode_mem, "ARKStepCreate", 0, report_error)) return 8396;

  /* Set the error handler */
  flag = ARKStepSetErrHandlerFn(arkode_mem, report_error, NULL);
  if (check_flag(&flag, "ARKStepSetErrHandlerFn", 1, report_error)) return 1093;

  double c_fixedstep = $(double c_fixedstep);
  if (c_fixedstep > 0.0) {
    flag = ARKStepSetFixedStep(arkode_mem, c_fixedstep);
    if (check_flag(&flag, "ARKStepSetFixedStep", 1, report_error)) return(1);
  }

  /* Set the user data */
  flag = ARKStepSetUserData(arkode_mem, $(UserData* c_rhs_userdata));
  if (check_flag(&flag, "ARKStepSetUserData", 1, report_error)) return(1949);

  tv = N_VNew_Serial(c_dim); /* Create serial vector for absolute tolerances */
  if (check_flag((void *)tv, "N_VNew_Serial", 0, report_error)) return 6471;
  /* Specify tolerances */
  for (i = 0; i < c_dim; i++) {
    NV_Ith_S(tv,i) = ($vec-ptr:(double *c_atol))[i];
  };

  flag = ARKStepSetMinStep(arkode_mem, $(double c_minstep));
  if (check_flag(&flag, "ARKStepSetMinStep", 1, report_error)) return 6433;
  flag = ARKStepSetMaxNumSteps(arkode_mem, $(sunindextype c_max_n_steps));
  if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1, report_error)) return 9904;
  flag = ARKStepSetMaxErrTestFails(arkode_mem, $(int c_max_err_test_fails));
  if (check_flag(&flag, "ARKStepSetMaxErrTestFails", 1, report_error)) return 2512;

  /* Specify the scalar relative tolerance and vector absolute tolerances */
  flag = ARKStepSVtolerances(arkode_mem, $(double c_rtol), tv);
  if (check_flag(&flag, "ARKStepSVtolerances", 1, report_error)) return(6212);

  /* Specify the root function */
  flag = ARKStepRootInit(arkode_mem, $(int c_n_event_specs), $fun:(int (* c_event_fn) (double t, SunVector y[], double gout[], void * params)));
  if (check_flag(&flag, "ARKStepRootInit", 1, report_error)) return(6290);

  /* Initialize a jacobian matrix and solver */
  int c_sparse_jac = $(int c_sparse_jac);
  if (c_sparse_jac) {
    A = SUNSparseMatrix(c_dim, c_dim, c_sparse_jac, CSC_MAT);
    if (check_flag((void *)A, "SUNSparseMatrix", 0, report_error)) return 9061;
    LS = SUNLinSol_KLU(y, A);
    if (check_flag((void *)LS, "SUNLinSol_KLU", 0, report_error)) return 9316;
  } else {
    A = SUNDenseMatrix(c_dim, c_dim);
    if (check_flag((void *)A, "SUNDenseMatrix", 0, report_error)) return 9061;
    LS = SUNDenseLinearSolver(y, A);
    if (check_flag((void *)LS, "SUNDenseLinearSolver", 0, report_error)) return 9316;
  }

  /* Attach matrix and linear solver */
  flag = ARKStepSetLinearSolver(arkode_mem, LS, A);
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1, report_error)) return 2625;

  /* Set the initial step size if there is one */
  if ($(int c_init_step_size_set)) {
    /* FIXME: We could check if the initial step size is 0 */
    /* or even NaN and then throw an error                 */
    flag = ARKStepSetInitStep(arkode_mem, $(double c_init_step_size));
    if (check_flag(&flag, "ARKStepSetInitStep", 1, report_error)) return 4010;
  }

  /* Set the Jacobian if there is one */
  if ($(int c_jac_set)) {
    flag = ARKStepSetJacFn(arkode_mem, $fun:(int (* c_jac) (double t, SunVector y[], SunVector fy[], SunMatrix Jac[], void * params, SunVector tmp1[], SunVector tmp2[], SunVector tmp3[])));
    if (check_flag(&flag, "ARKStepSetJacFn", 1, report_error)) return 3124;
  }

  /* Store initial conditions */
  ($vec-ptr:(double *c_output_mat))[0 * (c_dim + 1) + 0] = ($vec-ptr:(double *c_sol_time))[0];
  for (j = 0; j < c_dim; j++) {
    ($vec-ptr:(double *c_output_mat))[0 * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
  }

  /* Explicitly set the method */
  if ($(int c_method) >= MIN_DIRK_NUM) {
    /* Implicit */
    flag = ARKStepSetTableNum(arkode_mem, $(int c_method), -1);
  } else {
    /* Explicit */
    flag = ARKStepSetTableNum(arkode_mem, -1, $(int c_method));
  }
  if (check_flag(&flag, "ARKStepSetTableNum", 1, report_error)) return 26643;

  while (1) {
    double ti = ($vec-ptr:(double *c_sol_time))[input_ind];
    double next_time_event = ($fun:(double (*c_next_time_event)()))();
    if (next_time_event < t_start) {
      size_t msg_size = 1000;
      char *msg = alloca(msg_size);
      snprintf(msg, msg_size, "time-based event is in the past: next event time = %.4f while we are at %.4f", next_time_event, t_start);
      report_error(0, "hmatrix-sundials", "solveC", msg, NULL);
      retval = 5669;
      goto finish;
    }
    double next_stop_time = fmin(ti, next_time_event);
    flag = ARKStepEvolve(arkode_mem, next_stop_time, y, &t, ARK_NORMAL); /* call integrator */
    int root_based_event = flag == ARK_ROOT_RETURN;
    int time_based_event = t == next_time_event;
    if (flag == ARK_TOO_CLOSE) {
      /* See Note [CV_TOO_CLOSE]
         No solving was required; just set the time t manually and continue
         as if solving succeeded. */
      t = next_stop_time;
    }
    else
    if (t == next_stop_time && t == t_start && flag == ARK_ROOT_RETURN && !time_based_event) {
      /* See Note [CV_TOO_CLOSE]
         Probably the initial step size was set, and that's why we didn't
         get ARK_TOO_CLOSE.
         Pretend that the root didn't happen, lest we keep handling it
         forever. */
      flag = ARK_SUCCESS;
    }
    else
    if (!(flag == ARK_TOO_CLOSE && time_based_event) &&
      check_flag(&flag, "ARKStepEvolve", 1, report_error)) {

      N_Vector ele = N_VNew_Serial(c_dim);
      N_Vector weights = N_VNew_Serial(c_dim);
      flag = ARKStepGetEstLocalErrors(arkode_mem, ele);
      // ARK_SUCCESS is defined is 0, so we OR the flags
      flag = flag || ARKStepGetErrWeights(arkode_mem, weights);
      if (flag == ARK_SUCCESS) {
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
        flag = ARKStepGetRootInfo(arkode_mem, c_root_info);
        if (check_flag(&flag, "ARKStepGetRootInfo", 1, report_error)) return 2829;
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
        /* A corner case: if the time-based event triggers at the very beginning,
           then we don't want to duplicate the initial row, so rewind it back.
           Note that we do this only in the branch where record_events is true;
           otherwise we may end up erasing the initial row (see below). */
        if (t == ($vec-ptr:(double *c_sol_time))[0] && output_ind == 2) {
          output_ind--;
          /* c_n_rows will be updated below anyway */
        }

        ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + 0] = t;
        for (j = 0; j < c_dim; j++) {
          ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
        }
        event_ind++;
        output_ind++;
        ($vec-ptr:(int *c_n_rows))[0] = output_ind;
      } else {
        /* Remove the saved row — unless the event time also coincides with a requested time point */
        if (t != ti) {
          output_ind--;
          ($vec-ptr:(int *c_n_rows))[0] = output_ind;
        }
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
        /* Reinitialize */
        if ($(int c_method) < MIN_DIRK_NUM) {
          flag = ARKStepReInit(arkode_mem, c_rhs, NULL, t, y);
        } else {
          flag = ARKStepReInit(arkode_mem, NULL, c_rhs, t, y);
        }
        if (check_flag(&flag, "ARKStepReInit", 1, report_error)) return(1576);
      }
    }
    if (t == ti) {
      if (++input_ind >= $(int c_n_sol_times))
        goto finish;
    }
  }

  finish:

  /* The number of actual roots we found */
  ($vec-ptr:(int *c_n_events))[0] = event_ind;

  /* Get some final statistics on how the solve progressed */
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKStepGetNumSteps", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[0] = nst;

  flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKStepGetNumStepAttempts", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[1] = nst_a;

  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKStepGetNumRhsEvals", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[2] = nfe;
  ($vec-ptr:(sunindextype *c_diagnostics))[3] = nfi;

  flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[4] = nsetups;

  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKStepGetNumErrTestFails", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[5] = netf;

  flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[6] = nni;

  flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[7] = ncfn;

  flag = ARKStepGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKStepGetNumJacEvals", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[8] = ncfn;

  flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKStepGetNumRhsEvals", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[9] = ncfn;

  /* Clean up and return */
  N_VDestroy(y);            /* Free y vector          */
  N_VDestroy(tv);           /* Free tv vector         */
  ARKStepFree(&arkode_mem);  /* Free integrator memory */
  SUNLinSolFree(LS);        /* Free linear solver     */
  SUNMatDestroy(A);         /* Free A matrix          */

  return ARK_SUCCESS;
 } |]
