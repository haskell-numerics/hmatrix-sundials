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
C.include "<sundials/sundials_types.h>"
C.include "<sundials/sundials_math.h>"
C.include "Numeric/Sundials/Foreign_hsc.h"
C.include "../../helpers.h"

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

  ARKRhsFn c_rhs = $(int (*c_rhs)(double, SunVector*, SunVector*, UserData*));
  if ($(int c_method) < MIN_DIRK_NUM) {
    arkode_mem = ARKStepCreate(c_rhs, NULL, T0, y);
  } else {
    arkode_mem = ARKStepCreate(NULL, c_rhs, T0, y);
  }
  if (check_flag(arkode_mem, "ARKStepCreate", 0, report_error)) return 1;

  /* Set the error handler */
  flag = ARKStepSetErrHandlerFn(arkode_mem, report_error, NULL);
  if (check_flag(&flag, "ARKStepSetErrHandlerFn", 1, report_error)) return 1;

  /* Set the user data */
  flag = ARKStepSetUserData(arkode_mem, $(UserData* c_rhs_userdata));
  if (check_flag(&flag, "ARKStepSetUserData", 1, report_error)) return(1);

  tv = N_VNew_Serial(c_dim); /* Create serial vector for absolute tolerances */
  if (check_flag((void *)tv, "N_VNew_Serial", 0, report_error)) return 1;
  /* Specify tolerances */
  for (i = 0; i < c_dim; i++) {
    NV_Ith_S(tv,i) = ($vec-ptr:(double *c_atol))[i];
  };

  flag = ARKStepSetMinStep(arkode_mem, $(double c_minstep));
  if (check_flag(&flag, "ARKStepSetMinStep", 1, report_error)) return 1;
  flag = ARKStepSetMaxNumSteps(arkode_mem, $(sunindextype c_max_n_steps));
  if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1, report_error)) return 1;
  flag = ARKStepSetMaxErrTestFails(arkode_mem, $(int c_max_err_test_fails));
  if (check_flag(&flag, "ARKStepSetMaxErrTestFails", 1, report_error)) return 1;

  /* Specify the scalar relative tolerance and vector absolute tolerances */
  flag = ARKStepSVtolerances(arkode_mem, $(double c_rtol), tv);
  if (check_flag(&flag, "ARKStepSVtolerances", 1, report_error)) return 1;

  /* Specify the root function */
  flag = ARKStepRootInit(arkode_mem, $(int c_n_event_specs), $fun:(int (* c_event_fn) (double t, SunVector y[], double gout[], void * params)));
  if (check_flag(&flag, "ARKStepRootInit", 1, report_error)) return(1);

  /* Initialize dense matrix data structure and solver */
  A = SUNDenseMatrix(c_dim, c_dim);
  if (check_flag((void *)A, "SUNDenseMatrix", 0, report_error)) return 1;
  LS = SUNDenseLinearSolver(y, A);
  if (check_flag((void *)LS, "SUNDenseLinearSolver", 0, report_error)) return 1;

  /* Attach matrix and linear solver */
  flag = ARKStepSetLinearSolver(arkode_mem, LS, A);
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1, report_error)) return 1;

  /* Set the initial step size if there is one */
  if ($(int c_init_step_size_set)) {
    /* FIXME: We could check if the initial step size is 0 */
    /* or even NaN and then throw an error                 */
    flag = ARKStepSetInitStep(arkode_mem, $(double c_init_step_size));
    if (check_flag(&flag, "ARKStepSetInitStep", 1, report_error)) return 1;
  }

  /* Set the Jacobian if there is one */
  if ($(int c_jac_set)) {
    flag = ARKStepSetJacFn(arkode_mem, $fun:(int (* c_jac) (double t, SunVector y[], SunVector fy[], SunMatrix Jac[], void * params, SunVector tmp1[], SunVector tmp2[], SunVector tmp3[])));
    if (check_flag(&flag, "ARKStepSetJacFn", 1, report_error)) return 1;
  }

  /* Store initial conditions */
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
  if (check_flag(&flag, "ARKStepSetTableNum", 1, report_error)) return 1;

  while (1) {

    double ti = ($vec-ptr:(double *c_sol_time))[input_ind];
    flag = ARKStepEvolve(arkode_mem, ti, y, &t, ARK_NORMAL); /* call integrator */
    if (flag == ARK_TOO_CLOSE) {
      /* See Note [CV_TOO_CLOSE]
         No solving was required; just set the time t manually and continue
         as if solving succeeded. */
      t = ti;
    }
    else
    if (check_flag(&flag, "ARKode", 1, report_error)) {
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
      return 1;
    }

    /* Store the results for Haskell */
    ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + 0] = t;
    for (j = 0; j < c_dim; j++) {
      ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
    }
    output_ind++;
    ($vec-ptr:(int *c_n_rows))[0] = output_ind;

    if (flag == ARK_ROOT_RETURN) {
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
      int record_event;
      flag = ARKStepGetRootInfo(arkode_mem, ($vec-ptr:(int *c_root_info)));
      if (check_flag(&flag, "ARKStepGetRootInfo", 1, report_error)) return 1;
      for (i = 0; i < $(int c_n_event_specs); i++) {
        int ev = ($vec-ptr:(int *c_root_info))[i];
        int req_dir = ($vec-ptr:(const int *c_requested_event_direction))[i];
        if (ev != 0 && ev * req_dir >= 0) {
          good_event = 1;
          record_event = ($vec-ptr:(int *c_event_record))[i];
          stop_solver = ($vec-ptr:(int *c_event_stops_solver))[i];

          if (record_event) {
            ($vec-ptr:(int *c_actual_event_direction))[event_ind] = ev;
            ($vec-ptr:(int *c_event_index))[event_ind] = i;
            ($vec-ptr:(double *c_event_time))[event_ind] = t;
            event_ind++;
          }

          /* Update the state with the supplied function */
          $fun:(int (* c_apply_event) (int, double, SunVector y[], SunVector z[]))(i, t, y, y);
        }
      }

      if (good_event) {
        if (record_event) {
          ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + 0] = t;
          for (j = 0; j < c_dim; j++) {
            ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
          }
          output_ind++;
          ($vec-ptr:(int *c_n_rows))[0] = output_ind;
        }

        if (stop_solver) {
          break;
        }
        if (event_ind >= $(int c_max_events)) {
          /* We collected the requested number of events. Stop the solver. */
          ($vec-ptr:(sunindextype *c_diagnostics))[10] = 1;
          break;
        }

        /* Reinitialize */
        if ($(int c_method) < MIN_DIRK_NUM) {
          flag = ARKStepReInit(arkode_mem, c_rhs, NULL, t, y);
        } else {
          flag = ARKStepReInit(arkode_mem, NULL, c_rhs, t, y);
        }
        if (check_flag(&flag, "ARKStepReInit", 1, report_error)) return(1);
      }
      if (!good_event || !record_event) {
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
