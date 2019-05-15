#include <stdio.h>
#include <math.h>
#include <arkode/arkode.h>                 /* prototypes for ARKODE fcts., consts. */
#include <nvector/nvector_serial.h>        /* serial N_Vector types, fcts., macros */
#include <sunmatrix/sunmatrix_dense.h>     /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>     /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>       /* definition of type realtype          */
#include <sundials/sundials_math.h>

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer  
*/
int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return 1; }

  return 0;
}
