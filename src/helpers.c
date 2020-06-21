#include <arkode/arkode.h>

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
int check_flag(void *flagvalue, const char *funcname, int opt, ARKErrHandlerFn report_error)
{
  int *errflag;

  /* Check if function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    report_error(0, "hmatrix-sundials", funcname, "returned NULL pointer", NULL);
    return 1;
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;

    if (*errflag < 0) {
      report_error(*errflag, "hmatrix-sundials", funcname, "function failed", NULL);
      return 1;
    }
  }

  return 0;
}
