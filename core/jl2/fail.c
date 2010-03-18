#include <stdio.h>  /* sprintf, vfprintf, stderr */
#include <stdarg.h> /* va_list, va_start, ... */
#include <stdlib.h> /* exit */
#include <string.h> /* memcpy, and str* functions in comm_fail */
#include "name.h"
#include "fail.h"
#include "types.h"
#include "comm.h"

#define nek_exitt FORTRAN_UNPREFIXED(exitt,EXITT)

void fail(int status, const char *fmt, ...)
{
  int le, lf;
  static char extfmt[1024];
#ifdef MPI
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD,&p);
  sprintf(extfmt, "ERROR (proc %d): ",(int)p);
#else
  strcpy(extfmt, "ERROR: ");
#endif
  le = strlen(extfmt), lf = strlen(fmt);
  if(le+lf+1>1000)
    memcpy(extfmt+le,fmt,lf-(le+lf+1-1000)),
    extfmt[le+lf-(le+lf+1-1000)]='\0';
  else
    strcat(extfmt,fmt);
  strcat(extfmt,"\n");
  {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, extfmt, ap);
    va_end(ap);
  }
#ifdef NO_NEK_EXITT
  exit(status);
#else
  nek_exitt();
#endif  
}
