#include <stdio.h>  /* sprintf, vfprintf, stderr */
#include <stdarg.h> /* va_list, va_start, ... */
#include <stdlib.h> /* exit */
#include <string.h> /* memcpy, and str* functions in comm_fail */
#include "name.h"
#include "fail.h"
#include "types.h"
#include "comm.h"

#define nek_exitt FORTRAN_UNPREFIXED(exitt,EXITT)

void vdiagnostic(const char *prefix, const char *file, unsigned line,
                 const char *fmt, va_list ap)
{
  static char buf[2048];
  sprintf(buf,"%s(proc %04d, %s:%d): ",prefix,(int)comm_gbl_id,file,line);
  vsprintf(buf+strlen(buf),fmt,ap);
  strcat(buf,"\n");
  fwrite(buf,1,strlen(buf),stderr);
  fflush(stderr);
}

void diagnostic(const char *prefix, const char *file, unsigned line,
                const char *fmt, ...)
{
  va_list ap; va_start(ap,fmt);
  vdiagnostic(prefix,file,line,fmt,ap);
  va_end(ap);
}

void vfail(int status, const char *file, unsigned line,
           const char *fmt, va_list ap)
{
  vdiagnostic("ERROR ",file,line,fmt,ap);
#ifdef NO_NEK_EXITT
  exit(status);
#else
  nek_exitt();
#endif  
}

void fail(int status, const char *file, unsigned line,
          const char *fmt, ...)
{
  va_list ap; va_start(ap,fmt);
  vfail(status,file,line,fmt,ap);
  va_end(ap);
}
