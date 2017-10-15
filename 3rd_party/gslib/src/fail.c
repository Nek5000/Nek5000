#include <stdio.h>  /* sprintf, vfprintf, stdout */
#include <stdarg.h> /* va_list, va_start, ... */
#include <stdlib.h> /* exit */
#include <string.h> /* memcpy, and str* functions in comm_fail */
#include "name.h"
#include "fail.h"
#include "types.h"
#include "comm.h"

#ifdef USE_USR_EXIT
#define userExitHandler FORTRAN_NAME(userexithandler,USEREXITHANDLER)
#define USEREXIT 1
extern void userExitHandler(int status);
#else
#define USEREXIT 0
void userExitHandler(int status) {};
#endif

void die(int status)
{
  if (USEREXIT) {
  	userExitHandler(status);
    	while(1);
  } else {
    	exit(status); 
    	while(1);
  }
}

void vdiagnostic(const char *prefix, const char *file, unsigned line,
                 const char *fmt, va_list ap)
{
  static char buf[2048]; int n,na,i=0;
  sprintf(buf,"%s(proc %04d, %s:%d): ",prefix,(int)comm_gbl_id,file,line);
  vsprintf(buf+strlen(buf),fmt,ap);
  strcat(buf,"\n");
  n=strlen(buf);
  while(n && (na=fwrite(buf+i,1,n,stdout))) n-=na, i+=na;
  fflush(stdout);
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
  die(status);
}

void fail(int status, const char *file, unsigned line,
          const char *fmt, ...)
{
  va_list ap; va_start(ap,fmt);
  vfail(status,file,line,fmt,ap);
  va_end(ap);
}
