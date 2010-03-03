#ifndef FAIL_H
#define FAIL_H

#if !defined(NAME_H)
#warning "fail.h" requires "name.h"
#endif

#define nek_exitt FORTRAN_UNPREFIXED(exitt,EXITT)
#define eexit PREFIXED_NAME(eexit)
void nek_exitt(void);
void eexit(void);

#define fail PREFIXED_NAME(fail)

#ifdef __GNUC__
#  define FAILDEF() \
   void fail(int status, const char *fmt, ...) __attribute__ ((noreturn));
   FAILDEF()
#  undef FAILDEF
#else
   void fail(int status, const char *fmt, ...);
#endif

#endif
