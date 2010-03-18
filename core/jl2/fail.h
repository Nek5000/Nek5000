#ifndef FAIL_H
#define FAIL_H

#if !defined(NAME_H)
#warning "fail.h" requires "name.h"
#endif

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
