#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef PREFIX
#  define PASTE1(a,b) a##b
#  define PASTE2(a,b) PASTE1(a,b)
#  define fail PASTE2(PREFIX,fail)
#endif

void fail(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  exit(1);
}

