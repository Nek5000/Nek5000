#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "errmem.h"

void eexit(void) { nek_exitt(); } /* exit wrapper */

void fail(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  eexit();
}

