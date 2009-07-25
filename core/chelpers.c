#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include "jl/fname.h"



#define print_stack FORTRAN_NAME(print_stack, PRINT_STACK)

/* Obtain a backtrace and print it to stdout. */
void print_stack(void)
{
  void *bt[50];
  int i;
     
  int bt_size = backtrace(bt, 50);
  char **symbols = backtrace_symbols(bt, bt_size);
     
  printf ("backtrace(): obtained %zd stack frames.\n", bt_size);
  for (i=0; i<bt_size; i++) printf("%s\n", symbols[i]);
  free (symbols);
}
