#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <string.h>
#include "name.h"

#define print_stack FORTRAN_UNPREFIXED(print_stack, PRINT_STACK)
#define sizeOfLongInt FORTRAN_UNPREFIXED(sizeoflongint, SIZEOFLONGINT)
#define getmaxrss FORTRAN_UNPREFIXED(getmaxrss, GETMAXRSS)
#define set_stdout FORTRAN_UNPREFIXED(set_stdout, SET_STDOUT)

#if defined __GLIBC__

#include <execinfo.h>
/* Obtain a backtrace and print it to stdout. */
void print_stack(void)
{
  void *bt[50];
  int i;
  int bt_size = backtrace(bt, 50);
  char **symbols = backtrace_symbols(bt, bt_size);
     
  printf ("backtrace(): obtained %d stack frames.\n", bt_size);
  for (i=0; i<bt_size; i++) printf("%s\n", symbols[i]);
  free (symbols);
}
#else
void print_stack(){};
#endif

double getmaxrss()
{
  struct rusage r_usage;

  getrusage(RUSAGE_SELF,&r_usage);
#if defined(__APPLE__) && defined(__MACH__)
    return (double)r_usage.ru_maxrss;
#else
    return (double)(r_usage.ru_maxrss * 1024L);
#endif
}

int sizeOfLongInt()
{
  return sizeof(long int);
}

void set_stdout(char *f, int flen)
{
  char* logfile;
  int i;

  logfile = (char *) malloc((flen+1)*sizeof(char));
  strncpy(logfile, f, flen);

  for (i=flen-1; i>=0; i--) if (logfile[i] != ' ') break;
  logfile[i+1] = '\0';

  if (flen != 0) {
    printf("redirecting stdout to %s\n",logfile);
    freopen(logfile, "a+", stdout);
  } 
  else {
    const char* s = getenv("NEK_LOGFILE");
    if (s != NULL ) {
      printf("redirecting stdout to %s\n",s);
      freopen(s, "w+", stdout);
    }
  }
}
