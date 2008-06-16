#include <stdio.h>
#include <limits.h>
#include <stddef.h>
#include <stdlib.h>
#include "name.h"
#include "errmem.h"
#include "types.h"
#include "sarray_sort.h"

int main()
{
  typedef struct { double d; slong l; sint i; float f; } rec_t;
  buffer buf = {0,0,0};
  rec_t recs[500];
  uint i;
  
  for(i=0;i<500;++i) {
    sint num1 = rand() & 0xff;
    slong num2 = rand();
    /*num2<<=(CHAR_BIT)*sizeof(int)-1;
    num2|=rand();
    num2<<=(CHAR_BIT)*sizeof(int)-1;
    num2|=rand();*/
    recs[i].d = num2;
    recs[i].f = num2;
    recs[i].l = num2;
    recs[i].i = num1;
  }
  sarray_sort_two(rec_t,recs,500, i,0, l,1, &buf);
  for(i=0;i<500;++i)
    printf("%g\t%g\t%ld\t%d\n",
      recs[i].d,recs[i].f,(long)recs[i].l,(int)recs[i].i);

  printf("\n");
  sarray_sort(rec_t,recs,500, l,1, &buf);
  for(i=0;i<500;++i)
    printf("%g\t%g\t%ld\t%d\n",
      recs[i].d,recs[i].f,(long)recs[i].l,(int)recs[i].i);
  buffer_free(&buf);
  return 0;
}

