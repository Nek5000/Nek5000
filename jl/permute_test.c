#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include "errmem.h"
#include "types.h"
#include "sort.h"
#include "tuple_list.h"

const unsigned mi=2, ml=2, mr=2;

void test1()
{
  buffer buf;
  uint i,j;
  tuple_list tl;
  buffer_init(&buf,1024);
  tuple_list_init_max(&tl,mi,ml,mr,500);
  tl.n=tl.max;
  for(i=0;i<tl.n;++i) {
    int num1 = rand();
    slong num2 = rand();
    num2<<=(CHAR_BIT)*sizeof(int)-1;
    num2|=rand();
    num2<<=(CHAR_BIT)*sizeof(int)-1;
    num2|=rand();
    for(j=0;j<mi;++j) tl.vi[i*mi+j]=num1;
    for(j=0;j<ml;++j) tl.vl[i*ml+j]=num2;
    for(j=0;j<mr;++j) tl.vr[i*mr+j]=num1;
  }
  tuple_list_sort(&tl,0,&buf);
  for(i=0;i<tl.n;++i) {
    for(j=0;j<mi;++j) printf(" %016llx",(long long)tl.vi[i*mi+j]);
    for(j=0;j<ml;++j) printf(" %016llx",(long long)tl.vl[i*ml+j]);
    for(j=0;j<mr;++j) printf(" %g"     ,(double)   tl.vr[i*mr+j]);
    printf("\n");
  }
  printf("on the long:\n");
  tuple_list_sort(&tl,mi,&buf);
  for(i=0;i<tl.n;++i) {
    for(j=0;j<mi;++j) printf(" %016llx",(long long)tl.vi[i*mi+j]);
    for(j=0;j<ml;++j) printf(" %016llx",(long long)tl.vl[i*ml+j]);
    for(j=0;j<mr;++j) printf(" %g"     ,(double)   tl.vr[i*mr+j]);
    printf("\n");
  }
  buffer_free(&buf);
  tuple_list_free(&tl);
}

int main()
{
  test1();
  return 0;
}

