#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "types.h"
#include "sort.h"

#define SMALL 22
#define NUM   500
#define SI 9

uint idx[NUM];

#ifdef GLOBAL_INT
  ulong A[NUM][SI], As[NUM];
  uint  B[NUM][SI], Bs[NUM];
  sort_data_long work[2*NUM];
  sort_data workb[2*NUM];
#else
  uint  A[NUM][SI], As[NUM];
  sort_data work[2*NUM];
#endif

int main()
{
  uint i;
  printf("value: %d bytes    index: %d bytes    data: %d bytes\n",
         sizeof(A[0][0]), sizeof(uint), sizeof(work[0]));
  printf("\nsource:\n");
  for(i=0;i!=NUM;++i) {
    A[i][0]=rand();
    A[i][0]<<=CHAR_BIT*sizeof(int)-1;
    A[i][0]^=rand();
    A[i][0]<<=CHAR_BIT*sizeof(int)-1;
    A[i][0]^=rand();
    /* A[i][0]&=0xff0000ff00ULL; */
#ifdef GLOBAL_INT
    B[i][0]=A[i][0];
    printf("%016llx\t%016llx\n",(long long)A[i][0],(long long)B[i][0]);
#else
    printf("%016llx\n",(long long)A[i][0]);
#endif
  }
  printf("\n");
  printf("merge sort:\n");
#ifdef GLOBAL_INT
  sort_long      (&A[0][0],SMALL,SI,&As [0],(void*)work);
  index_sort_long(&A[0][0],SMALL,SI,&idx[0],work);
#else  
  sort           (&A[0][0],SMALL,SI,&As [0],(void*)work);
  index_sort     (&A[0][0],SMALL,SI,&idx[0],work);
#endif
  for(i=0;i!=SMALL;++i)
    printf("%u\t%016llx\t%016llx\n",(unsigned)idx[i],
           (long long)A[idx[i]][0],(long long)As[i]);
  printf("\n");
  printf("radix sort:\n");
#ifdef GLOBAL_INT
  sort_long      (&A[0][0],NUM,SI,&As [0],(void*)work);
  index_sort_long(&A[0][0],NUM,SI,&idx[0],work);
#else  
  sort           (&A[0][0],NUM,SI,&As [0],(void*)work);
  index_sort     (&A[0][0],NUM,SI,&idx[0],work);
#endif
  for(i=0;i!=NUM;++i)
    printf("%u\t%016llx\t%016llx\n",(unsigned)idx[i],
           (long long)A[idx[i]][0],(long long)As[i]);
#ifdef GLOBAL_INT
  printf("\nsmall integers:\n");
  printf("\n");
  printf("merge sort:\n");
  sort      (&B[0][0],SMALL,SI,&Bs [0],(void*)workb);
  index_sort(&B[0][0],SMALL,SI,&idx[0],workb);
  for(i=0;i!=SMALL;++i)
    printf("%u\t%016llx\t%016llx\n",(unsigned)idx[i],
           (long long)B[idx[i]][0],(long long)Bs[i]);
  printf("\n");
  printf("radix sort:\n");
  sort      (&B[0][0],NUM,SI,&Bs [0],(void*)workb);
  index_sort(&B[0][0],NUM,SI,&idx[0],workb);
  for(i=0;i!=NUM;++i)
    printf("%u\t%016llx\t%016llx\n",(unsigned)idx[i],
           (long long)B[idx[i]][0],(long long)Bs[i]);
#endif    
  return 0;
}

