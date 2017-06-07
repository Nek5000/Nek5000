#ifdef MPI

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include "errmem.h"
#include "types.h"
#include "crystal.h"
#include "tuple_list.h"
#include "transfer.h"

#define MI 3
#define ML 2
#define MR 3
#define PF 0

int main(int narg, char *arg[])
{
  tuple_list tl;
  int id=0,np=1; uint i;
  crystal_data crystal;
  MPI_Comm comm;
  MPI_Init(&narg,&arg);
  MPI_Comm_dup(MPI_COMM_WORLD,&comm);
  MPI_Comm_rank(comm,&id);
  MPI_Comm_size(comm,&np);
  
  srand(id+1);
  tuple_list_init_max(&tl,MI,ML,MR,np*3);
  tl.n = np*3;
  for(i=0;i<tl.n;++i) {
    int j;
    sint *ri = tl.vi+MI*i; slong *rl = tl.vl+ML*i; real *rr = tl.vr+MR*i;
    for(j=0;j<MI;++j) ri[j] = rand();
    ri[PF] %= np;
    for(j=0;j<ML;++j) rl[j] = rand();
    for(j=0;j<MR;++j) rr[j] = rand()/(real)rand();
  }
  for(i=0;i<tl.n;++i) {
    int j;
    sint *ri = tl.vi+MI*i; slong *rl = tl.vl+ML*i; real *rr = tl.vr+MR*i;
    printf("%d send: ", id);
    for(j=0;j<MI;++j) printf("%lld ",(long long)ri[j]);
    for(j=0;j<ML;++j) printf("%lld ",(long long)rl[j]);
    for(j=0;j<MR;++j) printf("%g ",(double)rr[j]);
    printf("\n");
  }

  crystal_init(&crystal,comm);
  transfer(1,&tl,PF,&crystal);
  crystal_free(&crystal);
  if(tl.n>tl.max) tl.n=tl.max, printf("%d lost some\n", id);

  for(i=0;i<tl.n;++i) {
    int j;
    sint *ri = tl.vi+MI*i; slong *rl = tl.vl+ML*i; real *rr = tl.vr+MR*i;
    printf("%d recv: ", id);
    for(j=0;j<MI;++j) printf("%lld ",(long long)ri[j]);
    for(j=0;j<ML;++j) printf("%lld ",(long long)rl[j]);
    for(j=0;j<MR;++j) printf("%g ",(double)rr[j]);
    printf("\n");
  }

  tuple_list_free(&tl);

  MPI_Finalize();
  return 0;
}
#else

#include <stdio.h>
int main()
{
  printf("Not compiled with -DMPI. Test is meaningless.\n");
  return 0;
}

#endif

