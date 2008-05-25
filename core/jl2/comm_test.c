#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "name.h"
#include "errmem.h"
#include "types.h"
#include "comm.h"

int main(int narg, char *arg[])
{
  comm_ext_t world; int np;
  comm_t comm;
  ulong sum[2];
#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);
#else
  world=0, np=1;
#endif

  comm_init_check(&comm,world,np);
  
  comm_partial_sum_ul(sum,&comm,comm.id+1);
  printf("%02d: %d %d\n",(int)comm.id,(int)sum[0],(int)sum[1]);

  comm_free(&comm);
  
#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
