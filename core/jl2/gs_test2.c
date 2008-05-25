#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "name.h"
#include "errmem.h"
#include "types.h"
#include "comm.h"
#include "gs.h"

typedef double T;
const gs_dom_t gs_dom = gs_double;

#define N 4

static void test(const comm_t *comm)
{
  uint i;
  slong id[3][N] = { { 1, 2, 3, 0 },
                     { 1, 0, 4, 5 },
                     { 1, 2, 3, 4 } };
  double v[N] = { .5, .5, .5 };
  
  gs_data *gs;

  gs = gs_setup(id[comm->id],N,comm);

  for(i=0;i<N;++i) v[i] = id[comm->id][i];
  printf("%d initial: [",comm->id);
  for(i=0;i<N;++i) printf(" %g",v[i]);
  printf(" ]\n");
  gs_op(v,gs_dom,gs_add,0,gs,0);
  printf("%d normal: [",comm->id);
  for(i=0;i<N;++i) printf(" %g",v[i]);
  printf(" ]\n");
  for(i=0;i<N;++i) v[i] = id[comm->id][i];
  gs_op(v,gs_dom,gs_add,1,gs,0);
  printf("%d transpose: [",comm->id);
  for(i=0;i<N;++i) printf(" %g",v[i]);
  printf(" ]\n");

  gs_free(gs);
}

int main(int narg, char *arg[])
{
  comm_ext_t world; int np;
  comm_t comm;
  
#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);
#else
  world=0, np=1;
#endif

  comm_init_check(&comm,world,np);

  test(&comm);
  
  comm_free(&comm);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
