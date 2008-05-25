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

static void test(const comm_t *comm)
{
  gs_data *gs;
  const uint np = comm->np;
  slong *id = tmalloc(slong,np+4);
  T *v = tmalloc(T,np+4);
  uint i;
  id[0] = -(slong)(np+10+3*comm->id);
  for(i=0;i<np;++i) id[i+1] = -(sint)(i+1);
  id[np+1] = comm->id+1;
  id[np+2] = comm->id+1;
  id[np+3] = np-comm->id;
  gs = gs_setup(id,np+4,comm);
  free(id);
  
  for(i=0;i<np+4;++i) v[i] = 1;
  gs_op(v,gs_dom,gs_add,0,gs,0);
  if(comm->id==0) for(i=0;i<np+4;++i) printf("%g\n",v[i]);
  if(comm->id==0) printf("\n");
  for(i=0;i<np+4;++i) v[i] = 1;
  gs_op(v,gs_dom,gs_add,1,gs,0);
  if(comm->id==0) for(i=0;i<np+4;++i) printf("%g\n",v[i]);

  gs_free(gs);
  free(v);
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
