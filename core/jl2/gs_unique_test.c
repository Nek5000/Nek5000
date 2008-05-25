#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "name.h"
#include "errmem.h"
#include "types.h"
#include "comm.h"
#include "gs.h"

static void test(const comm_t *comm)
{
  uint i,np=comm->np;
  slong *glindex = tmalloc(slong,np*2);
  char *out, *buf = tmalloc(char,80+np*2*30);
  
  for(i=0;i<np;++i) glindex[2*i+1]=glindex[2*i]=i+1;
  
  out = buf+sprintf(buf, "%03d bgn : [", (int)comm->id);
  for(i=0;i<np*2;++i) out += sprintf(out, " %+d", (int)glindex[i]);
  sprintf(out," ]"), puts(buf);
  
  gs_unique(glindex,np*2,comm);

  out = buf+sprintf(buf, "%03d end : [", (int)comm->id);
  for(i=0;i<np*2;++i) out += sprintf(out, " %+d", (int)glindex[i]);
  sprintf(out," ]"), puts(buf);

  free(buf);
  free(glindex);
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
