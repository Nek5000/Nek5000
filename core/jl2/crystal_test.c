#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "name.h"
#include "errmem.h"
#include "types.h"
#include "comm.h"
#include "crystal.h"

int main(int narg, char *arg[])
{
  comm_ext_t world; int np;
  comm_t comm;
  crystal_data crystal;
  uint i,sum, *data, *end;
#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);
#else
  world=0, np=1;
#endif

  comm_init_check(&comm,world,np);
  
  crystal_init(&crystal,&comm);

  crystal.n = (4+(comm.id&1))*comm.np;
  buffer_reserve(&crystal.data,crystal.n*sizeof(uint));
  data = crystal.data.ptr;
  for(i=0;i<comm.np;++i, data+=3+data[2]) {
    data[0] = i, data[1] = comm.id, data[2] = 1;
    data[3] = 2*comm.id;
    if(comm.id&1) data[2] = 2, data[4] = data[3]+1;
  }

#if 0
  data = crystal.data.ptr, end = data + crystal.n;
  for(;data!=end; data+=3+data[2]) {
    uint i;
    printf("%u -> %u:",data[1],data[0]);
    for(i=0;i<data[2];++i) printf(" %u",data[3+i]);
    printf("\n");
  }
#endif
  
  crystal_router(&crystal);

#if 0
  printf("\n");
  data = crystal.data.ptr, end = data + crystal.n;
  for(;data!=end; data+=3+data[2]) {
    uint i;
    printf("%u <- %u:",data[0],data[1]);
    for(i=0;i<data[2];++i) printf(" %u",data[3+i]);
    printf("\n");
  }
#endif
  
  if(crystal.n != comm.np*4 + (comm.np/2)) fail("failure on %u\n",comm.id);
  sum = 0;
  data = crystal.data.ptr, end = data + crystal.n;
  for(;data!=end; data+=3+data[2]) {
    sum+=data[1];
    if(data[3]!=data[1]*2) fail("failure on %u\n",comm.id);
    if(data[1]&1 && (data[2]!=2 || data[4]!=data[3]+1))
      fail("failure on %u\n",comm.id);
  }
  if(sum != comm.np*(comm.np-1)/2) fail("failure on %u\n",comm.id);

  crystal_free(&crystal);
  comm_free(&comm);

  printf("test successful %u/%u\n",(unsigned)comm.id,(unsigned)comm.np);
  
#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
