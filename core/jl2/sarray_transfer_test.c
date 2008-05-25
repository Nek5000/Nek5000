#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "name.h"
#include "errmem.h"
#include "types.h"
#include "comm.h"
#include "crystal.h"
#include "sarray_transfer.h"

typedef struct {
  uint i;
  ulong l;
  uint p;
  double d;
} rec_t;

int main(int narg, char *arg[])
{
  comm_ext_t world; int np;
  comm_t comm;
  crystal_data crystal;
  array A; rec_t *row;
  uint i;
#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);
#else
  world=0, np=1;
#endif

  comm_init_check(&comm,world,np);
  crystal_init(&crystal,&comm);

  array_init(rec_t,&A,np*3), A.n=np*3, row=A.ptr;
  for(i=0;i<A.n;++i) {
    row[i].i = rand();
    row[i].l = rand();
    row[i].p = rand()%np;
    row[i].d = rand()/(double)rand();
    printf("%d send: %x %x %d %g\n",
      (int)comm.id,(int)row[i].i,(int)row[i].l,(int)row[i].p,row[i].d);
  }
  
  sarray_transfer(rec_t,&A, p, &crystal);

  row=A.ptr;
  for(i=0;i<A.n;++i)
    printf("%d recv: %x %x %d %g\n",
      (int)comm.id,(int)row[i].i,(int)row[i].l,(int)row[i].p,row[i].d);

  array_free(&A);
  crystal_free(&crystal);
  comm_free(&comm);
  
#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
