#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#ifdef MPI
#  include <mpi.h>
#endif

#include "types.h"
#include "errmem.h"
#ifdef MPI
#  include "crystal.h"
#endif
#include "xxt.h"

#define M 3

int main(int narg, char* arg[])
{
  uint n; ulong *xid;
  uint nz; uint *Ai, *Aj; real *A;
  uint i;
  real *x, *b, *x2;

  xxt_data *xxt;
  uint id=0,np=1;
#ifndef MPI
  void *crystal=0;
#else
  crystal_data crystal_, *crystal=&crystal_;
  MPI_Comm comm;
  MPI_Init(&narg,&arg);
  MPI_Comm_dup(MPI_COMM_WORLD,&comm);
  { int i;
    MPI_Comm_rank(comm,&i); id=i;
    MPI_Comm_size(comm,&i); np=i;
  }
  crystal_init(crystal,comm);
#endif

  n=M+1; if(id==np-1) --n;
  xid = tmalloc(ulong,n); x=tmalloc(real,3*n), b=x+n, x2=b+n;
  for(i=0;i<n;++i) xid[i]=1+id*M+i;
  nz=2*M; if(id==np-1) --nz;
  Ai=tmalloc(uint,2*nz), Aj=Ai+nz;
  A =tmalloc(real,nz);
  for(i=0;i<M;++i) Ai[i]=i,Aj[i]=i,A[i]=2;
  if(id==0) A[0]=1;
  if(id==np-1) A[n-1]=1;
  for(i=M;i<nz;++i) Ai[i]=i-M,Aj[i]=i-M+1,A[i]=-1;

  xxt = xxt_setup(n,xid, nz,Ai,Aj,A, 1, crystal);
#ifdef MPI
  crystal_free(crystal);
#endif
  xxt_stats(xxt);
  
  {
    real tn = M*np, mean = 1*(tn+1)/2;
    real avg = 1*(tn-1)*(tn+1)/12;
    for(i=0;i<n;++i) x[i]=(xid[i]-mean)*(xid[i]-mean)-avg;
    for(i=0;i<n;++i) b[i]=A[i]*x[i]
                - (i>0   || id!=0    ? (xid[i]-1-mean)*(xid[i]-1-mean)-avg : 0)
                - (i+1<n || id!=np-1 ? (xid[i]+1-mean)*(xid[i]+1-mean)-avg : 0);
    if(id!=np-1) b[n-1]=0;
  }
  xxt_solve(x2,xxt,b);
  xxt_free(xxt);
  
  { real dif=0;
    for(i=0;i<n;++i) {
      real d=fabsr(x[i]-x2[i])/(.1+fabsr(x[i]));
      if(d>dif) dif=d;
    }
    printf("%d : max dif = %g\n",id,dif);
  }

#ifdef MPI
  MPI_Finalize();
#endif

  if(id==0) printf("test successful\n");
  
  return 0;
}

