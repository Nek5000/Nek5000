#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif

#include "errmem.h"
#include "types.h"
#include "poly.h"
#include "tuple_list.h"
#include "crystal.h"
#include "pfindpt.h"

#define N 16
#define NNN (N*N*N)

real lobz[N];
unsigned ord[3]={N,N,N};
real elx[NNN], ely[NNN], elz[NNN];
real *const xw[3] = {elx,ely,elz};

pfindpt_data *p;

int main(int narg, char* arg[])
{
  int id=0,np=1,i,j,k;
#ifndef MPI
  void *crystal=0;
#else
  crystal_data cd, *crystal=&cd;
  MPI_Comm comm;
  MPI_Init(&narg,&arg);
  MPI_Comm_dup(MPI_COMM_WORLD,&comm);
  MPI_Comm_rank(comm,&id);
  MPI_Comm_size(comm,&np);
  crystal_init(crystal,comm);
#endif

  lobatto_nodes(lobz,N);
  for(k=0;k<N;++k) {
    for(j=0;j<N;++j) {
      for(i=0;i<N;++i) {
        unsigned in = (k*N+j)*N+i;
        real theta = (lobz[i]+1)*(PI/12);
        real r = (lobz[j]+1)/2 + (id+1);
        xw[0][in] = r*cosr(theta);
        xw[1][in] = r*sinr(theta);
        xw[2][in] = lobz[k];
      }
    }
  }

  p = pfindpt_setup(3,(const real *const*)xw,ord,1,10,0.01,crystal);

  {
#   define NPT 6
    real pd[NPT][8];
    sint pi[NPT][3];
    tuple_list list = { 3,0,8, NPT,NPT, &pi[0][0],0,&pd[0][0] };
    real     r[NPT] = { 3.5    ,    1.  , .8     , 2.5    , 3.     , 2. };
    real theta[NPT] = { .1*PI/6, .1*PI/6, .1*PI/6, .1*PI/6, .1*PI/6, .1*PI/6};
    real     z[NPT] = { 0      , 0      , 0      , 0      , 0      , 0  };
    for(i=0;i<NPT;++i) {
      pd[i][1]=r[i]*cosr(theta[i]);
      pd[i][2]=r[i]*sinr(theta[i]);
      pd[i][3]=z[i];
    }
    pfindpt(p,&list,0);
    for(i=0;i<(int)list.n;++i) {
      if(pi[i][2]!=-1)
      printf("%d: code=%2d p=%u el=%u r=(%6g,%6g,%6g) x=(%g,%g,%g) dist=%6g\n",
             id,(int)pi[i][2],(int)pi[i][0],(int)pi[i][1],
             (double)pd[i][4],(double)pd[i][5],(double)pd[i][6],
             (double)pd[i][1],(double)pd[i][2],(double)pd[i][3],
             (double)pd[i][0]);
      else
        printf("%d: code=-1 x=(%g,%g,%g)\n",id,
             (double)pd[i][1],(double)pd[i][2],(double)pd[i][3]);
    }
  }
  pfindpt_free(p);
#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}
