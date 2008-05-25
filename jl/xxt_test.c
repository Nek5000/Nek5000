
#ifdef MPI

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "types.h"
#include "errmem.h"
#include "crystal.h"
#include "xxt.h"

/*

   +---+---+    1---7---4
   | 0 | 1 |    |   |   |
   +---+---+    3---8---6
   | 0 | 2 |    |   |   |
   +---+---+    2---9---5
   

   element laplacian matrix:
   
    4 -1 -1 -2
   -1  4 -2 -1
   -1 -2  4 -1
   -2 -1 -1  4
   
   assembled matrix:
   

    4    -1          -1 -2
       4 -1             -2 -1
   -1 -1  8          -2 -2 -2
             4    -1 -1 -2
                4 -1    -2 -1
            -1 -1  8 -2 -2 -2
   -1    -2 -1    -2  8 -2
   -2 -2 -2 -2 -2 -2 -2 16 -2
      -1 -2    -1 -2    -2  8
      
*/

const uint nx[3] = {8,4,4};
/*
const ulong x_id[3][8] = { {0,7,3,8, 3,8,2,9},
                           {7,4,8,6},
                           {8,6,9,5} };
*/
const ulong x_id[3][8] = { {0,2,4,5, 4,5,7,8},
                           {2,3,5,6},
                           {5,6,8,9} };

const real bv[3][8][8] = { { {0,1/2.,0,0,0,0,0,0},
                             {0,0,0,0,0,0,0,0},
                             {0,0,1/2.,0,1/2.,0,0,0},
                             {0,0,0,1/4.,0,1/4.,0,0},
                             {0,0,0,0,0,0,0,0},
                             {0,0,0,0,0,0,1,0},
                             {0,0,0,0,0,0,0,1/2.},
                             {0,0,0,0,0,0,0,0} },
                             
                           { {1/2.,0,0,0},
                             {0,1,0,0},
                             {0,0,0,0},
                             {0,0,1/4.,0},
                             {0,0,0,1/2.},
                             {0,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0} },

                           { {0,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0},
                             {1/4.,0,0,0},
                             {0,1/2.,0,0},
                             {0,0,0,0},
                             {0,0,1/2.,0},
                             {0,0,0,1} } };

const uint nz[3] = {32,16,16};
const
uint Ai[3][32] = { {0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3,
                    4,4,4,4, 5,5,5,5, 6,6,6,6, 7,7,7,7},
                   {0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3},
                   {0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3} };
const
uint Aj[3][32] = { {0,1,2,3, 0,1,2,3, 0,1,2,3, 0,1,2,3,
                    4,5,6,7, 4,5,6,7, 4,5,6,7, 4,5,6,7},
                   {0,1,2,3, 0,1,2,3, 0,1,2,3, 0,1,2,3},
                   {0,1,2,3, 0,1,2,3, 0,1,2,3, 0,1,2,3} };
const
real Ar[3][32] = { { 4,-1,-1,-2, -1,4,-2,-1, -1,-2,4,-1, -2,-1,-1,4,
                     4,-1,-1,-2, -1,4,-2,-1, -1,-2,4,-1, -2,-1,-1,4 },
                   { 4,-1,-1,-2, -1,4,-2,-1, -1,-2,4,-1, -2,-1,-1,4 },
                   { 4,-1,-1,-2, -1,4,-2,-1, -1,-2,4,-1, -2,-1,-1,4 } };

int main(int narg, char* arg[])
{
  xxt_data *xxt; crystal_data crystal;
  uint id=0,np=1;
  MPI_Comm comm;
  MPI_Init(&narg,&arg);
  MPI_Comm_dup(MPI_COMM_WORLD,&comm);
  { int i;
    MPI_Comm_rank(comm,&i); id=i;
    MPI_Comm_size(comm,&i); np=i;
  }
  if(np!=3) { puts("run with 3 procs"); exit(1); }

  crystal_init(&crystal,comm);

  xxt = xxt_setup(nx[id], &x_id[id][0],
                  nz[id], &Ai[id][0], &Aj[id][0], &Ar[id][0],
                  0, &crystal);

  crystal_free(&crystal);

  xxt_stats(xxt);
  
#if 1
  { uint i,j; real xv[8];
    for(i=0;i<8;++i) {
      xxt_solve(xv,xxt,&bv[id][i][0]);
      printf("%d col %d:",id,i);
      for(j=0;j<nx[id];++j) printf("\t%.4g",xv[j]);
      printf("\n");
    }
  }
#endif

  xxt_free(xxt);

  MPI_Finalize();
  return 0;
}

#else

int main()
{
  puts("not compiled with MPI support");
  return 1;
}

#endif
