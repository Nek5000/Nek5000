/* simple stand-alone test for parallel gather-scatter routines
   assumes gather-scatter routines were compiled with default names
   can compile to sequential version if MPI is not defined
   
   the test is as follows, where N is the number of procs:
     there are N physical nodes (vertices)
     each proc has 2 local/virtual nodes mapping to each physical node,
       for a total of 2*N*N virtual nodes
     virtual nodes are given values that correspond to a sequential ordering
       (so that they range from 0 to 2*N*N-1)
       the addition operation is performed and the result is checked,
       the correct result being known a priori
     the addition operation is also checked, in a similar manner, for
       both the cpgs_op_vec and cpgs_op_many routines with vector dimension 3
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef MPI
#  include <mpi.h>
#else
   typedef void MPI_Comm;
#endif
#include "types.h"
#include "fname.h"

#define crystal_new      FORTRAN_NAME(crystal_new,CRYSTAL_NEW)
#define crystal_done     FORTRAN_NAME(crystal_done,CRYSTAL_DONE)
#define crystal_transfer FORTRAN_NAME(crystal_transfer,CRYSTAL_TRANSFER)

void crystal_new(sint *h, const MPI_Comm *comm, const sint *np);
void crystal_done(sint *h);

#define cpgs_setup   FORTRAN_NAME(cpgs_setup  ,CPGS_SETUP  )
#define cpgs_op      FORTRAN_NAME(cpgs_op     ,CPGS_OP     )
#define cpgs_op_vec  FORTRAN_NAME(cpgs_op_vec ,CPGS_OP_VEC )
#define cpgs_op_many FORTRAN_NAME(cpgs_op_many,CPGS_OP_MANY)
#define cpgs_free    FORTRAN_NAME(cpgs_free   ,CPGS_FREE   )

void cpgs_setup(sint *handle, const sint *crystal_handle,
                const slong v[], const sint *vn, const sint *maxv);
void cpgs_op(const sint *handle, real u[], const sint *op);
void cpgs_op_vec(const sint *handle, real u[], const sint *n, const sint *op);
void cpgs_op_many(const sint *handle,
                  real u1[], real u2[], real u3[],
                  real u4[], real u5[], real u6[],
                  const sint *n, const sint *op);
void cpgs_free(sint *handle);

void assert_is_zero(real v)
{
  if(fabs(v) < 1e-20) return;
  printf("test failed\n");
  exit(1);
}

int main(int narg, char* arg[])
{
  sint id=0,np=1;
  sint i,handle,maxv=3;
  sint chandle;
  real *u;
  slong *glindex;
#ifndef MPI
  int comm;
#else
  MPI_Comm comm;
  MPI_Init(&narg,&arg);
  MPI_Comm_dup(MPI_COMM_WORLD,&comm);
  { int i;
    MPI_Comm_rank(comm,&i); id=i;
    MPI_Comm_size(comm,&i); np=i;
  }
#endif
  
  crystal_new(&chandle,&comm,&np);

  glindex = malloc(np*2*sizeof(slong));
  for(i=0;i<np;++i) glindex[2*i+1] = glindex[2*i] = i+1;
  i=np*2;
  cpgs_setup(&handle,&chandle,glindex,&i,&maxv);
  crystal_done(&chandle);
  free(glindex);
  
  u = malloc(np*2*sizeof(real));
  for(i=0;i<np;++i) u[2*i  ] = (real)( 2*np*id + 2*i ),
                    u[2*i+1] = (real)( 2*np*id + 2*i+1 );
  i=1, cpgs_op(&handle,u,&i);
  /* for(i=0;i<np;++i) printf(" (%g %g)", u[2*i], u[2*i+1]); */
  for(i=0;i<np;++i) assert_is_zero( np*(2*np*(np-1)+4*i+1) - u[2*i] ),
                    assert_is_zero( np*(2*np*(np-1)+4*i+1) - u[2*i+1]  );
  free(u);
  
  u = malloc(np*2*3*sizeof(real));
  for(i=0;i<np;++i)
    u[3*(2*i  )+0] = (real)( 3*(2*np*id + 2*i  ) + 0 ),
    u[3*(2*i  )+1] = (real)( 3*(2*np*id + 2*i  ) + 1 ),
    u[3*(2*i  )+2] = (real)( 3*(2*np*id + 2*i  ) + 2 ),
    u[3*(2*i+1)+0] = (real)( 3*(2*np*id + 2*i+1) + 0 ),
    u[3*(2*i+1)+1] = (real)( 3*(2*np*id + 2*i+1) + 1 ),
    u[3*(2*i+1)+2] = (real)( 3*(2*np*id + 2*i+1) + 2 );
  /* for(i=0;i<np;++i) {
    int j;
    printf("%d: ( ", id);
    for(j=3*(2*i);j<=3*(2*i+1)+2;++j) printf("%g ",u[j]);
    printf(")\n");
  } */
  i=1, maxv=3, cpgs_op_vec(&handle,u,&maxv,&i);
  /* for(i=0;i<np;++i) {
    int j;
    printf("%d: ( ", id);
    for(j=3*(2*i);j<=3*(2*i+1)+2;++j) printf("%g ",u[j]);
    printf(")\n");
  } */
  for(i=0;i<np;++i)
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*0) - u[3*(2*i  )+0] ),
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*1) - u[3*(2*i  )+1] ),
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*2) - u[3*(2*i  )+2] ),
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*0) - u[3*(2*i+1)+0] ),
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*1) - u[3*(2*i+1)+1] ),
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*2) - u[3*(2*i+1)+2] );
  free(u);

  u = malloc(np*2*3*sizeof(real));
  for(i=0;i<np;++i)
    u[2*np*0+(2*i  )] = (real)( 3*(2*np*id + 2*i  ) + 0 ),
    u[2*np*1+(2*i  )] = (real)( 3*(2*np*id + 2*i  ) + 1 ),
    u[2*np*2+(2*i  )] = (real)( 3*(2*np*id + 2*i  ) + 2 ),
    u[2*np*0+(2*i+1)] = (real)( 3*(2*np*id + 2*i+1) + 0 ),
    u[2*np*1+(2*i+1)] = (real)( 3*(2*np*id + 2*i+1) + 1 ),
    u[2*np*2+(2*i+1)] = (real)( 3*(2*np*id + 2*i+1) + 2 );
  i=1, maxv=3, cpgs_op_many(&handle,u,u+2*np,u+4*np,0,0,0,&maxv,&i);
  for(i=0;i<np;++i)
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*0) - u[2*np*0+(2*i  )] ),
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*1) - u[2*np*1+(2*i  )] ),
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*2) - u[2*np*2+(2*i  )] ),
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*0) - u[2*np*0+(2*i+1)] ),
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*1) - u[2*np*1+(2*i+1)] ),
    assert_is_zero( np*(6*np*(np-1)+12*i+3+2*2) - u[2*np*2+(2*i+1)] );
  free(u);
  
  cpgs_free(&handle);
  printf("test on node %d/%d succeeded\n", (int)id+1, (int)np);
#ifdef MPI  
  MPI_Finalize();
#endif
  return 0;
}
