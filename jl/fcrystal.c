/*------------------------------------------------------------------------------
  
  FORTRAN interface for crystal router
  
  integer h, np
  MPI_Comm comm
  call crystal_new(h,comm,np)  ! set h to handle to new instance
  ! it is a runtime error if MPI_Comm_size gives a value different than np
  call crystal_done(h)         ! release instance

  integer*? vi(mi,max)         ! these integer and real types
  integer*? vl(ml,max)         !   better match up with what is
  real      vr(mr,max)         !   in "types.h" 
  call crystal_transfer(h,n,max,vi,mi,vl,ml,vr,mr,p)
  
  - this treats  { vi(:,i), vl(:,i), vr(:,i) } , i in [1 ... n]
      as a list of n tuples with mi integers and md reals each
  - the parameter p indicates that the tuple
      { vi(:,i), vl(:,i), vr(:,i) } should be sent to proc vi(p,i),
      and that on return, vi(p,j) will be the source proc of tuple j
  - n will be set to the number of tuples that came in
  - if more tuples come in than max, n will be set to max+1,
      although only max tuples were stored (the rest are lost)

  ----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#ifdef MPI
#  include <mpi.h>
#endif

#include "fname.h"
#include "errmem.h"
#include "types.h"
#ifdef MPI
#  include "crystal.h"
#  include "tuple_list.h"
#  include "transfer.h"
#else
   typedef void MPI_Comm;
#endif

#define fcrystal_init     FORTRAN_NAME(crystal_new,CRYSTAL_NEW)
#define fcrystal_free     FORTRAN_NAME(crystal_done,CRYSTAL_DONE)
#define fcrystal_transfer FORTRAN_NAME(crystal_transfer,CRYSTAL_TRANSFER)

#ifdef MPI
  static crystal_data **handle=0;
  static int n=0, max=0;
#else
  typedef int MPI_Fint;
#endif

void fcrystal_init(sint *h, const MPI_Fint *comm, const sint *np)
{
#ifdef MPI
  MPI_Comm local_com;
  if(n==max) max+=max/2+1,handle=trealloc(crystal_data*,handle,max);
  handle[n] = tmalloc(crystal_data,1);

  MPI_Comm_dup(MPI_Comm_f2c(*comm),&local_com);

  crystal_init(handle[n],local_com);
  if(*np!=(sint)handle[n]->num)
    fail("crystal_new: passed P=%d, but MPI_Comm_size gives P=%d\n",
         *np,handle[n]->num);
  *h=n++;
#else
  if(*np!=1)
    fail("crystal_new: passed P=%d, but not compiled with -DMPI\n",*np);
  *h=-1;
#endif
}

#ifdef MPI
crystal_data *fcrystal_handle(sint h)
{
  if(h<0 || h>=n || handle[h]==0) failwith("invalid crystal router handle");
  return handle[h];
}
#endif

void fcrystal_free(sint *h)
{
#ifdef MPI
  crystal_data *p = fcrystal_handle(*h);
  handle[*h]=0;
  MPI_Comm_free(&p->comm);
  crystal_free(p);
  free(p);
#endif  
}

void fcrystal_transfer(const sint *h, sint *n, const sint *max,
                       sint  vi[], const sint *mi,
                       slong vl[], const sint *ml,
                       real  vr[], const sint *mr,
                       const sint *pp)
{
#ifdef MPI
  crystal_data *crystal = fcrystal_handle(*h);
  tuple_list tl = { *mi, *ml, *mr, *n, *max, vi, vl, vr };
  const sint p = *pp-1;      /* switch to 0-based index */
  transfer(0,&tl,p,crystal);
  *n = tl.n;
#endif
}

