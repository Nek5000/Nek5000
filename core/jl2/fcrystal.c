#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "comm.h"
#include "crystal.h"
#include "sort.h"
#include "sarray_sort.h"
#include "sarray_transfer.h"

/*--------------------------------------------------------------------------

  FORTRAN Interface to crystal router
   
  integer h, np
  MPI_Comm comm
  call crystal_setup(h,comm,np)  ! set h to handle to new instance
  ! it is a runtime error if MPI_Comm_size gives a value different than np
  call crystal_free(h)         ! release instance

  integer*? ituple(m,max)   ! integer type matching sint from "types.h"
  call crystal_ituple_transfer(h, ituple,m,n,max, kp)
    - moves each column ituple(:,i), 1 <= i <= n,
      to proc ituple(kp,i)
    - sets n to the number of columns received,
      which may be larger than max (indicating loss of n-max columns)
    - also sets ituple(kp,i) to the source proc of column ituple(:,i)

  call crystal_ituple_sort(h, ituple,m,n, key,nkey)
    - locally sorts columns ituple(:,1...n) in ascending order,
      ranked by ituple(key(1),i),
           then ituple(key(2),i),
           ...
           then ituple(key(nkey),i)
    - no communication; h used for scratch area
    - linear time
    - assumes nonnegative integers

  --------------------------------------------------------------------------*/

#undef   crystal_free
#define ccrystal_free  PREFIXED_NAME(crystal_free)

#define fcrystal_setup           \
  FORTRAN_NAME(crystal_setup          ,CRYSTAL_SETUP          )
#define fcrystal_ituple_sort     \
  FORTRAN_NAME(crystal_ituple_sort    ,CRYSTAL_ITUPLE_SORT    )
#define fcrystal_ituple_transfer \
  FORTRAN_NAME(crystal_ituple_transfer,CRYSTAL_ITUPLE_TRANSFER)
#define fcrystal_free            \
  FORTRAN_NAME(crystal_free           ,CRYSTAL_FREE           )

static crystal_data **handle_array = 0;
static int handle_max = 0;
static int handle_n = 0;

void fcrystal_setup(sint *handle, const MPI_Fint *comm, const sint *np)
{
  crystal_data *p;
  if(handle_n==handle_max)
    handle_max+=handle_max/2+1,
    handle_array=trealloc(crystal_data*,handle_array,handle_max);
  handle_array[handle_n]=p=tmalloc(crystal_data,1);
  comm_init_check(&p->comm, *comm, *np);
  buffer_init(&p->data,1000);
  buffer_init(&p->work,1000);
  p->n=0;
  *handle = handle_n++;
}

void fcrystal_ituple_sort(const sint *handle,
                          sint A[], const sint *m, const sint *n,
                          const sint keys[], const sint *nkey)
{
  const size_t size = (*m)*sizeof(sint);
  sint nk = *nkey;
  buffer *buf;
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle])
    fail(1,"invalid handle to crystal_ituple_sort");
  buf = &handle_array[*handle]->data;
  if(--nk>=0) {
    sortp(buf,0, (uint*)&A[keys[nk]-1],*n,size);
    while(--nk>=0)
      sortp(buf,1, (uint*)&A[keys[nk]-1],*n,size);
    sarray_permute_(ALIGNOF(sint),size,A,*n, buf);
  }
}

void fcrystal_ituple_transfer(const sint *handle,
                              sint A[], const sint *m, sint *n,
                              const sint *nmax, const sint *proc_key)
{
  array ar;
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle])
    fail(1,"invalid handle to crystal_ituple_transfer");
  ar.ptr=A, ar.n=*n, ar.max=*nmax;
  sarray_transfer_(&ar,(*m)*sizeof(sint),(*proc_key-1)*sizeof(sint),
                   1,handle_array[*handle]);
  *n=ar.n;
}

void fcrystal_free(sint *handle)
{
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle])
    fail(1,"invalid handle to crystal_free");
  ccrystal_free(handle_array[*handle]);
  free(handle_array[*handle]);
  handle_array[*handle] = 0;
}


