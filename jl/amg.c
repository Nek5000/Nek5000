#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef MPI
#  include <mpi.h>
#endif

#include "errmem.h"
#include "types.h"
#include "sort.h"
#ifdef MPI
#  include "crystal.h"
#else
  typedef void crystal_data;
#endif
#include "gs.h"
#include "tuple_list.h"
#include "transfer.h"

/* condensed dof tuple list
   ( index, proc, level,   id ) */
static const unsigned cdof_mi = 3, cdof_ml = 1, cdof_mr = 1;
static const unsigned cdof_index = 0, cdof_proc = 1, cdof_level = 2;

/* sorts an array of ids, removes 0's and duplicates;
   just returns the permutation */
static uint unique_ids(uint n, const ulong *id, sint *perm, buffer *buf)
{
  uint *p, i, un=0; ulong last=0;
  buffer_reserve(buf,2*n*sizeof(sort_data_long));
  p = buf->ptr;
  index_sort_long(id,n,1, p, buf->ptr);
  for(i=0;i<n;++i) {
    uint j = p[i]; ulong v = id[j];
    if(v==0) perm[j]=-1;
    else {
      if(v!=last) last=v, ++un;
      perm[j]=un-1;
    }
  }
  return un;
}

/* assign a master process for each id, in an arbitrary but consistent way */
static void setup_dofs(tuple_list *cdof, sint *perm, uint un, const ulong *id,
                       uint pid, crystal_data *crystal, buffer *buf)
{
  tuple_list shared; gs_data *gs; sint *out;
  uint i,ie,index, cn = unique_ids(un,id,perm,buf);
  tuple_list_init_max(cdof,cdof_mi,cdof_ml,cdof_mr,cn), cdof->n=cn;
  for(i=0;i<un;++i) if(perm[i]!=-1) cdof->vl[perm[i]]=id[i];
  gs = gs_data_setup(cn,(const ulong*)cdof->vl,1,crystal);
  gs_data_dump(&shared,gs);
  gs_data_free(gs);
  tuple_list_resize(&shared,shared.n+(cn+1));
  out = shared.vi+shared.mi*shared.n;
  for(i=0;i<cn;++i) *out++ = i, *out++ = pid, ++shared.n;
  *out++ = cn, *out++ = -1, ++shared.n;
  tuple_list_sort(&shared,1,buf);
  tuple_list_sort(&shared,0,buf);
  for(i=0;(index=shared.vi[2*i])<cn;i=ie) {
    ie=i+1;
    while(shared.vi[2*ie]==index) ++ie;
    cdof->vi[cdof_mi*index+cdof_proc] 
      = shared.vi[2*(i + cdof->vl[index]%(ie-i))+1];
  }
  tuple_list_free(&shared);
}

#ifdef AMG_DUMP
#  include "amg_dump.c"
#else
#  include "amg_run.c"
#endif

/*--------------------------------------------------------------------------
   FORTRAN Interface
  --------------------------------------------------------------------------*/

#include "fname.h"

#define xxtsetup   FORTRAN_NAME(crs_setup,CRS_SETUP)
#define xxtsolve   FORTRAN_NAME(crs_solve,CRS_SOLVE)
#define xxtstats   FORTRAN_NAME(crs_stats,CRS_STATS)
#define xxtfree    FORTRAN_NAME(crs_free ,CRS_FREE)

static amg_data **handle_array = 0;
static int handle_max = 0;
static int handle_n = 0;

#ifdef MPI
typedef MPI_Comm comm_ext_t;
#else
typedef int comm_ext_t;
#endif

void xxtsetup(sint *handle, const comm_ext_t *comm, const sint *np,
              const sint *n, const slong id[],
              const sint *nz, const sint Ai[], const sint Aj[], const real A[],
              const sint *null_space)
{
#ifndef MPI
  int crystal = 0;
#else
  crystal_data crystal;
  crystal_init(&crystal,*comm);
  if(crystal.num!=*np)
    fail("amg_setup: passed P=%d, but MPI_Comm_size gives P=%d\n",
         *np,crystal.num);
#endif
  if(handle_n==handle_max)
    handle_max+=handle_max/2+1,
    handle_array=trealloc(amg_data*,handle_array,handle_max);
  handle_array[handle_n]=amg_setup(*n,(const ulong*)id,
                                   *nz,(const uint*)Ai,(const uint*)Aj,A,
                                   *null_space,&crystal);
  *handle = handle_n++;
#ifdef MPI
  crystal_free(&crystal);
#endif
}

void xxtsolve(const sint *handle, real x[], const real b[])
{
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle])
    failwith("invalid handle to xxtsolve");
  amg_solve(x,handle_array[*handle],b);
}

void xxtstats(const sint *handle)
{
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle])
    failwith("invalid handle to xxtstats");
  amg_stats(handle_array[*handle]);
}

void xxtfree(sint *handle)
{
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle])
    failwith("invalid handle to xxtfree");
  amg_free(handle_array[*handle]);
  handle_array[*handle] = 0;
}

