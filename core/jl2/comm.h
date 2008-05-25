#ifndef COMM_H
#define COMM_H

#if !defined(NAME_H) || !defined(ERRMEM_H)  || !defined(TYPES_H)
#warning "comm.h" requires "name.h", "errmem.h" and "types.h"
#endif

#ifdef MPI
#include <mpi.h>
typedef MPI_Comm comm_ext_t;
typedef MPI_Request comm_req_t;
#else
typedef int comm_ext_t;
typedef int comm_req_t;
#endif

#ifdef PREFIX
#  define comm_partial_sum_ul TOKEN_PASTE(PREFIX,comm_partial_sum_ul)
#  define comm_allreduce      TOKEN_PASTE(PREFIX,comm_allreduce     )
#endif

typedef struct {
  uint id, np;
  comm_ext_t comm;
} comm_t;

static void comm_init(comm_t *c, comm_ext_t comm)
{
#ifdef MPI
  int i;
  MPI_Comm_dup(comm, &c->comm);
  MPI_Comm_rank(c->comm,&i), c->id=i;
  MPI_Comm_size(c->comm,&i), c->np=i;
#else
  c->id = 0, c->np = 1;
#endif
}

static void comm_init_check(comm_t *c, comm_ext_t comm, uint np)
{
  comm_init(c,comm);
  if(c->np != np) {
#ifdef MPI
    fail("comm_init_check: passed P=%u, but MPI_Comm_size gives P=%u\n",
         np,c->np);
#else
    fail("comm_init_check: passed P=%u, but not compiled with -DMPI\n",np);
#endif
  }
}

static void comm_dup_(comm_t *d, const comm_t *s, const char *file)
{
  d->id = s->id, d->np = s->np;
#ifdef MPI
  MPI_Comm_dup(s->comm,&d->comm);
#else
  if(s->np!=1) fail("%s not compiled with -DMPI\n",file);
#endif
}
#define comm_dup(d,s) comm_dup_(d,s,__FILE__)

static void comm_free(comm_t *c)
{
#ifdef MPI
  MPI_Comm_free(&c->comm);
#endif
}

static double comm_time(void)
{
#ifdef MPI
  return MPI_Wtime();
#else
  return 0;
#endif
}

static void comm_barrier(const comm_t *c)
{
#ifdef MPI
  MPI_Barrier(c->comm);
#endif
}

static void comm_recv(const comm_t *c, void *p, size_t n, uint src, int tag)
{
#ifdef MPI
  MPI_Recv(p,n,MPI_UNSIGNED_CHAR,src,tag,c->comm,MPI_STATUS_IGNORE);
#endif
}

static void comm_send(const comm_t *c, void *p, size_t n, uint dst, int tag)
{
#ifdef MPI
  MPI_Send(p,n,MPI_UNSIGNED_CHAR,dst,tag,c->comm);
#endif
}

static void comm_irecv(comm_req_t *req, const comm_t *c,
                       void *p, size_t n, uint src, int tag)
{
#ifdef MPI
  MPI_Irecv(p,n,MPI_UNSIGNED_CHAR,src,tag,c->comm,req);
#endif
}

static void comm_isend(comm_req_t *req, const comm_t *c,
                       void *p, size_t n, uint dst, int tag)
{
#ifdef MPI
  MPI_Isend(p,n,MPI_UNSIGNED_CHAR,dst,tag,c->comm,req);
#endif
}

static void comm_wait(comm_req_t *req, unsigned n)
{
#ifdef MPI
  MPI_Waitall(n,req,MPI_STATUSES_IGNORE);
#endif
}

void comm_partial_sum_ul(ulong sum[2], const comm_t *comm, ulong v);

#ifdef GS_DEFS_H
void comm_allreduce(void *v, void *buf, uint vn, const comm_t *comm,
                    gs_dom_t dom, gs_op_t op);
#endif

#endif
