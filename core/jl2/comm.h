#ifndef COMM_H
#define COMM_H

/* requires:
     <stddef.h>            for size_t
     <stdlib.h>            for exit
     "fail.h", "types.h"
     "gs_defs.h"           for comm_allreduce, comm_scan, comm_reduce_T
*/

#if !defined(FAIL_H) || !defined(TYPES_H)
#warning "comm.h" requires "fail.h" and "types.h"
#endif

#ifdef MPI
#include <mpi.h>
typedef MPI_Comm comm_ext;
typedef MPI_Request comm_req;
#else
typedef int comm_ext;
typedef int comm_req;
typedef int MPI_Fint;
#endif

#define comm_allreduce PREFIXED_NAME(comm_allreduce)
#define comm_scan      PREFIXED_NAME(comm_scan     )
#define comm_dot       PREFIXED_NAME(comm_dot      )

struct comm {
  uint id, np;
  comm_ext c;
};

static void comm_exit(const struct comm *c, int status);
static void comm_init(struct comm *c, comm_ext ce);
static void comm_init_check(struct comm *c, MPI_Fint ce, uint np);
/* (macro) static void comm_dup(struct comm *d, const struct comm *s); */
static void comm_free(struct comm *c);
static double comm_time(void);
static void comm_barrier(const struct comm *c);
static void comm_recv(const struct comm *c, void *p, size_t n,
                      uint src, int tag);
static void comm_send(const struct comm *c, void *p, size_t n,
                      uint dst, int tag);
static void comm_irecv(comm_req *req, const struct comm *c,
                       void *p, size_t n, uint src, int tag);
static void comm_isend(comm_req *req, const struct comm *c,
                       void *p, size_t n, uint dst, int tag);
static void comm_wait(comm_req *req, int n);

double comm_dot(const struct comm *comm, double *v, double *w, uint n);

#ifdef GS_DEFS_H
void comm_allreduce(const struct comm *com, gs_dom dom, gs_op op,
                          void *v, uint vn, void *buf);
void comm_scan(void *scan, const struct comm *com, gs_dom dom, gs_op op,
               const void *v, uint vn, void *buffer);

#define DEFINE_REDUCE(T) \
T PREFIXED_NAME(comm_reduce__##T)( \
    const struct comm *comm, gs_op op, const T *in, uint n); \
static T comm_reduce_##T(const struct comm *c, gs_op op, const T *v, uint vn) \
{ return PREFIXED_NAME(comm_reduce__##T)(c,op,v,vn); }
GS_FOR_EACH_DOMAIN(DEFINE_REDUCE)
#undef DEFINE_REDUCE

#define comm_reduce_sint \
    TYPE_LOCAL(comm_reduce_int,comm_reduce_long,comm_reduce_long_long)
#define comm_reduce_slong \
   TYPE_GLOBAL(comm_reduce_int,comm_reduce_long,comm_reduce_long_long)

#endif

/*----------------------------------------------------------------------------
  Code for static (inline) functions
  ----------------------------------------------------------------------------*/

static void comm_init(struct comm *c, comm_ext ce)
{
#ifdef MPI
  int i;
  MPI_Comm_dup(ce, &c->c);
  MPI_Comm_rank(c->c,&i), c->id=i;
  MPI_Comm_size(c->c,&i), c->np=i;
#else
  c->id = 0, c->np = 1;
#endif
}

static void comm_init_check(struct comm *c, MPI_Fint ce, uint np)
{
#ifdef MPI
  comm_init(c,MPI_Comm_f2c(ce));
  if(c->np != np)
    fail(1,"comm_init_check: passed P=%u, but MPI_Comm_size gives P=%u",
         np,c->np);
#else
  comm_init(c,0);
  if(np != 1)
    fail(1,"comm_init_check: passed P=%u, but not compiled with -DMPI",np);
#endif
}

static void comm_dup_(struct comm *d, const struct comm *s, const char *file)
{
  d->id = s->id, d->np = s->np;
#ifdef MPI
  MPI_Comm_dup(s->c,&d->c);
#else
  if(s->np!=1) fail(1,"%s not compiled with -DMPI\n",file);
#endif
}
#define comm_dup(d,s) comm_dup_(d,s,__FILE__)

static void comm_free(struct comm *c)
{
#ifdef MPI
  MPI_Comm_free(&c->c);
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

static void comm_barrier(const struct comm *c)
{
#ifdef MPI
  MPI_Barrier(c->c);
#endif
}

static void comm_exit(const struct comm *c, int status)
{
  comm_barrier(c);
#ifdef MPI
  MPI_Finalize();
#endif
  exit(status);
}


static void comm_recv(const struct comm *c, void *p, size_t n,
                      uint src, int tag)
{
#ifdef MPI
# ifndef MPI_STATUS_IGNORE
  MPI_Status stat;
  MPI_Recv(p,n,MPI_UNSIGNED_CHAR,src,tag,c->c,&stat);
# else  
  MPI_Recv(p,n,MPI_UNSIGNED_CHAR,src,tag,c->c,MPI_STATUS_IGNORE);
# endif
#endif
}

static void comm_send(const struct comm *c, void *p, size_t n,
                      uint dst, int tag)
{
#ifdef MPI
  MPI_Send(p,n,MPI_UNSIGNED_CHAR,dst,tag,c->c);
#endif
}

static void comm_irecv(comm_req *req, const struct comm *c,
                       void *p, size_t n, uint src, int tag)
{
#ifdef MPI
  MPI_Irecv(p,n,MPI_UNSIGNED_CHAR,src,tag,c->c,req);
#endif
}

static void comm_isend(comm_req *req, const struct comm *c,
                       void *p, size_t n, uint dst, int tag)
{
#ifdef MPI
  MPI_Isend(p,n,MPI_UNSIGNED_CHAR,dst,tag,c->c,req);
#endif
}

static void comm_wait(comm_req *req, int n)
{
#ifdef MPI
# ifndef MPI_STATUSES_IGNORE
  MPI_Status status[8];
  while(n>=8) MPI_Waitall(8,req,status), req+=8, n-=8;
  if(n>0) MPI_Waitall(n,req,status);
# else
  MPI_Waitall(n,req,MPI_STATUSES_IGNORE);
# endif  
#endif
}

#endif
