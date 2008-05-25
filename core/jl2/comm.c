#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "name.h"
#include "types.h"
#include "errmem.h"
#include "gs_defs.h"
#include "gs_local.h"
#include "comm.h"

/* sets sum[0] = sum_{p<id} v, sum[1] = sum_p v */
void comm_partial_sum_ul(ulong sum[2], const comm_t *comm, ulong v)
{
  comm_req_t req[2];
  const uint id=comm->id, np=comm->np;
  uint n = np, c=1, odd=0, base=0;
  ulong r[2];
  while(n>1) {
    odd=(odd<<1)|(n&1);
    c<<=1, n>>=1;
    if(id>=base+n) c|=1, base+=n, n+=(odd&1);
  }
  sum[0]=0, sum[1]=v;
  while(n<np) {
    if(c&1) n-=(odd&1), base-=n;
    c>>=1, n<<=1, n+=(odd&1);
    odd>>=1;
    if(base==id) {
      comm_irecv(&req[0],comm, &r[0]  ,sizeof(ulong), id+n/2,id+n/2);
      comm_isend(&req[1],comm, &sum[1],sizeof(ulong), id+n/2,id);
      comm_wait(req,2);
      sum[1]+=r[0];
    } else {
      comm_irecv(&req[0],comm, &sum[0],sizeof(ulong), base,base);
      comm_isend(&req[1],comm, &sum[1],sizeof(ulong), base,id);
      comm_wait(req,2);
      break;
    }
  }
  while(n>1) {
    if(base==id) {
      comm_send(comm, sum,2*sizeof(ulong), id+n/2,id);
    } else {
      comm_recv(comm, r,2*sizeof(ulong), base,base);
      sum[0]+=r[0], sum[1]=r[1];
    }
    odd=(odd<<1)|(n&1);
    c<<=1, n>>=1;
    if(id>=base+n) c|=1, base+=n, n+=(odd&1);
  }
}

GS_DEFINE_DOM_SIZES()

static void allreduce_imp(void *v, void *buf, uint vn, const comm_t *comm,
                          gs_dom_t dom, gs_op_t op)
{
  size_t total_size = vn*gs_dom_size[dom];
  const uint id=comm->id, np=comm->np;
  uint n = np, c=1, odd=0, base=0;
  while(n>1) {
    odd=(odd<<1)|(n&1);
    c<<=1, n>>=1;
    if(id>=base+n) c|=1, base+=n, n+=(odd&1);
  }
  while(n<np) {
    if(c&1) n-=(odd&1), base-=n;
    c>>=1, n<<=1, n+=(odd&1);
    odd>>=1;
    if(base==id) {
      comm_recv(comm, buf,total_size, id+n/2,id+n/2);
      gs_gather_array(v,buf,vn, dom,op);
    } else {
      comm_send(comm, v,total_size, base,id);
      break;
    }
  }
  while(n>1) {
    if(base==id)
      comm_send(comm, v,total_size, id+n/2,id);
    else
      comm_recv(comm, v,total_size, base,base);
    odd=(odd<<1)|(n&1);
    c<<=1, n>>=1;
    if(id>=base+n) c|=1, base+=n, n+=(odd&1);
  }
}

void comm_allreduce(void *v, void *buf, uint vn, const comm_t *comm,
                    gs_dom_t dom, gs_op_t op)
{
  if(vn==0) return;
#ifdef MPI
  {
    MPI_Datatype mpitype;
    MPI_Op mpiop;
    switch(dom) { case gs_double: mpitype=MPI_DOUBLE; break;
                  case gs_float:  mpitype=MPI_FLOAT;  break;
                  case gs_int:    mpitype=MPI_INT;    break;
                  case gs_long:   mpitype=MPI_LONG;   break;
                  default:        goto comm_allreduce_byhand;
    }
    switch(op) { case gs_add: mpiop=MPI_SUM;  break;
                 case gs_mul: mpiop=MPI_PROD; break;
                 case gs_min: mpiop=MPI_MIN;  break;
                 case gs_max: mpiop=MPI_MAX;  break;
                 default:        goto comm_allreduce_byhand;
    }
    MPI_Allreduce(v,buf,vn,mpitype,mpiop,comm->comm);
    memcpy(v,buf,vn*gs_dom_size[dom]);
    return;
  }
#endif
#ifdef MPI
comm_allreduce_byhand:
  allreduce_imp(v,buf,vn,comm,dom,op);
#endif
}
