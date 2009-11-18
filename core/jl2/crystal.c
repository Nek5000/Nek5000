/*------------------------------------------------------------------------------
  
  Crystal Router
  
  Accomplishes all-to-all communication in log P msgs per proc
  The routine is low-level; the format of the input/output is an
  array of integers, consisting of a sequence of messages with format:
  
      target proc
      source proc
      m
      integer
      integer
      ...
      integer  (m integers in total)

  Before crystal_router is called, the source of each message should be
  set to this proc id; upon return from crystal_router, the target of each
  message will be this proc id.
  
  Example Usage:
  
    crystal_data crystal;
    
    crystal_init(&crystal, &comm);  // makes an internal copy of comm
    
    crystal.n = ... ;  // total number of integers
    buffer_reserve(&crystal.data, crystal.n * sizeof(uint));
    ... // fill crystal.data.ptr with messages
    crystal_router(&crystal);
    
    crystal_free(&crystal);
    
  ----------------------------------------------------------------------------*/

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "name.h"
#include "fail.h"
#include "types.h"
#include "comm.h"
#include "mem.h"

#define crystal_init   PREFIXED_NAME(crystal_init  )
#define crystal_free   PREFIXED_NAME(crystal_free  )
#define crystal_router PREFIXED_NAME(crystal_router)

typedef struct {
  struct comm comm;
  buffer data, work;
  uint n;
} crystal_data;

void crystal_init(crystal_data *p, const struct comm *comm)
{
  comm_dup(&p->comm, comm);
  buffer_init(&p->data,1000);
  buffer_init(&p->work,1000);
  p->n=0;
}

void crystal_free(crystal_data *p)
{
  comm_free(&p->comm);
  buffer_free(&p->data);
  buffer_free(&p->work);
}

static void uintcpy(uint *dst, const uint *src, uint n)
{
  if(dst+n<=src)    memcpy (dst,src,n*sizeof(uint));
  else if(dst!=src) memmove(dst,src,n*sizeof(uint));
}

static uint crystal_move(crystal_data *p, uint cutoff, int send_hi)
{
  uint len, *src, *end;
  uint *keep = p->data.ptr, *send;
  buffer_reserve(&p->work,p->n*sizeof(uint)), send = p->work.ptr;
  if(send_hi) { /* send hi, keep lo */
    for(src=keep,end=keep+p->n; src<end; src+=len) {
      len = 3 + src[2];
      if(src[0]>=cutoff) memcpy (send,src,len*sizeof(uint)), send+=len;
      else               uintcpy(keep,src,len),              keep+=len;
    }
  } else      { /* send lo, keep hi */
    for(src=keep,end=keep+p->n; src<end; src+=len) {
      len = 3 + src[2];
      if(src[0]< cutoff) memcpy (send,src,len*sizeof(uint)), send+=len;
      else               uintcpy(keep,src,len),              keep+=len;
    }
  }
  p->n = keep - (uint*)p->data.ptr;
  return send - (uint*)p->work.ptr;
}

static void crystal_exchange(crystal_data *p, uint send_n, uint targ,
                             int recvn, int tag)
{
  comm_req req[3];
  uint count[2] = {0,0}, sum, *recv[2];

  if(recvn)   
    comm_irecv(&req[1],&p->comm, &count[0],sizeof(uint), targ        ,tag);
  if(recvn==2)
    comm_irecv(&req[2],&p->comm, &count[1],sizeof(uint), p->comm.id-1,tag);
  comm_isend(&req[0],&p->comm, &send_n,sizeof(uint), targ,tag);
  comm_wait(req,recvn+1);
  
  sum = p->n + count[0] + count[1];
  buffer_reserve(&p->data,sum*sizeof(uint));
  recv[0] = (uint*)p->data.ptr + p->n, recv[1] = recv[0] + count[0];
  p->n = sum;
  
  if(recvn)    comm_irecv(&req[1],&p->comm,
                          recv[0],count[0]*sizeof(uint), targ        ,tag+1);
  if(recvn==2) comm_irecv(&req[2],&p->comm,
                          recv[1],count[1]*sizeof(uint), p->comm.id-1,tag+1);
  comm_isend(&req[0],&p->comm, p->work.ptr,send_n*sizeof(uint), targ,tag+1);
  comm_wait(req,recvn+1);
}

void crystal_router(crystal_data *p)
{
  uint bl=0, bh, nl;
  uint id = p->comm.id, n=p->comm.np;
  uint send_n, targ, tag = 0;
  int send_hi, recvn;
  while(n>1) {
    nl = (n+1)/2, bh = bl+nl;
    send_hi = id<bh;
    send_n = crystal_move(p,bh,send_hi);
    recvn = 1, targ = n-1-(id-bl)+bl;
    if(id==targ) targ=bh, recvn=0;
    if(n&1 && id==bh) recvn=2;
    crystal_exchange(p,send_n,targ,recvn,tag);
    if(id<bh) n=nl; else n-=nl,bl=bh;
    tag += 2;
  }
}
