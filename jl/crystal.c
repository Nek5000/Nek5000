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

  Usage:
  
      MPI_Comm comm = ... ;
      crystal_data crystal;
      
      crystal_init(&crystal, comm);  // initialize the data structure
      // now crystal.id  = this proc
      // and crystal.num = num of procs
      
      // allocate space for at least MAX ints
      buffer_reserve(&crystal->all->buf, MAX*sizeof(uint));
      
      // fill up ((uint*)crystal->all->buf.ptr)[0 ... n-1]
      // and set crystal->all->n
      
      crystal_router(&crystal);
      
      // incoming messages available as
      // ((uint*)crystal->all->buf.ptr)[0 ... crystal->all->n-1]
      
      crystal_free(&crystal); // release acquired memory

  ----------------------------------------------------------------------------*/

#ifdef MPI

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <mpi.h>

#include "errmem.h"
#include "types.h"

typedef struct { uint n; buffer buf; } crystal_buf;

typedef struct {
  crystal_buf buffers[3];
  crystal_buf *all, *keep, *send;
  MPI_Comm comm;
  uint num, id;
} crystal_data;

#define crystal_free crystal_old_free

void crystal_init(crystal_data *p, MPI_Comm comm)
{
  int num,id;
  buffer_init(&p->buffers[0].buf,1024);
  buffer_init(&p->buffers[1].buf,1024);
  buffer_init(&p->buffers[2].buf,1024);
  p->all=&p->buffers[0];
  p->keep=&p->buffers[1];
  p->send=&p->buffers[2];
  memcpy(&p->comm,&comm,sizeof(MPI_Comm));
  MPI_Comm_rank(comm,&id ); p->id =id ;
  MPI_Comm_size(comm,&num); p->num=num;
}

void crystal_free(crystal_data *p)
{
  buffer_free(&p->buffers[0].buf);
  buffer_free(&p->buffers[1].buf);
  buffer_free(&p->buffers[2].buf);
}

static void crystal_partition(crystal_data *p, uint cutoff,
                              crystal_buf *lo, crystal_buf *hi)
{
  const uint *src = p->all->buf.ptr;
  const uint *end = src+p->all->n;
  uint *target, *lop, *hip;
  lo->n=hi->n=0;
  buffer_reserve(&lo->buf,p->all->n*sizeof(uint));
  buffer_reserve(&hi->buf,p->all->n*sizeof(uint));
  lop = lo->buf.ptr, hip = hi->buf.ptr;
  while(src!=end) {
    uint chunk_len = 3 + src[2];
    if(src[0]<cutoff) target=lop,lo->n+=chunk_len,lop+=chunk_len;
                 else target=hip,hi->n+=chunk_len,hip+=chunk_len;
    memcpy(target,src,chunk_len*sizeof(uint));
    src+=chunk_len;
  }
}

static void crystal_send(crystal_data *p, uint target, int recvn)
{
  MPI_Request req[3];
  MPI_Status status[3];
  uint count[2]={0,0},sum,*recv[2];
  crystal_buf *t;
  int i;
  
  for(i=0;i<recvn;++i)
    MPI_Irecv(&count[i],1,UINT_MPI,target+i,target+i,p->comm,&req[i+1]);
  MPI_Isend(&p->send->n,1,UINT_MPI,target  ,p->id   ,p->comm,&req[  0]);
  MPI_Waitall(recvn+1,req,status);
  sum = p->keep->n;
  for(i=0;i<recvn;++i) sum+=count[i];
  buffer_reserve(&p->keep->buf,sum*sizeof(uint));
  recv[0]=p->keep->buf.ptr;
  recv[0]+=p->keep->n;
  recv[1]=recv[0]+count[0];
  p->keep->n=sum;

  MPI_Isend(p->send->buf.ptr,p->send->n*sizeof(uint),
            MPI_UNSIGNED_CHAR,target,p->id,p->comm,&req[0]);
  if(recvn) {
    MPI_Irecv(recv[0],count[0]*sizeof(uint),MPI_UNSIGNED_CHAR,
              target,target,p->comm,&req[1]);
    if(recvn==2)
    MPI_Irecv(recv[1],count[1]*sizeof(uint),MPI_UNSIGNED_CHAR,
              target+1,target+1,p->comm,&req[2]);
  }
  MPI_Waitall(recvn+1,req,status);

  t=p->all,p->all=p->keep,p->keep=t;
}

void crystal_router(crystal_data *p)
{
  uint bl=0, bh, n=p->num, nl, target;
  int recvn;
  crystal_buf *lo, *hi;
  while(n>1) {
    nl = n/2, bh = bl+nl;
    if(p->id<bh)
      target=p->id+nl,recvn=(n&1 && p->id==bh-1)?2:1   ,lo=p->keep,hi=p->send;
    else
      target=p->id-nl,recvn=(target==bh)?(--target,0):1,hi=p->keep,lo=p->send;
    crystal_partition(p,bh,lo,hi);
    crystal_send(p,target,recvn);
    if(p->id<bh) n=nl; else n-=nl,bl=bh;
  }
}

#endif

