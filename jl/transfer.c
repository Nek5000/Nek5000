#ifdef MPI

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include "errmem.h"
#include "types.h"
#include "minmax.h"
#include "sort.h"
#include "tuple_list.h"
#include "crystal.h"

#define UINT_PER_X(X) ((sizeof(X)+sizeof(uint)-1)/sizeof(uint))
#define UINT_PER_REAL UINT_PER_X(real)
#define UINT_PER_LONG UINT_PER_X(slong)

/*------------------------------------------------------------------------------
  
  Transfer
 
  Treats one integer (not long) member of the tuple list as a target proc;
  Sends out tuples accordingly, using the crystal router.
  Target proc member overwritten with source proc.
  
  dynamic: non-zero if the tuple list should grow to accomodate arrivals
  tl:      the tuple list
  pf:      which tuple member specifies target proc
  crystal: an initialized crystal router structure (cf. crystal.h)

  ----------------------------------------------------------------------------*/

void transfer(int dynamic, tuple_list *tl,
              unsigned pf, crystal_data *crystal)
{
  const unsigned mi=tl->mi,ml=tl->ml,mr=tl->mr;
  const unsigned tsize = (mi-1) + ml*UINT_PER_LONG + mr*UINT_PER_REAL;
  sint p, lp = -1;
  sint *ri; slong *rl; real *rr;
  uint i, j, *buf, *len=0, *buf_end;

  /* sort to group by target proc */
  tuple_list_sort(tl,pf,&crystal->all->buf);

  /* pack into buffer for crystal router */
  buffer_reserve(&crystal->all->buf,(tl->n*(3+tsize))*sizeof(uint));
  crystal->all->n=0, buf = crystal->all->buf.ptr;
  ri=tl->vi,rl=tl->vl,rr=tl->vr;
  for(i=tl->n;i;--i) {
    p = ri[pf];
    if(p!=lp) {
      lp = p;
      *buf++ = p;           /* target */
      *buf++ = crystal->id; /* source */
      len = buf++; *len=0;  /* length */
      crystal->all->n += 3;
    }
    for(j=0;j<mi;++j,++ri) if(j!=pf) *buf++ = *ri;
    for(j=ml;j;--j,++rl)
      memcpy(buf,rl,sizeof(slong)), buf+=UINT_PER_LONG;
    for(j=mr;j;--j,++rr)
      memcpy(buf,rr,sizeof(real )), buf+=UINT_PER_REAL;
    *len += tsize, crystal->all->n += tsize;
  }
  
  crystal_router(crystal);
  
  /* unpack */
  buf = crystal->all->buf.ptr, buf_end = buf + crystal->all->n;
  tl->n = 0;
  ri=tl->vi,rl=tl->vl,rr=tl->vr;
  while(buf != buf_end) {
    sint p, len;
    buf++;        /* target ( == this proc ) */
    p = *buf++;   /* source */
    len = *buf++; /* length */
    while(len>0) {
      if(tl->n==tl->max) {
        if(!dynamic) { tl->n = tl->max + 1; return; }
        tuple_list_grow(tl);
        ri = tl->vi + mi*tl->n, rl = tl->vl + ml*tl->n, rr = tl->vr + mr*tl->n;
      }
      ++tl->n;
      for(j=0;j<mi;++j) if(j!=pf) *ri++ = *buf++; else *ri++ = p;
      for(j=ml;j;--j) memcpy(rl++,buf,sizeof(slong)), buf+=UINT_PER_LONG;
      for(j=mr;j;--j) memcpy(rr++,buf,sizeof(real )), buf+=UINT_PER_REAL;
      len-=tsize;
    }
  }
}

#endif

