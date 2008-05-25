#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "name.h"
#include "errmem.h"
#include "types.h"
#include "comm.h"
#include "crystal.h"
#include "sort.h"

#ifdef PREFIX
#  define sarray_transfer_ TOKEN_PASTE(PREFIX,sarray_transfer_)
#endif

static void pack(const char *input, uint n, crystal_data *cr,
                 size_t off, size_t size)
{
  const size_t row_size = ((size-sizeof(uint))+sizeof(uint)-1)/sizeof(uint);
  const size_t after = off + sizeof(uint), after_len = size-after;

  uint *perm, *out, *len=0;
  uint i, cn, p,lp = -(uint)1;
  
  /* get the permutation that orders by processor */
  buffer_reserve(&cr->work, n*sizeof(uint));
  perm = cr->work.ptr;
  cr->data.n = 0;
  sortp(perm,0, (const uint*)(input + off),n,size, &cr->data,1);
  
  /* pack into crystal router messages */
  buffer_reserve(&cr->data, n*(row_size+3)*sizeof(uint));
  out = cr->data.ptr; cn=0;
  for(i=0;i<n;++i) {
    const char *row = input + size*perm[i];
    p = *(uint*)(row+off);
    if(p!=lp) {
      lp = p;
      *out++ = p;           /* target */
      *out++ = cr->comm.id; /* source */
      len = out++; *len=0;  /* length */
      cn += 3;
    }
    memcpy(out,row,off);
    memcpy((char*)out + off,row+after,after_len);
    out += row_size, cn += row_size, *len += row_size;
  }
  cr->n = cn;
}

static void unpack(array *A, crystal_data *cr, size_t off, size_t size)
{
  const size_t row_size = ((size-sizeof(uint))+sizeof(uint)-1)/sizeof(uint);
  const size_t after = off + sizeof(uint), after_len = size-after;

  const uint *in = cr->data.ptr, *in_end = in + cr->n;
  char *out;
  uint n=0;
  while(in!=in_end) { uint len = in[2]; n+=len, in+=len+3; }
  n /= row_size;
  array_reserve_(A,n,size,__FILE__);
  A->n = n;
  out = A->ptr;
  in = cr->data.ptr;
  while(in!=in_end) {
    uint p=in[1], len=in[2];
    in += 3;
    for(;len;len-=row_size) {
      memcpy(out,in,off);
      memcpy(out+off,&p,sizeof(uint));
      memcpy(out+after,(char*)in+off,after_len);
      out += size;
      in += row_size;
    }
  }
}

void sarray_transfer_(array *A, crystal_data *cr, size_t off, size_t size)
{
  pack(A->ptr,A->n,cr,off,size);
  crystal_router(cr);
  unpack(A,cr,off,size);
}

