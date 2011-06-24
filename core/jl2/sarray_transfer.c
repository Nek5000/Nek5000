#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#include "comm.h"
#include "mem.h"
#include "crystal.h"
#include "sort.h"

#define sarray_transfer_     PREFIXED_NAME(sarray_transfer_    )
#define sarray_transfer_ext_ PREFIXED_NAME(sarray_transfer_ext_)

/* Case 1:

  The input is an array with element type
    struct { ...; uint proc; ... } .
  Each element is transfered to where its .proc field indicates,
  and its .proc field is changed where it came from. (Hence, two calls
  to sarray_transfer are equivalent to some local permutation of elements.)

*/

static void pack(struct crystal *const cr,
                 const char *const input, const uint n, const size_t size,
                 const size_t off)
{
  const size_t row_size = ((size-sizeof(uint))+sizeof(uint)-1)/sizeof(uint);
  const size_t after = off + sizeof(uint), after_len = size-after;

  uint *perm, *out, *len=0;
  uint i, cn, p,lp = -(uint)1;
  
  /* get the permutation that orders by processor */
  perm = sortp(&cr->work,0, (const uint*)(input+off),n,size);
  
  /* pack into crystal router messages */
  buffer_reserve(&cr->data, n*(row_size+3)*sizeof(uint));
  out = cr->data.ptr; cn=0;
  for(i=0;i<n;++i) {
    const char *row = input + size*perm[i];
    memcpy(&p,row+off,sizeof(uint));
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

static void unpack(struct array *const A, const size_t size, const size_t off,
                   const int fixed, const int set_src, struct crystal *const cr)
{
  const size_t row_size = ((size-sizeof(uint))+sizeof(uint)-1)/sizeof(uint);
  const size_t after = off + sizeof(uint), after_len = size-after;

  const uint *in = cr->data.ptr, *in_end = in + cr->n;
  char *out;
  uint n=0;
  if(!fixed)
    while(in!=in_end) { uint len = in[2]; n+=len, in+=len+3; }
  else {
    const uint *original_in_end=in_end;
    const uint maxn = A->max * row_size;
    uint *in = cr->data.ptr;
    while(in!=original_in_end) {
      uint len = in[2]; n+=len;
      if(n<maxn) in+=len+3;
      else { in[2]-=(maxn-n); in_end=in+in[2]+3; in+=len+3; break; }
    }
    while(in!=original_in_end) { uint len = in[2]; n+=len, in+=len+3; }
  }
  n /= row_size;
  if(!fixed) array_reserve_(A,n,size,__FILE__,__LINE__);
  A->n = n;
  out = A->ptr;
  in = cr->data.ptr;
  #define UNPACK_MESSAGE(p,len) do { \
    for(;len;len-=row_size) { \
      memcpy(out,in,off); \
      memcpy(out+off,&p,sizeof(uint)); \
      memcpy(out+after,(const char*)in+off,after_len); \
      out += size; \
      in += row_size; \
    } \
  } while(0)
  if(set_src) {
    while(in!=in_end) {
      uint p=in[1], len=in[2]; in += 3;
      UNPACK_MESSAGE(p,len);
    }
  } else {
    const uint p = cr->comm.id;
    while(in!=in_end) {
      uint          len=in[2]; in += 3;
      UNPACK_MESSAGE(p,len);
    }
  }
  #undef UNPACK_MESSAGE
}

void sarray_transfer_(struct array *const A, const size_t size,
                      const size_t off, const int fixed, const int set_src,
                      struct crystal *const cr)
{
  pack(cr, A->ptr,A->n,size,off);
  crystal_router(cr);
  unpack(A,size,off, fixed,set_src, cr);
}

/* Case 2:

  The user supplies the destination of element i in proc[i].
  Here proc may be an external array of uints, or a field within the input.
  No source information is available on output.

*/

static void pack_ext(struct crystal *const cr,
                     const char *const input, const uint n, const size_t size,
                     const uint *const proc, const unsigned proc_stride)
{
  const size_t row_size = (size+sizeof(uint)-1)/sizeof(uint);

  uint *perm, *out, *len=0;
  uint i, cn, p,lp = -(uint)1;
  
  /* get the permutation that orders by processor */
  perm = sortp(&cr->work,0, proc,n,proc_stride);
  
  /* pack into crystal router messages */
  buffer_reserve(&cr->data, n*(row_size+3)*sizeof(uint));
  out = cr->data.ptr; cn=0;
  for(i=0;i<n;++i) {
    const char *row = input + size*perm[i];
    p = *(const uint *)((const char*)proc + proc_stride*perm[i]);
    if(p!=lp) {
      lp = p;
      *out++ = p;           /* target */
      *out++ = cr->comm.id; /* source */
      len = out++; *len=0;  /* length */
      cn += 3;
    }
    memcpy(out,row,size);
    out += row_size, cn += row_size, *len += row_size;
  }
  cr->n = cn;
}

static void unpack_ext(struct array *const A, const size_t size,
                       struct crystal *const cr)
{
  const size_t row_size = (size+sizeof(uint)-1)/sizeof(uint);

  const uint *in = cr->data.ptr, *in_end = in + cr->n;
  char *out;
  uint n=0;
  while(in!=in_end) { uint len = in[2]; n+=len, in+=len+3; }
  n /= row_size;
  array_reserve_(A,n,size,__FILE__,__LINE__);
  A->n = n;
  out = A->ptr;
  in = cr->data.ptr;
  while(in!=in_end) {
    uint len=in[2];
    in += 3;
    for(;len;len-=row_size) {
      memcpy(out,in,size);
      out += size;
      in += row_size;
    }
  }
}

void sarray_transfer_ext_(struct array *const A, const size_t size,
                          const uint *const proc, const unsigned proc_stride,
                          struct crystal *const cr)
{
  pack_ext(cr, A->ptr,A->n,size, proc,proc_stride);
  crystal_router(cr);
  unpack_ext(A,size,cr);
}

