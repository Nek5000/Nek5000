#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "name.h"
#include "errmem.h"
#include "types.h"
#include "sort.h"

#ifdef PREFIX
#  define sarray_permute_ TOKEN_PASTE(PREFIX,sarray_permute_)
#  define sarray_sort_begin_ TOKEN_PASTE(PREFIX,sarray_sort_begin_)
#  define sarray_sort_cont_ TOKEN_PASTE(PREFIX,sarray_sort_cont_)
#  define sarray_sort_end_ TOKEN_PASTE(PREFIX,sarray_sort_end_)
#endif

static size_t sarray_permute_worksize_(size_t n, size_t work_off,
                                       size_t align, size_t size)
{
  return align_as_(align,work_off+n*size);
}

void sarray_permute_(void *A, size_t n, const uint *perm,
                     buffer *buf, int resize, size_t align, size_t size)
{
  char *work, *src, *dst;
  const uint *pe;
  size_t work_size = sarray_permute_worksize_(n,buf->n,align,size);
  if(resize) buffer_reserve(buf,work_size);
  else if(buf->max<work_size)
    failwith(__FILE__ ": sarray_permute_() not given enough workspace");
  work = (char*)buf->ptr + align_as_(align,buf->n);
  src = A, dst = work;
  for(pe=perm+n;perm!=pe;++perm,dst+=size) memcpy(dst,src+(*perm)*size,size);
  memcpy(A,work,n*size);
}

static size_t sarray_sort_worksize_(size_t n, size_t work_off,
                                    size_t align, size_t size)
{
  size_t perm_size = align_as(uint,work_off+n*sizeof(uint));
  size_t n1 = sortp_worksize(n,perm_size),
         n2 = sortp_long_worksize(n,perm_size),
         n3 = sarray_permute_worksize_(n,perm_size,align,size);
  return n1>n2 ? (n1>n3?n1:n3) : (n2>n3?n2:n3);
}

void sarray_sort_begin_(void *A, size_t n, int is_long, buffer *buf,
                        size_t align, size_t off, size_t size)
{
  size_t bufn = buf->n;
  uint *perm;
  buffer_reserve(buf,sarray_sort_worksize_(n,bufn,align,size));
  perm = align_ptr(uint,buf->ptr,bufn);
  buf->n = (char*)&perm[n] - (char*)buf->ptr;
  if(is_long)
    sortp_long(perm,0, (ulong*)((char*)A + off),n,size, buf,0);
  else
    sortp     (perm,0, (uint *)((char*)A + off),n,size, buf,0);
  buf->n = bufn;
}

void sarray_sort_cont_(void *A, size_t n, int is_long, buffer *buf,
                       size_t off, size_t size)
{
  size_t bufn = buf->n;
  uint *perm = align_ptr(uint,buf->ptr,bufn);
  buf->n = (char*)&perm[n] - (char*)buf->ptr;
  if(is_long)
    sortp_long(perm,1, (ulong*)((char*)A + off),n,size, buf,0);
  else
    sortp     (perm,1, (uint *)((char*)A + off),n,size, buf,0);
  buf->n = bufn;
}

void sarray_sort_end_(void *A, size_t n, buffer *buf, size_t align, size_t size)
{
  size_t bufn = buf->n;
  const uint *perm = align_ptr(uint,buf->ptr,bufn);
  buf->n = (char*)&perm[n] - (char*)buf->ptr;
  sarray_permute_(A,n,perm, buf,0, align,size);
  buf->n = bufn;
}
