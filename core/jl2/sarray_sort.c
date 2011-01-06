#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "c99.h"
#include "name.h"
#include "types.h"
#include "fail.h"
#include "mem.h"
#include "sort.h"

#define sarray_permute_ PREFIXED_NAME(sarray_permute_)

void sarray_permute_(size_t align, size_t size, void *A, size_t n, buffer *buf)
{
  const uint *perm, *pe;
  char *work, *src, *dst;
  buffer_reserve(buf,align_as_(align,n*sizeof(uint)+n*size));
  perm = buf->ptr;
  work = (char*)buf->ptr + align_as_(align,n*sizeof(uint));
  src = A, dst = work;
  for(pe=perm+n;perm!=pe;++perm,dst+=size) memcpy(dst,src+(*perm)*size,size);
  memcpy(A,work,n*size);
}
