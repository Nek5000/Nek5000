#ifndef SARRAY_SORT_H
#define SARRAY_SORT_H

#if !defined(SORT_H)
#warning "sarray_sort.h" requires "sort.h"
#endif

/*------------------------------------------------------------------------------
  
  Array of Structs Sort
  
  buffer *buf;
  typedef struct { ... } T;
  T A[n];

  sarray_sort(T,A,n, field_name,is_long, buf)
    - sort A according to the struct field "field_name",
      which is a ulong/uint field according as is_long is true/false

  sarray_sort_two(T,A,n, field1,is_long1, field2,is_long2, buf)
    - sort A by field1 then field2

  sarray_permute(T,A,n, buf);
    - permute A according to the permutation in buf
      A[0] <- A[perm[0]], etc.
      where uint *perm = buf->ptr   (see "sort.h")

  ----------------------------------------------------------------------------*/


#define sarray_permute_ PREFIXED_NAME(sarray_permute_)

void sarray_permute_(size_t align, size_t size, void *A, size_t n, buffer *buf);

#define sarray_permute(T,A,n, buf) \
  sarray_permute_(ALIGNOF(T),sizeof(T),A,n,buf)

#define sarray_sort(T,A,n, field,is_long, buf) do { \
  if(is_long) \
    sortp_long(buf,0, (ulong*)((char*)(A)+offsetof(T,field)),n,sizeof(T)); \
  else \
    sortp     (buf,0, (uint *)((char*)(A)+offsetof(T,field)),n,sizeof(T)); \
  sarray_permute(T,A,n, buf); \
} while (0)

#define sarray_sort_two(T,A,n, field1,is_long1, field2,is_long2, buf) do { \
  if(is_long2) \
    sortp_long(buf,0, (ulong*)((char*)(A)+offsetof(T,field2)),n,sizeof(T)); \
  else \
    sortp     (buf,0, (uint *)((char*)(A)+offsetof(T,field2)),n,sizeof(T)); \
  if(is_long1) \
    sortp_long(buf,1, (ulong*)((char*)(A)+offsetof(T,field1)),n,sizeof(T)); \
  else \
    sortp     (buf,1, (uint *)((char*)(A)+offsetof(T,field1)),n,sizeof(T)); \
  sarray_permute(T,A,n, buf); \
} while (0)

#endif
