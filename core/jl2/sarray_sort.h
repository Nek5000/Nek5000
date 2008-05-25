#ifndef SARRAY_SORT_H
#define SARRAY_SORT_H

#if !defined(ERRMEM_H) || !defined(TYPES_H)
#warning "sarray_sort.h" requires "errmem.h" and "types.h"
#endif

#ifdef PREFIX
#  define sarray_permute_ TOKEN_PASTE(PREFIX,sarray_permute_)
#  define sarray_sort_begin_ TOKEN_PASTE(PREFIX,sarray_sort_begin_)
#  define sarray_sort_cont_ TOKEN_PASTE(PREFIX,sarray_sort_cont_)
#  define sarray_sort_end_ TOKEN_PASTE(PREFIX,sarray_sort_end_)
#endif

void sarray_permute_(void *A, size_t n, const uint *perm,
                     buffer *buf, int resize, size_t align, size_t size);
void sarray_sort_begin_(void *A, size_t n, int is_long, buffer *buf,
                        size_t align, size_t off, size_t size);
void sarray_sort_cont_(void *A, size_t n, int is_long, buffer *buf,
                       size_t off, size_t size);
void sarray_sort_end_(void *A, size_t n, buffer *buf,
                      size_t align, size_t size);

#define sarray_sort_begin(T,A,n,field,is_long,buf) \
  sarray_sort_begin_(A,n,is_long,buf,ALIGNOF(T),offsetof(T,field),sizeof(T))

#define sarray_sort_cont(T,A,n,field,is_long,buf) \
  sarray_sort_cont_(A,n,is_long,buf,offsetof(T,field),sizeof(T))

#define sarray_sort_end(T,A,n,buf) \
  sarray_sort_end_(A,n,buf,ALIGNOF(T),sizeof(T))

#define sarray_sort(T,A,n, field,is_long, buf) \
  (sarray_sort_begin(T,A,n,field,is_long,buf), sarray_sort_end(T,A,n,buf))

#define sarray_sort_two(T,A,n, field1,is_long1, field2,is_long2, buf) \
  (sarray_sort_begin(T,A,n,field2,is_long2,buf), \
   sarray_sort_cont(T,A,n,field1,is_long1,buf), \
   sarray_sort_end(T,A,n,buf))

#endif
