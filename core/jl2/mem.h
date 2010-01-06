#ifndef MEM_H
#define MEM_H

/* requires:
     <stddef.h> for size_t, offsetof
     <stdlib.h> for malloc, calloc, realloc, free
     "fail.h"
*/

#ifndef FAIL_H
#warning "mem.h" requires "fail.h"
#endif

/*--------------------------------------------------------------------------
   Memory Allocation Wrappers to Catch Out-of-memory
  --------------------------------------------------------------------------*/

static void *smalloc(size_t size, const char *file, unsigned line)
{
  void *res = malloc(size);
  if(!res && size)
    fail(1,"%s(%u): allocation of %ld bytes failed\n",file,line,(long)size);
  return res;
}

static void *scalloc(size_t nmemb, size_t size, const char *file, unsigned line)
{
  void *res = calloc(nmemb, size);
  if(!res && nmemb)
    fail(1,"%s(%u): allocation of %ld bytes failed\n",file,line,
           (long)size*nmemb);
  return res;
}

static void *srealloc(void *ptr, size_t size, const char *file, unsigned line)
{
  void *res = realloc(ptr, size);
  if(!res && size)
    fail(1,"%s(%u): allocation of %ld bytes failed\n",file,line,(long)size);
  return res;
}

#define tmalloc(type, count) \
  ((type*) smalloc((count)*sizeof(type),__FILE__,__LINE__) )
#define tcalloc(type, count) \
  ((type*) scalloc((count),sizeof(type),__FILE__,__LINE__) )
#define trealloc(type, ptr, count) \
  ((type*) srealloc((ptr),(count)*sizeof(type),__FILE__,__LINE__) )

/*--------------------------------------------------------------------------
   A dynamic array
  --------------------------------------------------------------------------*/
typedef struct { void *ptr; size_t n,max; } array;
#define null_array {0,0,0}
static void array_init_(array *a, size_t max, size_t size,
                        const char *file, unsigned line)
{
  a->n=0, a->max=max, a->ptr=smalloc(max*size,file,line);
}
static void array_resize_(array *a, size_t max, size_t size,
                          const char *file, unsigned line)
{
  a->max=max, a->ptr=srealloc(a->ptr,max*size,file,line);
}
static void *array_reserve_(array *a, size_t min, size_t size,
                            const char *file, unsigned line)
{
  size_t max = a->max;
  if(max<min) {
    max+=max/2+1;
    if(max<min) max=min;
    array_resize_(a,max,size,file,line);
  }
  return a->ptr;
}

static void array_free(array *a) { free(a->ptr); }
#define array_init(T,a,max) array_init_(a,max,sizeof(T),__FILE__,__LINE__)
#define array_resize(T,a,max) array_resize_(a,max,sizeof(T),__FILE__,__LINE__)
#define array_reserve(T,a,min) array_reserve_(a,min,sizeof(T),__FILE__,__LINE__)

/*--------------------------------------------------------------------------
   Buffer = char array
  --------------------------------------------------------------------------*/
typedef array buffer;
#define null_buffer null_array
#define buffer_init(b,max) array_init(char,b,max)
#define buffer_resize(b,max) array_resize(char,b,max)
#define buffer_reserve(b,max) array_reserve(char,b,max)
#define buffer_free(b) array_free(b)

/*--------------------------------------------------------------------------
   Alignment routines
  --------------------------------------------------------------------------*/
#define ALIGNOF(T) offsetof(struct { char c; T x; }, x)
static size_t align_as_(size_t a, size_t n) { return (n+a-1)/a*a; }
#define align_as(T,n) align_as_(ALIGNOF(T),n)
#define align_ptr(T,base,offset) ((T*)((char*)(base)+align_as(T,offset)))
#endif

