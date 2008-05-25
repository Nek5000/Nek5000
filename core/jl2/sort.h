#ifndef SORT_H
#define SORT_H

#if !defined(ERRMEM_H) || !defined(TYPES_H)
#warning "sort.h" requires "errmem.h" and "types.h"
#endif

/*------------------------------------------------------------------------------
  
  Sort
  
  O(n) stable sort with good performance for all n

  sortv     (uint  *out,  const uint  *A, uint n, uint stride,
             void *work, size_t work_off)
  sortv_long(ulong *out,  const ulong *A, uint n, uint stride,
             void *work, size_t work_off)

  sortp     (uint *perm, int perm_start,
             const uint  *A, uint n, uint stride,
             void *work, size_t work_off)
  sortp_long(uint *perm, int perm_start,
             const ulong *A, uint n, uint stride,
             void *work, size_t work_off)

  A, n, stride : specifices the input (stride is in bytes!)
  perm_start: nonzero indicates that perm specifies the starting
              permutation, otherwise the starting permutation is the identity;
              the sort is stable w.r.t. the starting permutation
  work, work_off:  scratch area begins at (char*)work + work_off
  
  each routine X has an accompanying function
    X_worksize(uint n, size_t work_off)
  that gives the minimum size of the work area in bytes (including work_off)

  the sortv routines return the sorted values
  the sortp routines return a permutation

  example, demonstrating composability:
  
    uint v[N][M]; // or ulong v...
    uint perm[N];
    buffer buf;
    ...
    buffer_init(&buf, sortp_worksize(N,0));
    sortp(perm,0, &v[0][key2],N,sizeof(uint[M]), buf.ptr,0);
    sortp(perm,1, &v[0][key1],N,sizeof(uint[M]), buf.ptr,0);
    buffer_free(&buf);
    
    now the array can be accessed in sorted order as
       (v[perm[0]][key1], v[perm[0]][key2])
       (v[perm[1]][key1], v[perm[1]][key2])
       ...

  ----------------------------------------------------------------------------*/

#ifdef PREFIX
#  define sortv_ui TOKEN_PASTE(PREFIX,sortv_ui)
#  define sortv_ul TOKEN_PASTE(PREFIX,sortv_ul)
#  define sortp_ui TOKEN_PASTE(PREFIX,sortp_ui)
#  define sortp_ul TOKEN_PASTE(PREFIX,sortp_ul)
#endif

#ifndef USE_LONG
#  define sortv sortv_ui
#  define sortp sortp_ui
#else
#  define sortv sortv_ul
#  define sortp sortp_ul
#endif

#ifndef GLOBAL_LONG
#  define sortv_long sortv
#  define sortp_long sortp
#else
#  define sortv_long sortv_ul
#  define sortp_long sortp_ul
#endif

static size_t sortv_worksize(uint n, size_t work_off)
{
  return align_as(uint, work_off+n*sizeof(uint));
}

static size_t sortv_long_worksize(uint n, size_t work_off)
{
  return align_as(ulong, work_off+n*sizeof(ulong));
}

static size_t sortp_worksize(uint n, size_t work_off)
{
  typedef struct { uint v; uint i; } sortp_data;
  return align_as(sortp_data, work_off+n*2*sizeof(sortp_data));
}

static size_t sortp_long_worksize(uint n, size_t work_off)
{
  typedef struct { ulong v; uint i; } sortp_long_data;
  return align_as(sortp_long_data, work_off+n*2*sizeof(sortp_long_data));
}

void sortv_ui(unsigned *out, const unsigned *A, uint n, unsigned stride,
              buffer *buf, int resize);
void sortv_ul(unsigned long *out,
              const unsigned long *A, uint n, unsigned stride,
              buffer *buf, int resize);
void sortp_ui(uint *out, int start_perm,
              const unsigned *A, uint n, unsigned stride,
              buffer *buf, int resize);
void sortp_ul(uint *out, int start_perm,
              const unsigned long *A, uint n, unsigned stride,
              buffer *buf, int resize);

#endif
