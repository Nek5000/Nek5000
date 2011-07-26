#ifndef SORT_H
#define SORT_H

/* requires "types.h" */
#ifndef TYPES_H
#warning "sort.h" requires "types.h"
#endif

/*------------------------------------------------------------------------------
  
  Sort
  
  O(n) stable sort with good performance for all n

  A, n, stride : specifices the input
  
  sort:
     Value out[n] : the sorted values (output)
     Value work[n]: scratch area
  
  index_sort:
     uint idx[n]   : the sorted indices (output)
     Data work[2*n]: scratch area

  example:
  
    uint b[N][M];                      or ulong b...
    sort_data work[2*N];               or sort_data_long work...
    uint p[N];
    ...
    index_sort(&b[0][key],N,M, p, work);   or index_sort_long(...
    
    now the array can be accessed in sorted order as
       b[p[0]][key]
       b[p[1]][key]
       ...

  ----------------------------------------------------------------------------*/

#define sort jl_sort
void sort(const uint *A, uint n, uint stride, uint *out, uint *work);

typedef struct { uint v; uint i; } sort_data;
void index_sort(const uint *A, uint n, uint stride,
                uint *idx, sort_data *work);

#ifdef GLOBAL_INT
  void sort_long(const ulong *A, uint n, uint stride, ulong *out, ulong *work);
  typedef struct { ulong v; uint i; } sort_data_long;
  void index_sort_long(const ulong *A, uint n, uint stride,
                       uint *idx, sort_data_long *work);
#else
#  define sort_long       sort
#  define sort_data_long  sort_data
#  define index_sort_long index_sort
#endif

#endif

