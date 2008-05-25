#include <limits.h>
#include <string.h>

#include "types.h"

typedef uint Index;

#define sort jl_sort
#define Value uint
#define Data sort_data
typedef struct { Value v; Index i; } Data;
#include "sort_imp.c"

#undef Value
#undef Data

#ifdef GLOBAL_INT
#  define Value ulong
#  define Data sort_data_long
   typedef struct { Value v; Index i; } Data;
#  define radix_count         radix_count_long
#  define radix_offsets       radix_offsets_long
#  define radix_zeros         radix_zeros_long
#  define radix_pass          radix_pass_long
#  define radix_sort          radix_sort_long
#  define radix_index_pass_b  radix_index_pass_b_long
#  define radix_index_pass_m  radix_index_pass_m_long
#  define radix_index_pass_e  radix_index_pass_e_long
#  define radix_index_pass_be radix_index_pass_be_long
#  define radix_index_sort    radix_index_sort_long
#  define merge_sort          merge_sort_long
#  define merge_index_sort    merge_index_sort_long
#  define sort                sort_long
#  define index_sort          index_sort_long
#  include "sort_imp.c"
#endif

