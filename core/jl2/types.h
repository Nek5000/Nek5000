#ifndef TYPES_H
#define TYPES_H

#if defined(USE_LONG_LONG) || defined(GLOBAL_LONG_LONG)
typedef long long long_long;
#  define WHEN_LONG_LONG(x) x
#  if !defined(LLONG_MAX)
#    if defined(LONG_LONG_MAX)
#      define LLONG_MAX LONG_LONG_MAX
#    else
#      define LLONG_MAX 9223372036854775807
#    endif
#  endif
#  if !defined(LLONG_MIN)
#    if defined(LONG_LONG_MIN)
#      define LLONG_MIN LONG_LONG_MIN
#    else
#      define LLONG_MIN -9223372036854775807
#    endif
#  endif
#else
#  define WHEN_LONG_LONG(x)
#endif

#if !defined(USE_LONG) && !defined(USE_LONG_LONG)
#  define TYPE_LOCAL(i,l,ll) i
#elif defined(USE_LONG)
#  define TYPE_LOCAL(i,l,ll) l
#elif defined(USE_LONG_LONG)
#  define TYPE_LOCAL(i,l,ll) ll
#endif

#if !defined(GLOBAL_LONG) && !defined(GLOBAL_LONG_LONG)
#  define TYPE_GLOBAL(i,l,ll) i
#elif defined(GLOBAL_LONG)
#  define TYPE_GLOBAL(i,l,ll) l
#else
#  define TYPE_GLOBAL(i,l,ll) ll
#endif

/* local integer type: for quantities O(N/P) */
#define sint   signed TYPE_LOCAL(int,long,long long)
#define uint unsigned TYPE_LOCAL(int,long,long long)
#define iabs TYPE_LOCAL(abs,labs,llabs)

/* global integer type: for quantities O(N) */
#define slong   signed TYPE_GLOBAL(int,long,long long)
#define ulong unsigned TYPE_GLOBAL(int,long,long long)
#define iabsl TYPE_GLOBAL(abs,labs,llabs)

#endif

