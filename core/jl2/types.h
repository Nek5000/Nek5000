#ifndef TYPES_H
#define TYPES_H

/* apparently uint and ulong can be defined already in standard headers */
#define uint uint_
#define ulong ulong_
#define sint sint_
#define slong slong_

/* local integer type: for quantities O(N/P) */
#if !defined(USE_LONG)
   typedef   signed int sint;
   typedef unsigned int uint;
#  define iabs abs
#  define SINT_MIN INT_MIN
#  define SINT_MAX INT_MAX
#else
   typedef   signed long sint;
   typedef unsigned long uint;
#  define iabs labs
#  define SINT_MIN LONG_MIN
#  define SINT_MAX LONG_MAX
#endif

/* global integer type: for quantities O(N) */
#if !defined(GLOBAL_LONG)
   typedef sint slong;
   typedef uint ulong;
#  define iabsl iabs
#  define SLONG_MIN SINT_MIN
#  define SLONG_MAX SINT_MAX
#else
   typedef   signed long slong;
   typedef unsigned long ulong;
#  define iabsl labs
#  define SLONG_MIN LONG_MIN
#  define SLONG_MAX LONG_MAX
#endif

#endif

