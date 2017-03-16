#ifndef XXT_H
#define XXT_H

/* requires "types.h", and, when MPI is defined, "crystal.h" */
#if !defined(TYPES_H) || ( defined(MPI) && !defined(CRYSTAL_H) )
#warning "xxt.h" requires "types.h" and "crystal.h"
#endif

typedef struct xxt_data_ xxt_data;

#ifndef MPI
#  define crystal_data void
#endif

#define xxt_free xxt_jl_free
#define xxt_solve xxt_jl_solve
#define xxt_stats xxt_jl_stats

xxt_data *xxt_setup(uint n, const ulong *id,
                    uint nz, const uint *Ai, const uint *Aj, const real *A,
                    uint null_space, crystal_data *crystal);
void xxt_solve(real *x, xxt_data *data, const real *b);
void xxt_stats(xxt_data *data);
void xxt_free(xxt_data *data);

#ifndef MPI
#  undef crystal_data
#endif

#endif

