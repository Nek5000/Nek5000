#ifndef XXT_H
#define XXT_H

#if !defined(COMM_H)
#warning "xxt.h" requires "comm.h"
#endif

#ifdef PREFIX
#  define xxt_setup TOKEN_PASTE(PREFIX,xxt_setup)
#  define xxt_solve TOKEN_PASTE(PREFIX,xxt_solve)
#  define xxt_stats TOKEN_PASTE(PREFIX,xxt_stats)
#  define xxt_free  TOKEN_PASTE(PREFIX,xxt_free )
#endif

typedef struct xxt_data_ xxt_data;

xxt_data *xxt_setup(uint n, const ulong *id,
                    uint nz, const uint *Ai, const uint *Aj, const double *A,
                    uint null_space, comm_ext_t comm, uint np);
void xxt_solve(double *x, xxt_data *data, const double *b);
void xxt_stats(xxt_data *data);
void xxt_free(xxt_data *data);

#endif

