#ifndef CRS_H
#define CRS_H

#if !defined(COMM_H)
#warning "crs.h" requires "comm.h"
#endif

#ifdef PREFIX
#  define crs_setup TOKEN_PASTE(PREFIX,crs_setup)
#  define crs_solve TOKEN_PASTE(PREFIX,crs_solve)
#  define crs_stats TOKEN_PASTE(PREFIX,crs_stats)
#  define crs_free  TOKEN_PASTE(PREFIX,crs_free )
#endif

typedef struct crs_data_ crs_data;

crs_data *crs_setup(uint n, const ulong *id,
                    uint nz, const uint *Ai, const uint *Aj, const double *A,
                    uint null_space, const comm_t *comm);
void crs_solve(double *x, crs_data *data, const double *b);
void crs_stats(crs_data *data);
void crs_free(crs_data *data);

#endif

