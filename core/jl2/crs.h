#ifndef CRS_H
#define CRS_H

#if !defined(COMM_H)
#warning "crs.h" requires "comm.h"
#endif

#define crs_setup PREFIXED_NAME(crs_setup)
#define crs_solve PREFIXED_NAME(crs_solve)
#define crs_stats PREFIXED_NAME(crs_stats)
#define crs_free  PREFIXED_NAME(crs_free )

typedef struct crs_data_ crs_data;

crs_data *crs_setup(uint n, const ulong *id,
                    uint nz, const uint *Ai, const uint *Aj, const double *A,
                    uint null_space, const struct comm *comm);
void crs_solve(double *x, crs_data *data, const double *b);
void crs_stats(crs_data *data);
void crs_free(crs_data *data);

#endif

