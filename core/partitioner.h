#if !defined(__NEK5000_PARTITIONER_H__)
#define __NEK5000_PARTITIONER_H__

#include "gslib.h"

#if defined(PARRSB)
#include "parRSB.h"
#endif

#define MAXNV 8 /* maximum number of vertices per element */

typedef struct {
  long long vtx[MAXNV];
  ulong eid;
  int proc;
} edata;

int parMETIS_partMesh(int *part, long long *vl, int nel, int nv, double *opt,
                      MPI_Comm comm);

int Zoltan2_partMesh(int *part, long long *vl, unsigned nel, int nv,
                     double *opt, MPI_Comm comm);

#endif // __NEK5000_PARTITIONER_H__
