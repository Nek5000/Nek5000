#if !defined(__NEK5000_PARTITIONER_H__)
#define __NEK5000_PARTITIONER_H__

#include "gslib.h"

#if defined(PARRSB)
#include "parRSB.h"
#endif

#define MAXNV 8 /* maximum number of vertices per element */

typedef struct {
  long long vtx[MAXNV];
  ulong     eid;
  int       proc;
} edata;

int parMETIS_partMesh(int *part, long long *vl, int nel, int nv, int *opt,
                      MPI_Comm ce);

int Zoltan_partMesh(int *part, long long *vl, int nel, int nv, MPI_Comm comm,
                    int verbose);

int Zoltan2_partMesh(int *part, long long *vl, unsigned nel, int nv,
                     MPI_Comm comm, int verbose);

int parHIP_partMesh(int *part, long long *vl, int nel, int nv, MPI_Comm comm,
                    int verbose);

typedef struct {
  uint   num_vertices;
  ulong *vertex_ids;
  int   *neighbor_index;
  ulong *neighbor_ids;
  int   *neighbor_procs;
  int   *neighbor_weights;
} graph_t;

graph_t *graph_create(long long *vl, int nelt, int nv, struct comm *c);

void graph_print(graph_t *graph);

void graph_destroy(graph_t **graph);

#endif // __NEK5000_PARTITIONER_H__
