#include <stdbool.h>

#include "parhip_interface.h"
#include "partitioner.h"

int parHIP_partMesh(int *part, long long *vl, int nel, int nv, double *opt,
                    MPI_Comm ce) {
  sint ierrm = 1;
  sint ibfr;

  if (sizeof(idxtype) != sizeof(unsigned long long)) {
    fprintf(stderr, "ERROR: invalid sizeof(idxtype)!\n");
    goto wait_and_check;
  }

  int    verbose       = (int)opt[1];
  int    num_parts     = (int)opt[2];
  double imbalance_tol = opt[3];

  struct comm comm;
  MPI_Comm    active;

  int color = (nel > 0) ? 1 : MPI_UNDEFINED;
  MPI_Comm_split(ce, color, 0, &active);
  if (color == MPI_UNDEFINED) goto wait_and_check;

  comm_init(&comm, active);
  if (comm.id == 0 && verbose) printf("Running parHIP ... "), fflush(stdout);

  idxtype *nel_array = tcalloc(idxtype, comm.np);
  idxtype  nel_ull   = nel;
  MPI_Allgather(&nel_ull, 1, MPI_UNSIGNED_LONG_LONG, nel_array, 1,
                MPI_UNSIGNED_LONG_LONG, comm.c);

  idxtype *vtxdist = tcalloc(idxtype, comm.np + 1);
  for (int i = 0; i < comm.np; i++) vtxdist[i + 1] = vtxdist[i] + nel_array[i];
  free(nel_array);

  assert((vtxdist[comm.id + 1] - vtxdist[comm.id]) == nel);

  graph_t *graph = graph_create(vl, nel, nv, &comm);
  assert(graph->num_vertices == nel);

  idxtype *xadj = tcalloc(idxtype, nel + 1);
  for (uint i = 0; i < nel + 1; i++) xadj[i] = graph->neighbor_index[i];

  uint     num_neighbors = graph->neighbor_index[nel];
  idxtype *adjncy        = tcalloc(idxtype, num_neighbors);
  for (uint i = 0; i < num_neighbors; i++)
    adjncy[i] = graph->neighbor_ids[i] - 1;

  idxtype *vwgt = tcalloc(idxtype, nel);
  for (uint i = 0; i < nel; i++) vwgt[i] = 1;

  idxtype *adjwgt = tcalloc(idxtype, num_neighbors);
  for (uint i = 0; i < num_neighbors; i++)
    adjwgt[i] = graph->neighbor_weights[i];

  graph_destroy(&graph);

  bool     suppress_output = false;
  int      seed            = 0;
  int      mode            = ECOMESH;
  int      edgecut         = 0;
  idxtype *part_           = tcalloc(idxtype, nel);

  double time0 = comm_time();
  ParHIPPartitionKWay(vtxdist, xadj, adjncy, vwgt, adjwgt, &num_parts,
                      &imbalance_tol, suppress_output, seed, mode, &edgecut,
                      part_, &active);

  double time = comm_time() - time0;
  if (comm.id == 0 && verbose) printf("%lf sec\n", time), fflush(stdout);

  free(vtxdist), free(xadj), free(adjncy), free(vwgt), free(adjwgt);
  comm_free(&comm);
  MPI_Comm_free(&active);

  for (uint i = 0; i < nel; i++) part[i] = part_[i];
  free(part_);
  ierrm = 0;

wait_and_check:
  fflush(stderr);
  comm_init(&comm, ce);
  comm_allreduce(&comm, gs_int, gs_max, &ierrm, 1, &ibfr);
  comm_free(&comm);
  return (ierrm != 0);
}
