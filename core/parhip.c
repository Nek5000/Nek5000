#include <stdbool.h>

#include "parhip_interface.h"
#include "partitioner.h"

int parHIP_partMesh(int *part, long long *vl, int nel, int nv, MPI_Comm ce,
                    int verbose) {
  if (sizeof(idxtype) != sizeof(unsigned long long)) {
    printf("ERROR: invalid sizeof(idxtype)!\n");
    goto err;
  }
  if (nv != 4 && nv != 8) {
    printf("ERROR: nv is %d but only 4 and 8 are supported!\n", nv);
    goto err;
  }

  struct comm comm;
  MPI_Comm    active;
  sint        ierr = 0;
  {
    int color = (nel > 0) ? 1 : MPI_UNDEFINED;
    MPI_Comm_split(ce, color, 0, &active);
    if (color == MPI_UNDEFINED) goto end;
    comm_init(&comm, active);
  }

  if (comm.id == 0 && verbose)
    fprintf(stderr, "Running parHIP ... "), fflush(stderr);

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

  int      num_parts       = comm.np;
  double   imbalance_tol   = 1.0 / (vtxdist[comm.np] / comm.np);
  bool     suppress_output = false;
  int      seed            = 0;
  int      mode            = ECOMESH;
  int      edgecut         = 0;
  idxtype *part_           = tcalloc(idxtype, nel);

  if (comm.id == 0 && verbose) {
    fprintf(stderr, "imbalance_tol = %lf", imbalance_tol);
    fflush(stderr);
  }

  ParHIPPartitionKWay(vtxdist, xadj, adjncy, vwgt, adjwgt, &num_parts,
                      &imbalance_tol, suppress_output, seed, mode, &edgecut,
                      part_, &active);

  free(vtxdist), free(xadj), free(adjncy), free(vwgt), free(adjwgt);
  comm_free(&comm);
  MPI_Comm_free(&active);

  for (uint i = 0; i < nel; i++) part[i] = part_[i];
  free(part_);

end:
  comm_init(&comm, ce);
  sint ierr_max = ierr, ibfr;
  comm_allreduce(&comm, gs_int, gs_max, &ierr_max, 1, &ibfr);
  return (ierr_max != 0);

err:
  return 1;
}
