#include "parRSB.h"

int parRSB_partMesh(int *part, long long *vl, int nel, int nv, double *opt,
                      MPI_Comm comm, double *xyz) {
  parrsb_options options = parrsb_default_options;

  int partitioner = (int)opt[0];
  options.partitioner = partitioner;

  if (partitioner == 0) // If the partitioner is RSB
    options.rsb_algo = (int)opt[1];

  return parrsb_part_mesh(part, vl, xyz, NULL, nel, nv, &options, comm);
}
