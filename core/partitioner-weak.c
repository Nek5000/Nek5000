#include "partitioner.h"

static int exit_if_not_enabled(const char *name, const char *pplist,
                               MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank > 0) goto barrier_before_abort;
  fprintf(stderr,
          "Error: %s is not enabled during the build. Use "
          "\"%s\" in PPLIST to enable it.",
          name, pplist);
  fflush(stderr);

barrier_before_abort:
  MPI_Barrier(comm);
  MPI_Abort(comm, EXIT_FAILURE);
}

int parMETIS_partMesh(int *part, long long *vl, int nel, int nv, double *opt,
                      MPI_Comm comm) __attribute__((weak));

int parMETIS_partMesh(int *part, long long *vl, int nel, int nv, double *opt,
                      MPI_Comm comm) {
  exit_if_not_enabled("parMETIS", "PARMETIS", comm);
}

int Zoltan2_partMesh(int *part, long long *vl, unsigned nel, int nv,
                     double *opt, MPI_Comm comm) __attribute__((weak));

int Zoltan2_partMesh(int *part, long long *vl, unsigned nel, int nv,
                     double *opt, MPI_Comm comm) {
  exit_if_not_enabled("Zoltan2", "ZOLTAN2", comm);
}
