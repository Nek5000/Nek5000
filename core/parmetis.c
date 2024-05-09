#include <stddef.h>

#include "defs.h"
#include "parmetis.h"
#include "partitioner.h"

int parMETIS_partMesh(int *part, long long *vl, int nel, int nv, int *opt,
                      MPI_Comm ce) {
  int    i, j;
  int    ierrm;
  double time, time0;

  MPI_Comm    comms;
  struct comm comm;
  int         color;
  int         ibuf;

  struct crystal cr;
  struct array   A;
  edata         *row;

  long long  nell;
  long long *nelarray;
  idx_t     *elmdist;
  idx_t     *evlptr;
  idx_t     *part_;
  real_t    *tpwgts;
  idx_t      edgecut;
  real_t     ubvec;
  idx_t     *elmwgt;
  idx_t      wgtflag;
  idx_t      numflag;
  idx_t      ncon;
  idx_t      ncommonnodes;
  idx_t      nparts;
  idx_t      nelsm;
  idx_t      options[10];

  ierrm        = METIS_OK;
  nell         = nel;
  edgecut      = 0;
  wgtflag      = 0;
  numflag      = 0;
  ncon         = 1;
  ubvec        = 1.02;
  elmwgt       = NULL; /* no weights */
  ncommonnodes = 2;

  part_ = (idx_t *)malloc(nel * sizeof(idx_t));

  if (sizeof(idx_t) != sizeof(long long)) {
    printf("ERROR: invalid sizeof(idx_t)!\n");
    goto err;
  }
  if (nv != 4 && nv != 8) {
    printf("ERROR: nv is %d but only 4 and 8 are supported!\n", nv);
    goto err;
  }

  color = MPI_UNDEFINED;
  if (nel > 0) color = 1;
  MPI_Comm_split(ce, color, 0, &comms);
  if (color == MPI_UNDEFINED) goto end;

  comm_init(&comm, comms);
  if (comm.id == 0) printf("Running parMETIS ... "), fflush(stdout);

  nelarray = (long long *)malloc(comm.np * sizeof(long long));
  MPI_Allgather(&nell, 1, MPI_LONG_LONG_INT, nelarray, 1, MPI_LONG_LONG_INT,
                comm.c);
  elmdist    = (idx_t *)malloc((comm.np + 1) * sizeof(idx_t));
  elmdist[0] = 0;
  for (i = 0; i < comm.np; ++i)
    elmdist[i + 1] = elmdist[i] + (idx_t)nelarray[i];
  free(nelarray);

  evlptr    = (idx_t *)malloc((nel + 1) * sizeof(idx_t));
  evlptr[0] = 0;
  for (i = 0; i < nel; ++i) evlptr[i + 1] = evlptr[i] + nv;
  nelsm = elmdist[comm.id + 1] - elmdist[comm.id];
  evlptr[nelsm]--;

  if (nv == 8) ncommonnodes = 4;
  nparts = comm.np;

  options[0]                  = 1;
  options[PMV3_OPTION_DBGLVL] = 0;
  options[PMV3_OPTION_SEED]   = 0;
  if (opt[0] != 0) {
    options[PMV3_OPTION_DBGLVL] = opt[1];
    if (opt[2] != 0) {
      options[3] = PARMETIS_PSR_UNCOUPLED;
      nparts     = opt[2];
    }
  }

  tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
  for (i = 0; i < ncon * nparts; ++i) tpwgts[i] = 1. / (real_t)nparts;

  if (options[3] == PARMETIS_PSR_UNCOUPLED)
    for (i = 0; i < nel; ++i) part_[i] = comm.id;

  comm_barrier(&comm);
  time0 = comm_time();
  ierrm =
      ParMETIS_V3_PartMeshKway(elmdist, evlptr, (idx_t *)vl, elmwgt, &wgtflag,
                               &numflag, &ncon, &ncommonnodes, &nparts, tpwgts,
                               &ubvec, options, &edgecut, part_, &comm.c);

  time = comm_time() - time0;
  if (comm.id == 0) printf("%lf sec\n", time), fflush(stdout);

  for (i = 0; i < nel; ++i) part[i] = part_[i];

  free(elmdist);
  free(evlptr);
  free(tpwgts);
  MPI_Comm_free(&comms);
  comm_free(&comm);

end:
  comm_init(&comm, ce);
  comm_allreduce(&comm, gs_int, gs_min, &ierrm, 1, &ibuf);
  if (ierrm != METIS_OK) goto err;
  return 0;

err:
  return 1;
}
