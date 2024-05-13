#include "partitioner.h"

#if defined(ZOLTAN2)

extern int Zoltan2_partMesh(int *part, long long *vl, unsigned nel, int nv,
                            MPI_Comm comm_, int verbose);
#endif // ZOLTAN2

#if defined(PARRSB)
#include "parRSB.h"
#endif

#if defined(PARHIP)
#include <stdbool.h>

#include <parhip_interface.h>

static int ParHIP_partMesh(int *part, long long *vl, int nel, int nv,
                           MPI_Comm ce, int verbose) {
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
    fprintf(stderr, "Running ParHIP ... "), fflush(stderr);

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
  if (ierr_max != 0) goto err;
  return 0;

err:
  return 1;
}

#endif // PARHIP

static void print_part_stat(long long *vtx, int nel, int nv, MPI_Comm ce) {
  int i, j;

  struct comm comm;
  int         np, id;

  int  Nmsg;
  int *Ncomm;

  int       nelMin, nelMax;
  long long nelSum;
  int       ncMin, ncMax, ncSum;
  int       nsMin, nsMax, nsSum;
  int       nssMin, nssMax;
  long long nssSum;

  struct gs_data *gsh;
  int             b;
  long long       b_long_long;

  int        numPoints;
  long long *data;

  comm_init(&comm, ce);
  np = comm.np;
  id = comm.id;

  if (np == 1) return;

  numPoints = nel * nv;
  data      = (long long *)malloc(numPoints * sizeof(long long));
  for (i = 0; i < numPoints; i++) data[i] = vtx[i];

  gsh = gs_setup(data, numPoints, &comm, 0, gs_pairwise, 0);

  pw_data_nmsg(gsh, &Nmsg);
  Ncomm = (int *)malloc(Nmsg * sizeof(int));
  pw_data_size(gsh, Ncomm);

  gs_free(gsh);
  free(data);

  ncMax = Nmsg;
  ncMin = Nmsg;
  ncSum = Nmsg;
  comm_allreduce(&comm, gs_int, gs_max, &ncMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &ncMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &ncSum, 1, &b);

  nsMax = nsSum = 0;
  nsMin         = INT_MAX;
  for (i = 0; i < Nmsg; ++i) {
    nsMax = Ncomm[i] > nsMax ? Ncomm[i] : nsMax;
    nsMin = Ncomm[i] < nsMin ? Ncomm[i] : nsMin;
    nsSum += Ncomm[i];
  }
  comm_allreduce(&comm, gs_int, gs_max, &nsMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nsMin, 1, &b);

  nssMin = nsSum;
  nssMax = nsSum;
  nssSum = nsSum;
  comm_allreduce(&comm, gs_int, gs_max, &nssMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nssMin, 1, &b);
  comm_allreduce(&comm, gs_long_long, gs_add, &nssSum, 1, &b_long_long);

  if (Nmsg > 0) nsSum = nsSum / Nmsg;
  else nsSum = 0;
  comm_allreduce(&comm, gs_int, gs_add, &nsSum, 1, &b);

  nelMax = nel;
  nelMin = nel;
  nelSum = nel;
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);
  comm_allreduce(&comm, gs_long_long, gs_add, &nelSum, 1, &b_long_long);

  sint npp = (Nmsg > 0);
  comm_allreduce(&comm, gs_int, gs_add, &npp, 1, &b);

  if (id > 0) goto comm_free_and_exit;
  if (nelMin > 0)
    printf(" nElements   max/min/avg: %d %d %.2f\n", nelMax, nelMin,
           (double)nelMax / nelMin);
  else printf(" nElements   max/min/avg: %d %d INFINITY\n", nelMax, nelMin);
  printf(" nMessages   max/min/avg: %d %d %.2f\n", ncMax, ncMin,
         (double)ncSum / npp);
  printf(" msgSize     max/min/avg: %d %d %.2f\n", nsMax, nsMin,
         (double)nsSum / npp);
  printf(" msgSizeSum  max/min/avg: %d %d %.2f\n", nssMax, nssMin,
         (double)nssSum / npp);
  fflush(stdout);

comm_free_and_exit:
  comm_free(&comm);
}

static int redistribute_data(int *nel_, long long *vl, long long *el, int *part,
                             int nv, int lelt, struct comm *comm) {
  int nel = *nel_;

  struct array eList;
  array_init(edata, &eList, nel), eList.n = nel;

  int    e, n;
  edata *data;
  for (data = eList.ptr, e = 0; e < nel; ++e) {
    data[e].proc = part[e];
    data[e].eid  = el[e];
    for (n = 0; n < nv; ++n) { data[e].vtx[n] = vl[e * nv + n]; }
  }

  struct crystal cr;
  crystal_init(&cr, comm);
  sarray_transfer(edata, &eList, proc, 0, &cr);
  crystal_free(&cr);

  buffer bfr;
  buffer_init(&bfr, 1024);
  sarray_sort(edata, eList.ptr, eList.n, eid, 1, &bfr);
  buffer_free(&bfr);

  *nel_ = nel = eList.n;
  int ibfr, count = (nel > lelt) ? 1 : 0;
  comm_allreduce(comm, gs_int, gs_add, &count, 1, &ibfr);
  if (count == 0) goto success;

  count = nel;
  comm_allreduce(comm, gs_int, gs_max, &count, 1, &ibfr);
  if (comm->id == 0)
    printf("ERROR: resulting parition requires lelt = %d!\n", count);
  return 1;

success:
  for (data = eList.ptr, e = 0; e < nel; ++e) {
    el[e] = data[e].eid;
    for (n = 0; n < nv; ++n) { vl[e * nv + n] = data[e].vtx[n]; }
  }

  array_free(&eList);

  return 0;
}

#define check_error(ierr)                                                      \
  {                                                                            \
    int ierr_ = (ierr);                                                        \
    if (ierr_ != 0) {                                                          \
      *rtval = 1;                                                              \
      return;                                                                  \
    }                                                                          \
  }

#define fpartmesh FORTRAN_UNPREFIXED(fpartmesh, FPARTMESH)
void fpartmesh(int *nell, long long *el, long long *vl, double *xyz,
               const int *const lelm, const int *const nve,
               const int *const fcomm, const int *const fpartitioner,
               const int *const falgo, const int *const loglevel, int *rtval) {
  int nel         = *nell;
  int nv          = *nve;
  int lelt        = *lelm;
  int partitioner = *fpartitioner;
  int algo        = *falgo;

  struct comm comm;
#if defined(MPI)
  MPI_Comm cext = MPI_Comm_f2c(*fcomm);
#else
  MPI_Comm cext = 0;
#endif
  comm_init(&comm, cext);

  int  ierr = 1;
  int *part = (int *)malloc(lelt * sizeof(int));

  if (*loglevel > 2) print_part_stat(vl, nel, nv, cext);

  if (partitioner == 0 || partitioner == 1) {
#if defined(PARRSB)
    parrsb_options options = parrsb_default_options;
    options.partitioner    = partitioner;
    if (partitioner == 0) // RSB
      options.rsb_algo = algo;

    ierr = parrsb_part_mesh(part, vl, xyz, NULL, nel, nv, &options, comm.c);
#endif
  } else if (partitioner == 8) {
#if defined(PARMETIS)
    int opt[3];
    opt[0] = 1;
    opt[1] = 0; /* verbosity */
    opt[2] = comm.np;

    ierr = parMETIS_partMesh(part, vl, nel, nv, opt, comm.c);
#endif
  } else if (partitioner == 16) {
#if defined(ZOLTAN2)
    ierr = Zoltan2_partMesh(part, vl, nel, nv, comm.c, 1);
#endif
  } else if (partitioner == 32) {
#if defined(ZOLTAN)
    ierr = Zoltan_partMesh(part, vl, nel, nv, comm.c, 1);
#endif
  } else if (partitioner == 64) {
#if defined(PARHIP)
    ierr = ParHIP_partMesh(part, vl, nel, nv, comm.c, 1);
#endif
  }

  check_error(ierr);

  ierr = redistribute_data(&nel, vl, el, part, nv, lelt, &comm);
  check_error(ierr);

  if (*loglevel > 2) print_part_stat(vl, nel, nv, cext);

  free(part), comm_free(&comm);
  *nell = nel, *rtval = 0;
}

#define fpartmesh_greedy FORTRAN_UNPREFIXED(fpartmesh_greedy, FPARTMESH_GRREDY)
void fpartmesh_greedy(int *const nel2, long long *const el2,
                      long long *const vl2, const int *const nel1,
                      const long long *const vl1, const int *const lelm_,
                      const int *const nv, const int *const fcomm,
                      int *const rtval) {
#if defined(PARRSB)
  const int lelm = *lelm_;

  struct comm comm;
#if defined(MPI)
  MPI_Comm cext = MPI_Comm_f2c(*fcomm);
#else
  MPI_Comm cext = 0;
#endif
  comm_init(&comm, cext);

  int *const part = (int *)malloc(lelm * sizeof(int));
  parrsb_part_solid(part, vl2, *nel2, vl1, *nel1, *nv, comm.c);

  int ierr = redistribute_data(nel2, vl2, el2, part, *nv, lelm, &comm);
  if (ierr != 0) goto err;
  *rtval = 0;

  free(part);
  comm_free(&comm);
  return;

err:
  fflush(stdout);
  *rtval = 1;
#endif
}

#define fprintpartstat FORTRAN_UNPREFIXED(printpartstat, PRINTPARTSTAT)
void fprintpartstat(long long *vtx, int *nel, int *nv, int *comm) {

#if defined(MPI)
  MPI_Comm c = MPI_Comm_f2c(*comm);
#else
  MPI_Comm c    = 0;
#endif

  print_part_stat(vtx, *nel, *nv, c);
}
