#include "partitioner.h"

static void print_part_stat(long long *vtx, int nel, int nv, MPI_Comm ce) {
  int i, j;

  struct comm comm;
  int np, id;

  int Nmsg;
  int *Ncomm;

  int nelMin, nelMax;
  long long nelSum;
  int ncMin, ncMax, ncSum;
  int nsMin, nsMax, nsSum;
  int nssMin, nssMax;
  long long nssSum;

  struct gs_data *gsh;
  int b;
  long long b_long_long;

  int numPoints;
  long long *data;

  comm_init(&comm, ce);
  np = comm.np;
  id = comm.id;

  if (np == 1) return;

  numPoints = nel * nv;
  data = (long long *)malloc(numPoints * sizeof(long long));
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
  nsMin = INT_MAX;
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

  if (Nmsg > 0)
    nsSum = nsSum / Nmsg;
  else
    nsSum = 0;
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
  else
    printf(" nElements   max/min/avg: %d %d INFINITY\n", nelMax, nelMin);
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

  int e, n;
  edata *data;
  for (data = eList.ptr, e = 0; e < nel; ++e) {
    data[e].proc = part[e];
    data[e].eid = el[e];
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

#define check_error(error)                                                     \
  {                                                                            \
    if (error != 0) goto check_global_error;                                   \
  }

#define fpartmesh FORTRAN_UNPREFIXED(fpartmesh, FPARTMESH)
void fpartmesh(int *nell, long long *el, long long *vl, double *xyz,
               const int *const lelm, const int *const nve,
               const int *const fcomm, const int *const fpartitioner,
               const int *const falgo, const int *const loglevel, int *rtval) {
  int nel = *nell;
  int nv = *nve;
  int lelt = *lelm;
  int partitioner = *fpartitioner;
  int algo = *falgo;
  int verbose = *loglevel;
  sint ierr = 1;

  if (nv != 4 && nv != 8) {
    fprintf(stderr, "ERROR: nv is %d but only 4 and 8 are supported!\n", nv);
    goto check_global_error;
  }

  struct comm comm;
#if defined(MPI)
  MPI_Comm cext = MPI_Comm_f2c(*fcomm);
#else
  int cext = 0;
#endif
  comm_init(&comm, cext);

  if (verbose >= 2) print_part_stat(vl, nel, nv, cext);

  double opt[10] = {0};
  opt[0] = 1;
  opt[1] = 0;       /* verbosity */
  opt[2] = comm.np; /* number of partitions */
  opt[3] = 1.05;    /* imbalance tolerance */

  int *part = (int *)malloc(lelt * sizeof(int));
  if (partitioner == 0 || partitioner == 1) {
    opt[0] = partitioner;
    opt[1] = algo;
    ierr = parRSB_partMesh(part, vl, nel, nv, opt, comm.c, xyz);
  } else if (partitioner == 8) {
    ierr = parMETIS_partMesh(part, vl, nel, nv, opt, comm.c);
  } else if (partitioner == 16) {
    ierr = Zoltan2_partMesh(part, vl, nel, nv, opt, comm.c);
  } else if (partitioner == 32) {
    ierr = Zoltan_partMesh(part, vl, nel, nv, opt, comm.c);
  } else if (partitioner == 64) {
    ierr = parHIP_partMesh(part, vl, nel, nv, opt, comm.c);
  }
  check_error(ierr);

  ierr = redistribute_data(&nel, vl, el, part, nv, lelt, &comm);
  check_error(ierr);
  *nell = nel;

  if (verbose >= 2) print_part_stat(vl, nel, nv, cext);

  free(part);

  sint b;
check_global_error:
  fflush(stderr);
  comm_allreduce(&comm, gs_int, gs_max, &ierr, 1, &b);
  comm_free(&comm);
  *rtval = ierr;
}

#define fprintpartstat FORTRAN_UNPREFIXED(printpartstat, PRINTPARTSTAT)
void fprintpartstat(long long *vtx, int *nel, int *nv, int *comm) {

#if defined(MPI)
  MPI_Comm c = MPI_Comm_f2c(*comm);
#else
  MPI_Comm c = 0;
#endif

  print_part_stat(vtx, *nel, *nv, c);
}

#undef fprintpartstat
#undef fpartmesh_greedy
#undef fpartmesh
#undef check_error
