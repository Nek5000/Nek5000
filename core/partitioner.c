#include "gslib.h"

#if defined(PARRSB)
#include "parRSB.h"
#endif

#include "name.h"

#define MAXNV 8 /* maximum number of vertices per element */

typedef struct {
  long long vtx[MAXNV];
  ulong eid;
  int proc;
} edata;

void print_part_stat(long long *vtx, int nel, int nv, comm_ext ce) {
  int i, j;

  struct comm comm;
  int np, id;

  int Nmsg;
  int *Ncomm;

  int nelMin, nelMax;
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

  if (np == 1)
    return;

  numPoints = nel * nv;
  data = (long long *)malloc(numPoints * sizeof(long long));
  for (i = 0; i < numPoints; i++)
    data[i] = vtx[i];

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

  nsMax = Ncomm[0];
  nsMin = Ncomm[0];
  nsSum = Ncomm[0];
  for (i = 1; i < Nmsg; ++i) {
    nsMax = Ncomm[i] > Ncomm[i - 1] ? Ncomm[i] : Ncomm[i - 1];
    nsMin = Ncomm[i] < Ncomm[i - 1] ? Ncomm[i] : Ncomm[i - 1];
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

  nsSum = nsSum / Nmsg;
  comm_allreduce(&comm, gs_int, gs_add, &nsSum, 1, &b);

  nelMax = nel;
  nelMin = nel;
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);

  if (id == 0) {
    printf(" nElements   max/min/bal: %d %d %.2f\n", nelMax, nelMin,
           (double)nelMax / nelMin);
    printf(" nMessages   max/min/avg: %d %d %.2f\n", ncMax, ncMin,
           (double)ncSum / np);
    printf(" msgSize     max/min/avg: %d %d %.2f\n", nsMax, nsMin,
           (double)nsSum / np);
    printf(" msgSizeSum  max/min/avg: %d %d %.2f\n", nssMax, nssMin,
           (double)nssSum / np);
    fflush(stdout);
  }

  comm_free(&comm);
}

int redistribute_data(int *nel_, long long *vl, long long *el, int *part,
                      int nv, int lelt, struct comm *comm) {
  int nel = *nel_;

  struct crystal cr;
  struct array eList;
  edata *data;

  int count, e, n, ibuf;

  /* redistribute data */
  array_init(edata, &eList, nel), eList.n = nel;
  for (data = eList.ptr, e = 0; e < nel; ++e) {
    data[e].proc = part[e];
    data[e].eid = el[e];
    for (n = 0; n < nv; ++n) {
      data[e].vtx[n] = vl[e * nv + n];
    }
  }

  crystal_init(&cr, comm);
  sarray_transfer(edata, &eList, proc, 0, &cr);
  crystal_free(&cr);

  buffer bfr;
  buffer_init(&bfr, 1024);
  sarray_sort(edata, eList.ptr, eList.n, eid, 1, &bfr);
  buffer_free(&bfr);

  *nel_ = nel = eList.n;
  count = 0;
  if (nel > lelt)
    count = 1;
  comm_allreduce(comm, gs_int, gs_add, &count, 1, &ibuf);
  if (count > 0) {
    count = nel;
    comm_allreduce(comm, gs_int, gs_max, &count, 1, &ibuf);
    if (comm->id == 0)
      printf("ERROR: resulting parition requires lelt=%d!\n", count);
    return 1;
  }

  for (data = eList.ptr, e = 0; e < nel; ++e) {
    el[e] = data[e].eid;
    for (n = 0; n < nv; ++n) {
      vl[e * nv + n] = data[e].vtx[n];
    }
  }

  array_free(&eList);

  return 0;
}

#define fpartmesh FORTRAN_UNPREFIXED(fpartmesh, FPARTMESH)

void fpartmesh(int *nell, long long *el, long long *vl, double *xyz,
               const int *const lelm, const int *const nve,
               const int *const fcomm, const int *const fpartitioner,
               const int *const falgo, const int *const loglevel, int *rtval) {
  struct comm comm;

  int nel, nv, lelt, partitioner, algo;
  int e, n;
  int count, ierr, ibuf;
  int *part;
  int opt[3];

#if defined(PARRSB)
  lelt = *lelm;
  nel = *nell;
  nv = *nve;
  partitioner = *fpartitioner;
  algo = *falgo;

#if defined(MPI)
  comm_ext cext = MPI_Comm_f2c(*fcomm);
#else
  comm_ext cext = 0;
#endif
  comm_init(&comm, cext);

  ierr = 1;
  part = (int *)malloc(lelt * sizeof(int));

  parrsb_options options = parrsb_default_options;
  options.partitioner = partitioner;
  if (partitioner == 0) // RSB
    options.rsb_algo = algo;

  if (*loglevel > 2)
    print_part_stat(vl, nel, nv, cext);

  ierr = parrsb_part_mesh(part, vl, xyz, NULL, nel, nv, &options, comm.c);
  if (ierr != 0)
    goto err;

  ierr = redistribute_data(&nel, vl, el, part, nv, lelt, &comm);
  if (ierr != 0)
    goto err;

  if (*loglevel > 2)
    print_part_stat(vl, nel, nv, cext);

  free(part);
  comm_free(&comm);

  *nell = nel;
  *rtval = 0;
  return;

err:
  fflush(stdout);
  *rtval = 1;
#endif
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
  comm_ext cext = MPI_Comm_f2c(*fcomm);
#else
  comm_ext cext = 0;
#endif
  comm_init(&comm, cext);

  int *const part = (int *)malloc(lelm * sizeof(int));
  parrsb_part_solid(part, vl2, *nel2, vl1, *nel1, *nv, comm.c);

  int ierr = redistribute_data(nel2, vl2, el2, part, *nv, lelm, &comm);
  if (ierr != 0)
    goto err;
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
  comm_ext c = MPI_Comm_f2c(*comm);
#else
  comm_ext c = 0;
#endif

  print_part_stat(vtx, *nel, *nv, c);
}
