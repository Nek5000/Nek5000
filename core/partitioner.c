#include "gslib.h"

#define MAXNV 8 /* maximum number of vertices per element */

typedef struct {
  long long vtx[MAXNV];
  ulong eid;
  int proc;
} edata;

#if defined(PARRSB)
#include "parRSB.h"
#endif

#if defined(PARMETIS)
#include "defs.h"
#include "parmetis.h"

static int parMETIS_partMesh(int *part, long long *vl, int nel, int nv,
                             int *opt, comm_ext ce) {
  int i, j;
  int ierrm;
  double time, time0;

  MPI_Comm comms;
  struct comm comm;
  int color;
  int ibuf;

  struct crystal cr;
  struct array A;
  edata *row;

  long long nell;
  long long *nelarray;
  idx_t *elmdist;
  idx_t *evlptr;
  idx_t *part_;
  real_t *tpwgts;
  idx_t edgecut;
  real_t ubvec;
  idx_t *elmwgt;
  idx_t wgtflag;
  idx_t numflag;
  idx_t ncon;
  idx_t ncommonnodes;
  idx_t nparts;
  idx_t nelsm;
  idx_t options[10];

  ierrm = METIS_OK;
  nell = nel;
  edgecut = 0;
  wgtflag = 0;
  numflag = 0;
  ncon = 1;
  ubvec = 1.02;
  elmwgt = NULL; /* no weights */
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
  elmdist = (idx_t *)malloc((comm.np + 1) * sizeof(idx_t));
  elmdist[0] = 0;
  for (i = 0; i < comm.np; ++i)
    elmdist[i + 1] = elmdist[i] + (idx_t)nelarray[i];
  free(nelarray);

  evlptr = (idx_t *)malloc((nel + 1) * sizeof(idx_t));
  evlptr[0] = 0;
  for (i = 0; i < nel; ++i) evlptr[i + 1] = evlptr[i] + nv;
  nelsm = elmdist[comm.id + 1] - elmdist[comm.id];
  evlptr[nelsm]--;

  if (nv == 8) ncommonnodes = 4;
  nparts = comm.np;

  options[0] = 1;
  options[PMV3_OPTION_DBGLVL] = 0;
  options[PMV3_OPTION_SEED] = 0;
  if (opt[0] != 0) {
    options[PMV3_OPTION_DBGLVL] = opt[1];
    if (opt[2] != 0) {
      options[3] = PARMETIS_PSR_UNCOUPLED;
      nparts = opt[2];
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
#endif // PARMETIS

#if defined(ZOLTAN)
#include "zoltan.h"

typedef struct {
  uint num_vertices;
  ZOLTAN_ID_TYPE *vertex_ids;
  int *neighbor_index;
  ZOLTAN_ID_TYPE *neighbor_ids;
  int *neighbor_procs;
} graph_t;

static graph_t *graph_create(long long *vl, int nelt, int nv, struct comm *c) {
  // Calculate the total number of elements and the offset of elements
  // for each process.
  slong out[2], wrk[2], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  ulong start = out[0], nelg = out[1];

  // Fill the vertices array with information about the vertices.
  typedef struct {
    ulong vid, eid;
    uint dest, np;
  } vertex_t;

  struct array vertices;
  {
    array_init(vertex_t, &vertices, nelt * nv);
    vertex_t vertex;
    for (uint i = 0; i < nelt; i++) {
      vertex.eid = i + start + 1;
      for (uint j = 0; j < nv; j++) {
        vertex.vid = vl[i * nv + j];
        vertex.dest = vertex.vid % c->np;
        array_cat(vertex_t, &vertices, &vertex, 1);
      }
    }
    // Sanity check.
    assert(vertices.n == nelt * nv);
  }

  // Send vertices to the correct process so we can match them up.
  struct crystal cr;
  {
    crystal_init(&cr, c);
    sarray_transfer(vertex_t, &vertices, dest, 1, &cr);
  }

  buffer bfr;
  {
    buffer_init(&bfr, 1024);
    sarray_sort_2(vertex_t, vertices.ptr, vertices.n, vid, 1, eid, 1, &bfr);
  }

  // Find all the elements shared by a given vertex.
  struct array neighbors;
  {
    array_init(vertex_t, &neighbors, vertices.n);

    const vertex_t *const pv = (const vertex_t *const)vertices.ptr;
    uint vn = vertices.n;
    uint s = 0;
    while (s < vn) {
      uint e = s + 1;
      while (e < vn && pv[s].vid == pv[e].vid) e++;

      for (uint j = s; j < e; j++) {
        vertex_t v = pv[j];
        for (uint k = s; k < e; k++) {
          v.vid = pv[k].eid;
          v.np = pv[k].dest;
          array_cat(vertex_t, &neighbors, &v, 1);
        }
      }
      s = e;
    }
  }
  array_free(&vertices);

  {
    sarray_transfer(vertex_t, &neighbors, dest, 0, &cr);
    sarray_sort_2(vertex_t, neighbors.ptr, neighbors.n, eid, 1, vid, 1, &bfr);
  }
  crystal_free(&cr);
  buffer_free(&bfr);

  // Compre the neighbors array to find the number of neighbors for each
  // element.
  typedef struct {
    ulong nid, eid;
    uint np;
  } neighbor_t;

  struct array compressed;
  {
    array_init(neighbor_t, &compressed, neighbors.n);

    const vertex_t *const pn = (const vertex_t *const)neighbors.ptr;
    uint nn = neighbors.n;
    neighbor_t nbr;
    uint s = 0;
    while (s < nn) {
      uint e = s + 1;
      while (e < nn && pn[s].eid == pn[e].eid && pn[s].vid == pn[e].vid) e++;

      // Sanity check.
      for (uint j = s + 1; j < e; j++) assert(pn[s].np == pn[j].np);

      // Add to compressed array.
      nbr.eid = pn[s].eid;
      nbr.nid = pn[s].vid;
      nbr.np = pn[s].np;
      if (nbr.eid != nbr.nid) array_cat(neighbor_t, &compressed, &nbr, 1);
      s = e;
    }
  }
  array_free(&neighbors);

  graph_t *graph = tcalloc(graph_t, 1);
  graph->num_vertices = nelt;
  graph->vertex_ids = tcalloc(ZOLTAN_ID_TYPE, nelt);
  for (uint i = 0; i < nelt; i++) graph->vertex_ids[i] = i + start + 1;

  graph->neighbor_index = tcalloc(int, nelt + 1);
  {
    neighbor_t *pc = (neighbor_t *)compressed.ptr;
    uint count = 0;
    uint s = 0;
    while (s < compressed.n) {
      uint e = s + 1;
      while (e < compressed.n && pc[s].eid == pc[e].eid) e++;
      count++;
      graph->neighbor_index[count] = e;
      s = e;
    }
    assert(count == nelt);
    assert(graph->neighbor_index[count] == compressed.n);
  }

  graph->neighbor_ids = tcalloc(ZOLTAN_ID_TYPE, compressed.n);
  graph->neighbor_procs = tcalloc(int, compressed.n);
  {
    neighbor_t *pc = (neighbor_t *)compressed.ptr;
    for (uint i = 0; i < compressed.n; i++) {
      graph->neighbor_ids[i] = pc[i].nid;
      graph->neighbor_procs[i] = pc[i].np;
    }
  }

  array_free(&compressed);

  return graph;
}

static void graph_print(graph_t *graph) {
  for (uint i = 0; i < graph->num_vertices; i++) {
    fprintf(stderr, "%lld : ", graph->vertex_ids[i]);
    for (uint j = graph->neighbor_index[i]; j < graph->neighbor_index[i + 1];
         j++)
      fprintf(stderr, " %lld (%d)", graph->neighbor_ids[j],
              graph->neighbor_procs[j]);
    fprintf(stderr, "\n");
    fflush(stderr);
  }
}

static void graph_destroy(graph_t **graph_) {
  graph_t *graph = *graph_;
  free(graph->vertex_ids);
  free(graph->neighbor_index);
  free(graph->neighbor_ids);
  free(graph->neighbor_procs);
  free(graph);
  graph_ = NULL;
}

static int get_number_of_vertices(void *data, int *ierr) {
  graph_t *graph = (graph_t *)data;
  *ierr = ZOLTAN_OK;
  return graph->num_vertices;
}

static void get_vertex_list(void *data, int size_gid, int size_lid,
                            ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                            int wgt_dim, float *obj_wgts, int *ierr) {
  graph_t *graph = (graph_t *)data;
  *ierr = ZOLTAN_FATAL;

  for (int i = 0; i < graph->num_vertices; i++) {
    global_id[i] = graph->vertex_ids[i];
    local_id[i] = i;
  }

  *ierr = ZOLTAN_OK;
}

static void get_edge_size_list(void *data, int size_gid, int size_lid,
                               int num_obj, ZOLTAN_ID_PTR global_id,
                               ZOLTAN_ID_PTR local_id, int *num_edges,
                               int *ierr) {
  graph_t *graph = (graph_t *)data;
  *ierr = ZOLTAN_FATAL;

  if (size_gid != 1 || size_lid != 1 || num_obj != graph->num_vertices) return;

  for (int i = 0; i < num_obj; i++) {
    int id = local_id[i];
    num_edges[i] = graph->neighbor_index[id + 1] - graph->neighbor_index[id];
  }

  *ierr = ZOLTAN_OK;
}

static void get_edge_list(void *data, int size_gid, int size_lid, int num_obj,
                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                          int *num_edges, ZOLTAN_ID_PTR nbr_global_id,
                          int *nbr_procs, int wgt_dim, float *ewgts,
                          int *ierr) {
  graph_t *graph = (graph_t *)data;
  *ierr = ZOLTAN_FATAL;

  if ((size_gid != 1) || (size_lid != 1) || (wgt_dim != 0) ||
      (num_obj != graph->num_vertices))
    return;

  for (int i = 0; i < num_obj; i++) {
    int id = local_id[i];
    if (num_edges[i] !=
        (graph->neighbor_index[id + 1] - graph->neighbor_index[id]))
      return;

    for (int j = graph->neighbor_index[id]; j < graph->neighbor_index[id + 1];
         j++) {
      nbr_global_id[j] = graph->neighbor_ids[j];
      nbr_procs[j] = graph->neighbor_procs[j];
    }
  }

  *ierr = ZOLTAN_OK;
}

#define check_zoltan(status, msg)                                              \
  {                                                                            \
    int rc_ = (rc);                                                            \
    const char *msg_ = msg;                                                    \
    if (rc_ != ZOLTAN_OK) {                                                    \
      fprintf(stderr, msg);                                                    \
      fflush(stderr);                                                          \
      goto err;                                                                \
    }                                                                          \
  }

int Zoltan_partMesh(int *part, long long *vl, int nel, int nv, MPI_Comm comm,
                    int verbose) {
  float ver;
  int rc = Zoltan_Initialize(0, NULL, &ver);

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  check_zoltan(rc, "ERROR: Zoltan_Initialize failed!\n");

  struct Zoltan_Struct *zz = Zoltan_Create(comm);
  // General parameters:
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "2");
  Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.01");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  // Graph parameters:
  Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
  Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG");
  Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");

  // parMETIS options:
  Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "PARMETIS");
  // Zoltan_Set_Param(zz, "", "PARTKWAY");
  // Zoltan_Set_Param(zz,"GRAPH_PACKAGE","SCOTCH");

  // We are going to create the dual graph by hand and input it to
  // Zoltan.
  struct comm c;
  comm_init(&c, comm);
  graph_t *graph = graph_create(vl, nel, nv, &c);

  // Query functions.
  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, graph);
  Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, graph);
  Zoltan_Set_Num_Edges_Multi_Fn(zz, get_edge_size_list, graph);
  Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, graph);

  int changes, num_gid_entries, num_lid_entries, num_import, num_export;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids, export_global_ids;
  ZOLTAN_ID_PTR export_local_ids;
  int *import_procs, *import_to_part, *export_procs, *export_to_part;

  // Now we can partition the graph.
  rc = Zoltan_LB_Partition(zz, &changes, &num_gid_entries, &num_lid_entries,
                           &num_import, &import_global_ids, &import_local_ids,
                           &import_procs, &import_to_part, &num_export,
                           &export_global_ids, &export_local_ids, &export_procs,
                           &export_to_part);

  check_zoltan(rc, "ERROR: Zoltan_LB_Partition failed!\n");

  for (uint i = 0; i < graph->num_vertices; i++) part[i] = rank;
  for (uint i = 0; i < num_export; i++)
    part[export_local_ids[i]] = export_to_part[i];

  Zoltan_LB_Free_Part(&import_global_ids, &import_local_ids, &import_procs,
                      &import_to_part);
  Zoltan_LB_Free_Part(&export_global_ids, &export_local_ids, &export_procs,
                      &export_to_part);
  Zoltan_Destroy(&zz);

  graph_destroy(&graph);

  return 0;
err:
  return 1;
}
#endif

#if defined(ZOLTAN2)

extern int Zoltan2_partMesh(int *part, long long *vl, unsigned nel, int nv,
                            MPI_Comm comm_, int verbose);
#endif // ZOLTAN2

#if defined(KAHIP)
#include <stdbool.h>
#include <kaHIP_interface.h>
#endif // KAHIP

static void print_part_stat(long long *vtx, int nel, int nv, comm_ext ce) {
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
  int nel = *nell;
  int nv = *nve;
  int lelt = *lelm;
  int partitioner = *fpartitioner;
  int algo = *falgo;

  struct comm comm;
#if defined(MPI)
  comm_ext cext = MPI_Comm_f2c(*fcomm);
#else
  comm_ext cext = 0;
#endif
  comm_init(&comm, cext);

  int ierr = 1;
  int *part = (int *)malloc(lelt * sizeof(int));

  if (*loglevel > 2) print_part_stat(vl, nel, nv, cext);

  if (partitioner == 0 || partitioner == 1) {
#if defined(PARRSB)
    parrsb_options options = parrsb_default_options;
    options.partitioner = partitioner;
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
  comm_ext cext = MPI_Comm_f2c(*fcomm);
#else
  comm_ext cext = 0;
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
  comm_ext c = MPI_Comm_f2c(*comm);
#else
  comm_ext c = 0;
#endif

  print_part_stat(vtx, *nel, *nv, c);
}
