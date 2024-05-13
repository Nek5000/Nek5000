#include "zoltan.h"
#include "partitioner.h"

typedef struct {
  uint   num_vertices;
  ulong *vertex_ids;
  int   *neighbor_index;
  ulong *neighbor_ids;
  int   *neighbor_procs;
  int   *neighbor_weights;
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
    uint  dest, np;
  } vertex_t;

  struct array vertices;
  {
    array_init(vertex_t, &vertices, nelt * nv);
    vertex_t vertex;
    for (uint i = 0; i < nelt; i++) {
      vertex.eid = i + start + 1;
      for (uint j = 0; j < nv; j++) {
        vertex.vid  = vl[i * nv + j];
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
    uint                  vn = vertices.n;
    uint                  s  = 0;
    while (s < vn) {
      uint e = s + 1;
      while (e < vn && pv[s].vid == pv[e].vid) e++;

      for (uint j = s; j < e; j++) {
        vertex_t v = pv[j];
        for (uint k = s; k < e; k++) {
          v.vid = pv[k].eid;
          v.np  = pv[k].dest;
          array_cat(vertex_t, &neighbors, &v, 1);
        }
      }
      s = e;
    }
  }
  array_free(&vertices);

  sarray_transfer(vertex_t, &neighbors, dest, 0, &cr);
  sarray_sort_2(vertex_t, neighbors.ptr, neighbors.n, eid, 1, vid, 1, &bfr);
  crystal_free(&cr);
  buffer_free(&bfr);

  // Compre the neighbors array to find the number of neighbors for each
  // element.
  typedef struct {
    ulong nid, eid;
    uint  np;
    int   weight;
  } neighbor_t;

  struct array compressed;
  {
    array_init(neighbor_t, &compressed, neighbors.n);

    const vertex_t *const pn = (const vertex_t *const)neighbors.ptr;
    uint                  nn = neighbors.n;
    neighbor_t            nbr;
    uint                  s = 0;
    while (s < nn) {
      uint e = s + 1;
      while (e < nn && pn[s].eid == pn[e].eid && pn[s].vid == pn[e].vid) e++;

      // Sanity check.
      for (uint j = s + 1; j < e; j++) assert(pn[s].np == pn[j].np);

      // Add to compressed array.
      nbr.eid    = pn[s].eid;
      nbr.nid    = pn[s].vid;
      nbr.np     = pn[s].np;
      nbr.weight = (int)(e - s);
      if (nbr.eid != nbr.nid) array_cat(neighbor_t, &compressed, &nbr, 1);
      s = e;
    }
  }
  array_free(&neighbors);

  graph_t *graph      = tcalloc(graph_t, 1);
  graph->num_vertices = nelt;
  graph->vertex_ids   = tcalloc(ulong, nelt);
  for (uint i = 0; i < nelt; i++) graph->vertex_ids[i] = i + start + 1;

  graph->neighbor_index = tcalloc(int, nelt + 1);
  {
    neighbor_t *pc    = (neighbor_t *)compressed.ptr;
    uint        count = 0;
    uint        s     = 0;
    while (s < compressed.n) {
      uint e = s + 1;
      while (e < compressed.n && pc[s].eid == pc[e].eid) e++;
      count++;
      graph->neighbor_index[count] = e;
      s                            = e;
    }
    assert(count == nelt);
    assert(graph->neighbor_index[count] == compressed.n);
  }

  graph->neighbor_ids     = tcalloc(ulong, compressed.n);
  graph->neighbor_procs   = tcalloc(int, compressed.n);
  graph->neighbor_weights = tcalloc(int, compressed.n);
  {
    neighbor_t *pc = (neighbor_t *)compressed.ptr;
    for (uint i = 0; i < compressed.n; i++) {
      graph->neighbor_ids[i]     = pc[i].nid;
      graph->neighbor_procs[i]   = pc[i].np;
      graph->neighbor_weights[i] = pc[i].weight;
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
  free(graph->neighbor_weights);
  free(graph);
  graph_ = NULL;
}

static int get_number_of_vertices(void *data, int *ierr) {
  graph_t *graph = (graph_t *)data;
  *ierr          = ZOLTAN_OK;
  return graph->num_vertices;
}

static void get_vertex_list(void *data, int size_gid, int size_lid,
                            ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                            int wgt_dim, float *obj_wgts, int *ierr) {
  graph_t *graph = (graph_t *)data;
  *ierr          = ZOLTAN_FATAL;

  for (int i = 0; i < graph->num_vertices; i++) {
    global_id[i] = graph->vertex_ids[i];
    local_id[i]  = i;
  }

  *ierr = ZOLTAN_OK;
}

static void get_edge_size_list(void *data, int size_gid, int size_lid,
                               int num_obj, ZOLTAN_ID_PTR global_id,
                               ZOLTAN_ID_PTR local_id, int *num_edges,
                               int *ierr) {
  graph_t *graph = (graph_t *)data;
  *ierr          = ZOLTAN_FATAL;

  if (size_gid != 1 || size_lid != 1 || num_obj != graph->num_vertices) return;

  for (int i = 0; i < num_obj; i++) {
    int id       = local_id[i];
    num_edges[i] = graph->neighbor_index[id + 1] - graph->neighbor_index[id];
  }

  *ierr = ZOLTAN_OK;
}

static void get_edge_list(void *data, int size_gid, int size_lid, int num_obj,
                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                          int *num_edges, ZOLTAN_ID_PTR nbr_global_id,
                          int *nbr_procs, int wgt_dim, float *neighbor_weights,
                          int *ierr) {
  graph_t *graph = (graph_t *)data;
  *ierr          = ZOLTAN_FATAL;

  if ((size_gid != 1) || (size_lid != 1) || (wgt_dim != 1) ||
      (num_obj != graph->num_vertices))
    return;

  for (int i = 0; i < num_obj; i++) {
    int id = local_id[i];
    if (num_edges[i] !=
        (graph->neighbor_index[id + 1] - graph->neighbor_index[id]))
      return;

    for (int j = graph->neighbor_index[id]; j < graph->neighbor_index[id + 1];
         j++) {
      nbr_global_id[j]    = graph->neighbor_ids[j];
      nbr_procs[j]        = graph->neighbor_procs[j];
      neighbor_weights[j] = graph->neighbor_weights[j];
    }
  }

  *ierr = ZOLTAN_OK;
}

#define check_zoltan(status, msg)                                              \
  {                                                                            \
    int         rc_  = (rc);                                                   \
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
  int   rc = Zoltan_Initialize(0, NULL, &ver);

  int    rank, size;
  double imbalance_tol;
  {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    long long nelg = 0, nel_ = nel;
    MPI_Allreduce(&nel_, &nelg, 1, MPI_LONG_LONG, MPI_SUM, comm);
    nel_          = nelg / size;
    imbalance_tol = (nel_ + 1.0) / nel_;

    if (rank == 0 && verbose) {
      fprintf(stderr, "\nnum_global_elements = %lld imbalance_tol = %lf\n",
              nelg, imbalance_tol);
      fflush(stderr);
    }
  }

  check_zoltan(rc, "ERROR: Zoltan_Initialize failed!\n");

  struct Zoltan_Struct *zz = Zoltan_Create(comm);
  // General parameters:
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "2");
  Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");

  char imbalance_tol_str[32] = {0};
  snprintf(imbalance_tol_str, 32, "%lf", imbalance_tol);
  Zoltan_Set_Param(zz, "IMBALANCE_TOL", imbalance_tol_str);

  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
  Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");

  // Graph parameters:
  Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");

  Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG");
  Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");

  // parMETIS options:
  const char *enable_parmetis = getenv("ENABLE_ZOLTAN_PARMETIS");
  if (!enable_parmetis || atoi(enable_parmetis) == 0) goto enable_scotch;
  Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "PARMETIS");
  Zoltan_Set_Param(zz, "PARMETIS_METHOD", "PARTKWAY");

  // Scotch options:
  const char *enable_scotch = NULL;
enable_scotch:
  enable_scotch = getenv("ENABLE_ZOLTAN_SCOTCH");
  if (!enable_scotch || atoi(enable_scotch) == 0) goto create_dual_graph;
  Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "SCOTCH");

  // We are going to create the dual graph by hand and input it to Zoltan.
  struct comm c;
create_dual_graph:
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
  int          *import_procs, *import_to_part, *export_procs, *export_to_part;

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

#undef check_zoltan
