#include "zoltan.h"
#include "partitioner.h"

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
