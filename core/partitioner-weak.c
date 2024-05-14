#include "partitioner.h"

graph_t *graph_create(long long *vl, int nelt, int nv, struct comm *c) {
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

void graph_print(graph_t *graph) {
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

void graph_destroy(graph_t **graph_) {
  graph_t *graph = *graph_;
  free(graph->vertex_ids);
  free(graph->neighbor_index);
  free(graph->neighbor_ids);
  free(graph->neighbor_procs);
  free(graph->neighbor_weights);
  free(graph);
  graph_ = NULL;
}

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

int Zoltan_partMesh(int *part, long long *vl, int nel, int nv, double *opt,
                    MPI_Comm comm) __attribute__((weak));

int Zoltan_partMesh(int *part, long long *vl, int nel, int nv, double *opt,
                    MPI_Comm comm) {
  exit_if_not_enabled("Zoltan", "ZOLTAN", comm);
}

int Zoltan2_partMesh(int *part, long long *vl, unsigned nel, int nv,
                     double *opt, MPI_Comm comm) __attribute__((weak));

int Zoltan2_partMesh(int *part, long long *vl, unsigned nel, int nv,
                     double *opt, MPI_Comm comm) {
  exit_if_not_enabled("Zoltan2", "ZOLTAN2", comm);
}

int parHIP_partMesh(int *part, long long *vl, int nel, int nv, double *opt,
                    MPI_Comm comm) __attribute__((weak));

int parHIP_partMesh(int *part, long long *vl, int nel, int nv, double *opt,
                    MPI_Comm comm) {
  exit_if_not_enabled("parHIP", "PARHIP", comm);
}
