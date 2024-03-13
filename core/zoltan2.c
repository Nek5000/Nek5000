#include <gslib.h>

struct laplacian_private {
  uint size;
  ulong *rows, *columns;
  double *values;
};
typedef struct laplacian_private *laplacian_t;

static void free_(void **ptr) { free(*ptr), *ptr = NULL; }
#define tfree(ptr) free_((void **)(ptr))

laplacian_t laplacian_weighted(long long *vl, unsigned nel, unsigned nv,
                               MPI_Comm comm, int verbose) {
  struct comm c;
  comm_init(&c, comm);

  slong element_offset = 0;
  {
    slong out[2], wrk[2], in = nel;
    comm_scan(out, &c, gs_long, gs_add, &in, 1, &wrk);
    element_offset = out[0];
  }

  typedef struct {
    ulong e, v;
    uint p;
  } pair_t;

  struct array arr;
  {
    array_init(pair_t, &arr, nel * nv + 1);

    pair_t p;
    for (uint i = 0; i < nel; i++) {
      p.e = element_offset + i;
      for (uint j = 0; j < nv; j++) {
        p.v = vl[i * nv + j];
        p.p = p.v % c.np;
        array_cat(pair_t, &arr, &p, 1);
      }
    }
  }

  typedef struct {
    ulong e1, e2;
    uint p;
    int w;
  } neighbor_t;

  struct array neighbors;
  array_init(neighbor_t, &neighbors, 27 * nel + 1);

  struct crystal cr;
  {
    crystal_init(&cr, &c);
    sarray_transfer(pair_t, &arr, p, 1, &cr);

    buffer bfr;
    buffer_init(&bfr, arr.n + 1);

    sarray_sort(pair_t, arr.ptr, arr.n, v, 1, &bfr);

    pair_t *pa = (pair_t *)arr.ptr;
    uint s = 0;
    while (s < arr.n) {
      uint e = s + 1;
      for (; e < arr.n && pa[e].v == pa[s].v; e++)
        ;
      for (uint i = s; i < e; i++) {
        neighbor_t nbr = {.e1 = pa[i].e, .p = pa[i].p};
        for (uint j = s; j < e; j++) {
          nbr.e2 = pa[j].e;
          array_cat(neighbor_t, &neighbors, &nbr, 1);
        }
      }
      s = e;
    }

    sarray_transfer(neighbor_t, &neighbors, p, 0, &cr);
    sarray_sort_2(neighbor_t, neighbors.ptr, neighbors.n, e1, 1, e2, 1, &bfr);
    buffer_free(&bfr);
  }
  array_free(&arr);
  crystal_free(&cr);
  comm_free(&c);

  {
    neighbor_t *pn = (neighbor_t *)neighbors.ptr;
    uint s = 0, num_unique = 0;
    while (s < neighbors.n) {
      uint e = s + 1;
      for (; e < neighbors.n && pn[e].e1 == pn[s].e1 && pn[e].e2 == pn[s].e2;
           e++)
        ;
      pn[num_unique] = pn[s];
      pn[num_unique++].w = -(int)(e - s);
      s = e;
    }
    neighbors.n = num_unique;
  }

  laplacian_t L = tcalloc(struct laplacian_private, 1);
  {
    neighbor_t *pn = (neighbor_t *)neighbors.ptr;
    uint s = 0;
    while (s < neighbors.n) {
      int sum = (pn[s].e1 == pn[s].e2) ? 0 : pn[s].w;
      int diag = (pn[s].e1 == pn[s].e2) ? s : -1;
      uint e = s + 1;
      for (; e < neighbors.n && pn[e].e1 == pn[s].e1; e++) {
        if (pn[e].e1 != pn[e].e2) {
          sum += pn[e].w;
          continue;
        }
        diag = e;
      }
      assert(diag >= 0 && sum < 0);
      pn[diag].w = -sum;
      s = e;
    }

    L->size = neighbors.n;
    L->rows = tcalloc(ulong, neighbors.n);
    L->columns = tcalloc(ulong, neighbors.n);
    L->values = tcalloc(double, neighbors.n);
    for (uint i = 0; i < L->size; i++) {
      L->rows[i] = pn[i].e1;
      L->columns[i] = pn[i].e2;
      L->values[i] = pn[i].w;
    }
  }

  array_free(&neighbors);

  return L;
}

uint laplacian_size(laplacian_t L) { return L->size; }

ulong *laplacian_rows(laplacian_t L) { return L->rows; }

ulong *laplacian_columns(laplacian_t L) { return L->columns; }

double *laplacian_values(laplacian_t L) { return L->values; }

void laplacian_free(laplacian_t *L_) {
  laplacian_t L = *L_;
  if (!L) return;
  tfree(&L->rows), tfree(&L->columns), tfree(&L->values);
  tfree(L_);
}

#undef tfree
