#ifndef FINDPT_H
#define FINDPT_H

/* requires "types.h", "poly.h", "tensor.h" */
#if !defined(TYPES_H) || !defined(POLY_H) || !defined(TENSOR_H)
#warning "findpt.h" requires "types.h", "poly.h", "tensor.h"
#endif

typedef struct {
  const real *xw[2];   /* geometry data */
  real *z[2];          /* lobatto nodes */
  lagrange_data ld[2]; /* interpolation, derivative weights & data */
  unsigned nptel;      /* nodes per element */
  struct findpt_hash_data_2 *hash;   /* geometric hashing data */
  struct findpt_listel *list, **sorted, **end;
                                        /* pre-allocated list of elements to
                                           check (found by hashing), and
                                           pre-allocated list of pointers into
                                           the first list for sorting */
  struct findpt_opt_data_2 *od; /* data for the optimization algorithm */
  real *od_work;
} findpt_data_2;

typedef struct {
  const real *xw[3];   /* geometry data */
  real *z[3];          /* lobatto nodes */
  lagrange_data ld[3]; /* interpolation, derivative weights & data */
  unsigned nptel;      /* nodes per element */
  struct findpt_hash_data_3 *hash;   /* geometric hashing data */
  struct findpt_listel *list, **sorted, **end;
                                        /* pre-allocated list of elements to
                                           check (found by hashing), and
                                           pre-allocated list of pointers into
                                           the first list for sorting */
  struct findpt_opt_data_3 *od; /* data for the optimization algorithm */
  real *od_work;
} findpt_data_3;

findpt_data_2 *findpt_setup_2(
          const real *const xw[2], const unsigned n[2], uint nel,
          uint max_hash_size, real bbox_tol);
findpt_data_3 *findpt_setup_3(
          const real *const xw[3], const unsigned n[3], uint nel,
          uint max_hash_size, real bbox_tol);

void findpt_free_2(findpt_data_2 *p);
void findpt_free_3(findpt_data_3 *p);

const real *findpt_allbnd_2(const findpt_data_2 *p);
const real *findpt_allbnd_3(const findpt_data_3 *p);

typedef int (*findpt_func)(void *, const real *, int, uint *, real *, real *);
int findpt_2(findpt_data_2 *p, const real x[2], int guess,
             uint *el, real r[2], real *dist);
int findpt_3(findpt_data_3 *p, const real x[3], int guess,
             uint *el, real r[3], real *dist);

static void findpt_weights_2(findpt_data_2 *p, const real r[2])
{
  lagrange_0(&p->ld[0],r[0]);
  lagrange_0(&p->ld[1],r[1]);
}

static void findpt_weights_3(findpt_data_3 *p, const real r[3])
{
  lagrange_0(&p->ld[0],r[0]);
  lagrange_0(&p->ld[1],r[1]);
  lagrange_0(&p->ld[2],r[2]);
}

static double findpt_eval_2(findpt_data_2 *p, const real *u)
{
  return tensor_i2(p->ld[0].J,p->ld[0].n,
                   p->ld[1].J,p->ld[1].n,
                   u, p->od_work);
}

static double findpt_eval_3(findpt_data_3 *p, const real *u)
{
  return tensor_i3(p->ld[0].J,p->ld[0].n,
                   p->ld[1].J,p->ld[1].n,
                   p->ld[2].J,p->ld[2].n,
                   u, p->od_work);
}

#endif

