#ifndef POLY_H
#define POLY_H

/* requires "types.h" */
#ifndef TYPES_H
#warning "poly.h" requires "types.h"
#endif

/* 
  For brevity's sake, some names have been shortened
  Quadrature rules
    Gauss   -> Gauss-Legendre quadrature (open)
    Lobatto -> Gauss-Lobatto-Legendre quadrature (closed at both ends)
  Polynomial bases
    Legendre -> Legendre basis
    Gauss    -> Lagrangian basis using Gauss   quadrature nodes
    Lobatto  -> Lagrangian basis using Lobatto quadrature nodes
*/

/*--------------------------------------------------------------------------
   Legendre Polynomial Matrix Computation
   (compute P_i(x_j) for i = 0, ..., n and a given set of x)
  --------------------------------------------------------------------------*/

/* precondition: n >= 1
   inner index is x index (0 ... m-1);
   outer index is Legendre polynomial number (0 ... n)
 */
void legendre_matrix(const real *x, int m, real *P, int n);

/* precondition: n >= 1
   inner index is Legendre polynomial number (0 ... n)
   outer index is x index (0 ... m-1);
 */
void legendre_matrix_t(const real *x, int m, real *P, int n);

/* precondition: n >= 1
   compute P_i(x) with i = 0 ... n
 */
void legendre_row(real x, real *P, int n);


/*--------------------------------------------------------------------------
   Quadrature Nodes and Weights Calculation
   
   call the _nodes function before calling the _weights function
  --------------------------------------------------------------------------*/

void gauss_nodes(real *z, int n);   /* n nodes (order = 2n-1) */
void lobatto_nodes(real *z, int n); /* n nodes (order = 2n-3) */

void gauss_weights(const real *z, real *w, int n);
void lobatto_weights(const real *z, real *w, int n);

/*--------------------------------------------------------------------------
   Lagrangian to Legendre Change-of-basis Matrix
  --------------------------------------------------------------------------*/

/* precondition: n >= 2
   given the Gauss quadrature rule (z,w,n), compute the square matrix J
   for transforming from the Gauss basis to the Legendre basis:
   
      u_legendre(i) = sum_j J(i,j) u_gauss(j)

   computes J   = .5 (2i+1) w  P (z )
             ij              j  i  j
             
   in column major format (inner index is i, the Legendre index)
 */
void gauss_to_legendre(const real *z, const real *w, int n, real *J);

/* precondition: n >= 2
   same as above, but
   in row major format (inner index is j, the Gauss index)
 */
void gauss_to_legendre_t(const real *z, const real *w, int n, real *J);

/* precondition: n >= 3
   given the Lobatto quadrature rule (z,w,n), compute the square matrix J
   for transforming from the Lobatto basis to the Legendre basis:
   
      u_legendre(i) = sum_j J(i,j) u_lobatto(j)

   in column major format (inner index is i, the Legendre index)
 */
void lobatto_to_legendre(const real *z, const real *w, int n, real *J);

/*--------------------------------------------------------------------------
   Lagrangian basis function evaluation
  --------------------------------------------------------------------------*/

/* given the Lagrangian nodes (z,n) and evaluation points (x,m)
   evaluate all Lagrangian basis functions at all points x
   
   inner index of output J is the basis function index (row-major format)
   provide work array with space for 4*n doubles
 */
void lagrange_weights(const real *z, unsigned n,
                      const real *x, unsigned m,
                      real *J, real *work);

/* given the Lagrangian nodes (z,n) and evaluation points (x,m)
   evaluate all Lagrangian basis functions and their derivatives
   
   inner index of outputs J,D is the basis function index (row-major format)
   provide work array with space for 6*n doubles
 */
void lagrange_weights_deriv(const real *z, unsigned n,
                            const real *x, unsigned m,
                            real *J, real *D, real *work);

/*--------------------------------------------------------------------------
   Speedy Lagrangian Interpolation
   
   Usage:
   
     lagrange_data p;
     lagrange_setup(&p,z,n);    // setup for nodes z[0 ... n-1]
     
     the weights
       p->J [0 ... n-1]     interpolation weights
       p->D [0 ... n-1]     1st derivative weights
       p->D2[0 ... n-1]     2nd derivative weights
     are computed for a given x with:
       lagrange_0(p,x);  // compute p->J
       lagrange_1(p,x);  // compute p->J, p->D
       lagrange_2(p,x);  // compute p->J, p->D, p->D2
       lagrange_2u(p);   // compute p->D2 after call of lagrange_1(p,x);
     These functions use the z array supplied to setup
       (that pointer should not be freed between calls)
     Weights for x=z[0] and x=z[n-1] are computed during setup; access as:
       p->J_z0, etc. and p->J_zn, etc.

     lagrange_free(&p);  // deallocate memory allocated by setup
  --------------------------------------------------------------------------*/

typedef struct {
  unsigned n;                /* number of Lagrange nodes            */
  const real *z;             /* Lagrange nodes (user-supplied)      */
  real *J, *D, *D2;          /* weights for 0th,1st,2nd derivatives */
  real *J_z0, *D_z0, *D2_z0; /* ditto at z[0]   (computed at setup) */
  real *J_zn, *D_zn, *D2_zn; /* ditto at z[n-1] (computed at setup) */
  real *w, *d, *u0, *v0, *u1, *v1, *u2, *v2; /* work data           */
} lagrange_data;

void lagrange_setup(lagrange_data *p, const real *z, unsigned n);
void lagrange_free(lagrange_data *p);

static void lagrange_0(lagrange_data *p, real x)
{
  unsigned i, n=p->n;
  for(i=0  ; i<n  ; ++i) p->d[i] = x-p->z[i];
  for(i=0  ; i<n-1; ++i) p->u0[i+1] = p->d[i]*p->u0[i];
  for(i=n-1; i    ; --i) p->v0[i-1] = p->d[i]*p->v0[i];
  for(i=0  ; i<n  ; ++i) p->J[i] = p->w[i]*p->u0[i]*p->v0[i];
}

static void lagrange_1(lagrange_data *p, real x)
{
  unsigned i, n=p->n;
  for(i=0  ; i<n  ; ++i) p->d[i] = x-p->z[i];
  for(i=0  ; i<n-1; ++i)
    p->u0[i+1] = p->d[i]*p->u0[i],
    p->u1[i+1] = p->d[i]*p->u1[i] + p->u0[i];
  for(i=n-1; i    ; --i)
    p->v0[i-1] = p->d[i]*p->v0[i],
    p->v1[i-1] = p->d[i]*p->v1[i] + p->v0[i];
  for(i=0  ; i<n  ; ++i)
    p->J[i] = p->w[i]*p->u0[i]*p->v0[i],
    p->D[i] = p->w[i]*(p->u1[i]*p->v0[i]+p->u0[i]*p->v1[i]);
}

static void lagrange_2(lagrange_data *p, real x)
{
  unsigned i,n=p->n;
  for(i=0  ; i<n  ; ++i) p->d[i] = x-p->z[i];
  for(i=0  ; i<n-1; ++i)
    p->u0[i+1]=p->d[i]*p->u0[i],
    p->u1[i+1]=p->d[i]*p->u1[i]+p->u0[i],
    p->u2[i+1]=p->d[i]*p->u2[i]+2*p->u1[i];
  for(i=n-1; i    ; --i)
    p->v0[i-1]=p->d[i]*p->v0[i],
    p->v1[i-1]=p->d[i]*p->v1[i]+p->v0[i],
    p->v2[i-1]=p->d[i]*p->v2[i]+2*p->v1[i];
  for(i=0  ; i<n  ; ++i)
    p->J [i]=p->w[i]*p->u0[i]*p->v0[i],
    p->D [i]=p->w[i]*(p->u1[i]*p->v0[i]+p->u0[i]*p->v1[i]),
    p->D2[i]=p->w[i]*(p->u2[i]*p->v0[i]+2*p->u1[i]*p->v1[i]+p->u0[i]*p->v2[i]);
}

static void lagrange_2u(lagrange_data *p)
{
  unsigned i,n=p->n;
  for(i=0  ; i<n-1; ++i)
    p->u2[i+1]=p->d[i]*p->u2[i]+2*p->u1[i];
  for(i=n-1; i    ; --i)
    p->v2[i-1]=p->d[i]*p->v2[i]+2*p->v1[i];
  for(i=0  ; i<n  ; ++i)
    p->D2[i]=p->w[i]*(p->u2[i]*p->v0[i]+2*p->u1[i]*p->v1[i]+p->u0[i]*p->v2[i]);
}

#endif

