#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>    /* for cos, fabs */
#include <string.h>  /* for memcpy */
#include <float.h>

#include "errmem.h"
#include "types.h"

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
void legendre_matrix(const real *x, int m, real *P, int n)
{
  int i,j;
  real *Pjm1 = P, *Pj = Pjm1+m, *Pjp1 = Pj+m;
  for(i=0; i<m; ++i) Pjm1[i] = 1;
  for(i=0; i<m; ++i) Pj[i] = x[i];
  for(j=1; j<n; ++j) {
    real c = 1/(real)(j+1), a = c*(2*j+1), b = c*j;
    for(i=0; i<m; ++i) Pjp1[i] = a*x[i]*Pj[i]-b*Pjm1[i];
    Pjp1+=m, Pj+=m, Pjm1+=m;
  }
}

/* precondition: n >= 1, n even */
static void legendre_row_even(real x, real *P, int n)
{
  real p[2] = {1, x};
  int i;
  P[0] = 1, P[1] = x;
  for(i=1; i<=n-2; i+=2) {
    p[0] = ((2*i+1)*x*p[1]- i   *p[0])/(i+1);
    p[1] = ((2*i+3)*x*p[0]-(i+1)*p[1])/(i+2);
    P[i+1] = p[0];
    P[i+2] = p[1];
  }
  P[n] = ((2*n-1)*x*p[1]-(n-1)*p[0])/n;
}

/* precondition: n >= 1, n odd */
static void legendre_row_odd(real x, real *P, int n)
{
  real p[2] = {1, x};
  int i;
  P[0] = 1, P[1] = x;
  for(i=1; i<=n-2; i+=2) {
    p[0] = ((2*i+1)*x*p[1]- i   *p[0])/(i+1);
    p[1] = ((2*i+3)*x*p[0]-(i+1)*p[1])/(i+2);
    P[i+1] = p[0];
    P[i+2] = p[1];
  }
}

/* precondition: n >= 1
   compute P_i(x) with i = 0 ... n
 */
void legendre_row(real x, real *P, int n)
{
  if(n&1) legendre_row_odd(x,P,n); else legendre_row_even(x,P,n);
}

/* precondition: n >= 1
   inner index is Legendre polynomial number (0 ... n)
   outer index is x index (0 ... m-1);
 */
void legendre_matrix_t(const real *x, int m, real *P, int n)
{
  int i;
  if(n&1) for(i=0;i<m;++i,P+=n+1) legendre_row_odd(x[i],P,n);
     else for(i=0;i<m;++i,P+=n+1) legendre_row_even(x[i],P,n);
}

/*--------------------------------------------------------------------------
   Legendre Polynomial Computation
   compute P_n(x) or P_n'(x) or P_n''(x)
  --------------------------------------------------------------------------*/

/* precondition: n >= 0 */
static real legendre(int n, real x)
{
  real p[2] = {1, x};
  int i;
  for(i=1; i<n; i+=2) {
    p[0] = ((2*i+1)*x*p[1]- i   *p[0])/(i+1);
    p[1] = ((2*i+3)*x*p[0]-(i+1)*p[1])/(i+2);
  }
  return p[n&1];
}

/* precondition: n > 0 */
static real legendre_d1(int n, real x)
{
  real p[2] = {3*x, 1};
  int i;
  for(i=2; i<n; i+=2) {
    p[1] = ((2*i+1)*x*p[0]-(i+1)*p[1])/i;
    p[0] = ((2*i+3)*x*p[1]-(i+2)*p[0])/(i+1);
  }
  return p[n&1];
}

/* precondition: n > 1 */
static real legendre_d2(int n, real x)
{
  real p[2] = {3, 15*x};
  int i;
  for(i=3; i<n; i+=2) {
    p[0] = ((2*i+1)*x*p[1]-(i+2)*p[0])/(i-1);
    p[1] = ((2*i+3)*x*p[0]-(i+3)*p[1])/i;
  }
  return p[n&1];
}

/*--------------------------------------------------------------------------
   Quadrature Nodes and Weights Calculation
   compute the n Gauss-Legendre nodes and weights or
           the n Gauss-Lobatto-Legendre nodes and weights
  --------------------------------------------------------------------------*/

/* n nodes */
void gauss_nodes(real *z, int n)
{
  int i,j;
  for(i=0; i<=n/2-1; ++i) {
    real ox, x = cosr( (2*n-2*i-1)*(PI/2)/n );
    do {
      ox = x;
      x -= legendre(n,x)/legendre_d1(n,x);
    } while(fabsr(x-ox)>-x*EPS);
    z[i] = x - legendre(n,x)/legendre_d1(n,x);
  }
  if(n&1) z[n/2]=0;
  for(j=(n+1)/2,i=n/2-1; j<n; ++j,--i) z[j]=-z[i];
}

/* n inner lobatto nodes (excluding -1,1) */
static void lobatto_nodes_aux(real *z, int n)
{
  int i,j,np=n+1;
  for(i=0; i<=n/2-1; ++i) {
    real ox, x = cosr( (n-i)*PI/np );
    do {
      ox = x;
      x -= legendre_d1(np,x)/legendre_d2(np,x);
    } while(fabsr(x-ox)>-x*EPS);
    z[i] = x - legendre_d1(np,x)/legendre_d2(np,x);
  }
  if(n&1) z[n/2]=0;
  for(j=(n+1)/2,i=n/2-1; j<n; ++j,--i) z[j]=-z[i];
}

/* n lobatto nodes */
void lobatto_nodes(real *z, int n)
{
  z[0] = -1, z[n-1] = 1;
  lobatto_nodes_aux(&z[1],n-2);
}

void gauss_weights(const real *z, real *w, int n)
{
  int i,j;
  for(i=0; i<=(n-1)/2; ++i) {
    real d = (n+1)*legendre(n+1,z[i]);
    w[i] = 2*(1-z[i]*z[i])/(d*d);
  }
  for(j=(n+1)/2,i=n/2-1; j<n; ++j,--i) w[j]=w[i];
}

void lobatto_weights(const real *z, real *w, int n)
{
  int i,j;
  for(i=0; i<=(n-1)/2; ++i) {
    real d = legendre(n-1,z[i]);
    w[i] = 2/((n-1)*n*d*d);
  }
  for(j=(n+1)/2,i=n/2-1; j<n; ++j,--i) w[j]=w[i];
}

/*--------------------------------------------------------------------------
   Lagrangian to Legendre Change-of-basis Matrix
   where the nodes of the Lagrangian basis are the GL or GLL quadrature nodes
  --------------------------------------------------------------------------*/

/* precondition: n >= 2
   given the Gauss quadrature rule (z,w,n), compute the square matrix J
   for transforming from the Gauss basis to the Legendre basis:
   
      u_legendre(i) = sum_j J(i,j) u_gauss(j)

   computes J   = .5 (2i+1) w  P (z )
             ij              j  i  j
             
   in column major format (inner index is i, the Legendre index)
 */
void gauss_to_legendre(const real *z, const real *w, int n, real *J)
{
  int i,j;
  legendre_matrix_t(z,n,J,n-1);
  for(j=0;j<n;++j) {
    real ww = w[j];
    for(i=0; i<n; ++i) *J++ *= (2*i+1) * ww/2;
  }
}

/* precondition: n >= 2
   same as above, but
   in row major format (inner index is j, the Gauss index)
 */
void gauss_to_legendre_t(const real *z, const real *w, int n, real *J)
{
  int i,j;
  legendre_matrix(z,n,J,n-1);
  for(i=0;i<n;++i) {
    real ii = (real)(2*i+1)/2;
    for(j=0; j<n; ++j) *J++ *= ii * w[j];
  }
}

/* precondition: n >= 3
   given the Lobatto quadrature rule (z,w,n), compute the square matrix J
   for transforming from the Gauss basis to the Legendre basis:
   
      u_legendre(i) = sum_j J(i,j) u_lobatto(j)

   in column major format (inner index is i, the Legendre index)
 */
void lobatto_to_legendre(const real *z, const real *w, int n, real *J)
{
  int i,j,m=(n+1)/2;
  real *p = J, *q;
  real ww, sum;
  if(n&1)
    for(j=0;j<m;++j,p+=n) legendre_row_odd(z[j],p,n-2);
  else
    for(j=0;j<m;++j,p+=n) legendre_row_even(z[j],p,n-2);
  p = J;
  for(j=0;j<m;++j) {
    ww = w[j], sum = 0;
    for(i=0; i<n-1; ++i) *p *= (2*i+1) * ww/2, sum += *p++;
    *p++ = -sum;
  }
  q = J+(n/2-1)*n;
  if(n&1)
    for(;j<n;++j,p+=n,q-=n) {
      for(i=0; i<n-1; i+=2) p[i]=q[i], p[i+1]=-q[i+1];
      p[i]=q[i];
    }
  else
    for(;j<n;++j,p+=n,q-=n) {
      for(i=0; i<n-1; i+=2) p[i]=q[i], p[i+1]=-q[i+1];
    }
}

/*--------------------------------------------------------------------------
   Lagrangian to Lagrangian change-of-basis matrix, and derivative matrix
  --------------------------------------------------------------------------*/

/* given the Lagrangian nodes (z,n) and evaluation points (x,m)
   evaluate all Lagrangian basis functions at all points x
   
   inner index of output J is the basis function index (row-major format)
   provide work array with space for 4*n reals
 */
void lagrange_weights(const real *z, unsigned n,
                      const real *x, unsigned m,
                      real *J, real *work)
{
  unsigned i,j;
  real *w = work, *d = w+n, *u = d+n, *v = u+n;
  for(i=0; i<n; ++i) {
    real ww = 1, zi = z[i];
    for(j=0; j<i; ++j) ww *= zi-z[j];
    for(++j; j<n; ++j) ww *= zi-z[j];
    w[i] = 1/ww;
  }
  u[0] = v[n-1] = 1;
  for(i=0; i<m; ++i) {
    real xi = x[i];
    for(j=0; j<n; ++j) d[j] = xi - z[j];
    for(j=0; j<n-1; ++j) u[j+1] = d[j] * u[j];
    for(j=n-1; j; --j) v[j-1] = d[j] * v[j];
    for(j=0; j<n; ++j) *J++ = w[j] * u[j] * v[j];
  }
}

/* given the Lagrangian nodes (z,n) and evaluation points (x,m)
   evaluate all Lagrangian basis functions and their derivatives
   
   inner index of outputs J,D is the basis function index (row-major format)
   provide work array with space for 6*n reals
 */
void lagrange_weights_deriv(const real *z, unsigned n,
                            const real *x, unsigned m,
                            real *J, real *D, real *work)
{
  unsigned i,j;
  real *w = work, *d = w+n, *u = d+n, *v = u+n, *up= v+n, *vp=up+n;
  for(i=0; i<n; ++i) {
    real ww = 1, zi = z[i];
    for(j=0; j<i; ++j) ww *= zi-z[j];
    for(++j; j<n; ++j) ww *= zi-z[j];
    w[i] = 1/ww;
  }
  u [0] = v [n-1] = 1;
  up[0] = vp[n-1] = 0;
  for(i=0; i<m; ++i) {
    real xi = x[i];
    for(j=0; j<n; ++j) d[j] = xi - z[j];
    for(j=0; j<n-1; ++j) u[j+1]=d[j]*u[j], up[j+1]=d[j]*up[j]+u[j];
    for(j=n-1; j; --j)   v[j-1]=d[j]*v[j], vp[j-1]=d[j]*vp[j]+v[j];
    for(j=0; j<n; ++j) *J++ = w[j]*u[j]*v[j], 
                       *D++ = w[j]*(up[j]*v[j]+u[j]*vp[j]);
  }
}

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
  real *w, *d, *u0, *v0, *u1, *v1, *u2, *v2; /* work data            */
} lagrange_data;

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

void lagrange_setup(lagrange_data *p, const real *z, unsigned n)
{
  unsigned i,j;
  p->n=n, p->z=z;
  p->w = tmalloc(real, 17*n);
  p->d = p->w+n;
  p->J = p->d+n, p->D = p->J+n, p->D2 = p->D+n;
  p->u0=p->D2+n, p->v0=p->u0+n;
  p->u1=p->v0+n, p->v1=p->u1+n;
  p->u2=p->v1+n, p->v2=p->u2+n;
  p->J_z0=p->v2+n, p->D_z0=p->J_z0+n, p->D2_z0=p->D_z0+n;
  p->J_zn=p->D2_z0+n, p->D_zn=p->J_zn+n, p->D2_zn=p->D_zn+n;
  for(i=0; i<n; ++i) {
    real ww = 1, zi = z[i];
    for(j=0; j<i; ++j) ww *= zi-z[j];
    for(++j; j<n; ++j) ww *= zi-z[j];
    p->w[i] = 1/ww;
  }
  p->u0[0] = p->v0[n-1] = 1;
  p->u1[0] = p->v1[n-1] = 0;
  p->u2[0] = p->v2[n-1] = 0;
  lagrange_2(p,z[0  ]); memcpy(p->J_z0,p->J,3*n*sizeof(real));
  lagrange_2(p,z[n-1]); memcpy(p->J_zn,p->J,3*n*sizeof(real));
}

void lagrange_free(lagrange_data *p)
{
  free(p->w);
}

