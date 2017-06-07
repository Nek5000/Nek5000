#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>    /* for cos, fabs */
#include <float.h>
#include <string.h>  /* for memcpy */

#include "errmem.h"
#include "types.h"
#include "minmax.h"
#include "poly.h"
#include "tensor.h"

/*--------------------------------------------------------------------------
   Lobatto Polynomial Bounds

   Needed inputs are the Gauss-Lobatto quadrature nodes and weights:   
     unsigned nr = ..., ns = ...;
     real zr[nr], wr[nr];
     real zs[ns], ws[ns];
     
     lobatto_nodes(zr,nr); lobatto_weights(zr,wr,nr);
     lobatto_nodes(zs,ns); lobatto_weights(zs,ws,ns);

   The number of points in the constructed piecewise (bi-)linear bounds
   is a parameter; more points give tighter bounds
   
     unsigned mr = 2*nr, ms = 2*ns;
   
   The necessary setup is accomplished via:
     lob_bnd_base b_data_r;
     lob_bnd_ext  e_data_s;
     
     lob_bnd_base_alloc(&b_data_r,nr,mr);
     lob_bnd_base_setup(&b_data_r,zr,wr);
     lob_bnd_ext_alloc(&e_data_s,ns,ms);
     lob_bnd_ext_setup(&e_data_s,zs,ws);
   
   Bounds may then be computed via:
     real work1r[2*mr], work1s[2*ms], work2[2*mr + 2*mr*ns + 2*mr*ms];
     real ur[nr], us[ns];    // 1-d polynomials on the zr[] and zs[] nodes
     real u[ns][nr];         // 2-d polynomial on zr[] (x) zs[]
     real bound[2];          // = { min, max } (to be computed)
     
     lob_bnd_1(&b_data_r  ,ur,bound,work1r); // compute bounds on ur
     lob_bnd_1(&e_data_s.b,us,bound,work1s); // compute bounds on us
     lob_bnd_2(&b_data_r, &e_data_s,
               (const double*)&u[0][0],bound,work2); // compute bounds on u
   The above routines access the zr,zs arrays passed to *_setup
     (so do not delete them between calls)
   
   Memory allocated in *_setup is freed with
     lob_bnd_base_free(&b_data_r);
     lob_bnd_ext_free(&e_data_s);

  --------------------------------------------------------------------------*/
  
typedef struct {
  unsigned n; /* number of Lobatto nodes in input */
  unsigned m; /* number of Chebyshev nodes used to calculate bounds */
  real *Q0, *Q1; /* Q0[n], Q1[n] -- first two rows of change of basis matrix
                    from Lobatto node Lagrangian to Legendre */
  const real *z; /* z[n] -- external; Lobatto nodes */
  real *h;       /* h[m] -- Chebyshev nodes */
  real *uv, *ov; /* uv[n][m], ov[n][m] --
                      uv[j][:] is a piecewise linear function in the nodal
                               basis with nodes h[m] that is everywhere less
                               than or equal to the jth Lagrangian basis
                               function (with the Lobatto nodes z[n])
                      ov[j][:] is everywhere greater than or equal */
} lob_bnd_base;

typedef struct {
  lob_bnd_base b;
  real *uvp, *uvn, *ovp, *ovn; /* [n][m] -- uv and ov split into
                                    positive and negative parts */
} lob_bnd_ext;

static void lob_bnd_base_alloc(lob_bnd_base *p, unsigned n, unsigned m)
{
  p->n = n, p->m = m;
  p->Q0 = tmalloc(real,2*n+m+2*n*m);
  p->Q1 = p->Q0+n;
  p->h  = p->Q1+n;
  p->uv = p->h +m;
  p->ov = p->uv+n*m;
}

static void lob_bnd_base_free(lob_bnd_base *p)
{
  free(p->Q0);
}

static void lob_bnd_ext_alloc(lob_bnd_ext *p, unsigned n, unsigned m)
{
  p->b.n = n, p->b.m = m;
  p->b.Q0 = tmalloc(real,2*n+m+6*n*m);
  p->b.Q1 = p->b.Q0+n;
  p->b.h  = p->b.Q1+n;
  p->b.uv = p->b.h +m;
  p->b.ov = p->b.uv+n*m;
  p->uvp  = p->b.ov+n*m;
  p->uvn  = p->uvp +n*m;
  p->ovp  = p->uvn +n*m;
  p->ovn  = p->ovp +n*m;
}

static void lob_bnd_ext_free(lob_bnd_ext *p)
{
  free(p->b.Q0);
}

static void lob_bnd_base_setup(lob_bnd_base *p, const real *z, const real *w)
{
  unsigned i,j,m=p->m,n=p->n,mm=2*m-1;
  real *q = tmalloc(real,(2*n+1)*mm+6*n),
       *J = q+mm, *D = J+n*mm, *work = D+n*mm;
  p->z = z;
  for(i=0;i<n;++i) p->Q0[i]=w[i]/2, p->Q1[i] = 3*p->Q0[i]*z[i];
  p->h[0] = -1, p->h[m-1] = 1;
  for(j=1;j<m-1;++j) p->h[j] = cosr((m-j-1)*PI/(m-1));
  for(j=0;j<m-1;++j) q[2*j] = p->h[j], q[2*j+1] = (p->h[j]+p->h[j+1])/2;
  q[mm-1] = p->h[m-1];
  lagrange_weights_deriv(z,n,q,mm,J,D,work);
  for(i=0;i<n;++i) {
    real *uv = p->uv+i*m, *ov = p->ov+i*m;
    ov[0]   = uv[0]   = J[i];
    ov[m-1] = uv[m-1] = J[(mm-1)*n+i];
    for(j=1;j<m-1;++j) {
      unsigned jj = 2*j;
      real c0 = J[(jj-1)*n+i] + (q[jj]-q[jj-1])*D[(jj-1)*n+i];
      real c1 = J[(jj+1)*n+i] + (q[jj]-q[jj+1])*D[(jj+1)*n+i];
      real c2 = J[jj*n+i];
      rminmax_3(&uv[j],&ov[j],c0,c1,c2);
    }
  }
  free(q);
}

static void lob_bnd_ext_setup(lob_bnd_ext *p, const real *z, const real *w)
{
  unsigned i, mn = p->b.m * p->b.n;
  lob_bnd_base_setup(&p->b,z,w);
  for(i=0;i<mn;++i) {
    real uvi = p->b.uv[i], ovi = p->b.ov[i];
    p->uvp[i] = p->uvn[i] = p->ovp[i] = p->ovn[i] = 0;
    if(uvi > 0) p->uvp[i]=uvi; else p->uvn[i]=uvi;
    if(ovi > 0) p->ovp[i]=ovi; else p->ovn[i]=ovi;
  }
}

static void lob_bnd_lines(const lob_bnd_base *p, const real *u,
                          real *a, real *b)
{
  unsigned i,j;
  real a0=0, a1=0;
  const real *uv = p->uv, *ov = p->ov;
  for(i=0;i<p->n;++i) a0 += p->Q0[i]*u[i], a1 += p->Q1[i]*u[i];
  for(j=0;j<p->m;++j) b[j] = a[j] = a0 + a1*p->h[j];
  for(i=0;i<p->n;++i) {
    real w = u[i] - (a0 + a1*p->z[i]);
    if(w>=0) 
      for(j=0;j<p->m;++j) a[j]+=w*(*uv++), b[j]+=w*(*ov++);
    else
      for(j=0;j<p->m;++j) a[j]+=w*(*ov++), b[j]+=w*(*uv++);
  }
}

/* work holds p->m * 2 doubles */
static void lob_bnd_1(const lob_bnd_base *p, const real *u, real bnd[2],
                      real *work)
{
  unsigned j;
  real *a = work, *b = work+p->m;
  lob_bnd_lines(p,u,a,b);
  bnd[0] = a[0], bnd[1] = b[0];
  for(j=1;j<p->m;++j) {
    if(a[j]<bnd[0]) bnd[0]=a[j];
    if(b[j]>bnd[1]) bnd[1]=b[j];
  }
}

/* work holds 2*mr + 2*mr*ns + 2*mr*ms doubles */
static void lob_bnd_2(const lob_bnd_base *pr, const lob_bnd_ext *ps,
                      const real *u, real bnd[2], real *work)
{
  unsigned nr = pr->n, mr = pr->m, ns = ps->b.n, ms = ps->b.m;
  real *a0 = work, *a1 = a0+mr,
       *ar_= a1+mr, *ar=ar_,
       *br_= ar+mr*ns, *br=br_,
       *a_ = br+mr*ns, *a =a_,
       *b_ = a +mr*ms, *b =b_,
       *uvp,*ovp,*uvn,*ovn;
  real b0,b1;
  unsigned i,j,k;
  for(i=0;i<mr;++i) a0[i]=a1[i]=0;
  for(j=0;j<ns;++j,u+=nr) {
    real q0 = ps->b.Q0[j], q1 = ps->b.Q1[j];
    lob_bnd_lines(pr,u,ar,br);
    for(i=0;i<mr;++i) {
      real t = (*ar++ + *br++)/2;
      a0[i]+=q0*t, a1[i]+=q1*t;
    }
  }
  for(i=0;i<mr;++i) {
    real a0i = a0[i], a1i = a1[i];
    for(k=0;k<ms;++k)
      *b++ = *a++ = a0i + a1i*ps->b.h[k];
  }
  ar = ar_, br = br_;
  uvp=ps->uvp, ovp=ps->ovp, uvn=ps->uvn, ovn=ps->ovn;
  for(j=0;j<ns;++j,uvp+=ms,uvn+=ms,ovp+=ms,ovn+=ms) {
    real zj = ps->b.z[j];
    a = a_, b = b_;
    for(i=0;i<mr;++i) {
      real t = a0[i] + a1[i]*zj;
      real uw = *ar++ - t;
      real ow = *br++ - t;
      if(uw>=0)      /* 0  <= uw <= ow */
        for(k=0;k<ms;++k)
          *a++ += uw * uvp[k] + ow * uvn[k],
          *b++ += ow * ovp[k] + uw * ovn[k];
      else if(ow<=0) /* uw <= ow <= 0  */
        for(k=0;k<ms;++k)
          *a++ += uw * ovp[k] + ow * ovn[k],
          *b++ += ow * uvp[k] + uw * uvn[k];
      else           /* uw <  0  <  ow */
        for(k=0;k<ms;++k)
          *a++ += uw * ovp[k] + ow * uvn[k],
          *b++ += ow * ovp[k] + uw * uvn[k];
    }
  }
  b0 = a_[0], b1 = b_[0];
  for(i=1;i<mr*ms;++i) {
    if(a_[i]<b0) b0=a_[i];
    if(b_[i]>b1) b1=b_[i];
  }
  bnd[0] = b0, bnd[1] = b1;
}

/*--------------------------------------------------------------------------
   Small Matrix Inverse 
  --------------------------------------------------------------------------*/

static void mat_inv_2(const real A[4], real inv[4])
{
  const real idet = 1/(A[0]*A[3]-A[1]*A[2]);
  inv[0] =   idet*A[3];
  inv[1] = -(idet*A[1]);
  inv[2] = -(idet*A[2]);
  inv[3] =   idet*A[0];
}

static void mat_inv_3(const real A[9], real inv[9])
{
  const real a = A[4]*A[8]-A[5]*A[7],
             b = A[5]*A[6]-A[3]*A[8],
             c = A[3]*A[7]-A[4]*A[6],
          idet = 1/(A[0]*a+A[1]*b+A[2]*c);
  inv[0] = idet*a;
  inv[1] = idet*(A[2]*A[7]-A[1]*A[8]);
  inv[2] = idet*(A[1]*A[5]-A[2]*A[4]);
  inv[3] = idet*b;
  inv[4] = idet*(A[0]*A[8]-A[2]*A[6]);
  inv[5] = idet*(A[2]*A[3]-A[0]*A[5]);
  inv[6] = idet*c;
  inv[7] = idet*(A[1]*A[6]-A[0]*A[7]);
  inv[8] = idet*(A[0]*A[4]-A[1]*A[3]);
}

static void mat_app_2r(real y[2], const real A[4], const real x[2])
{
  y[0] = A[0]*x[0] + A[1]*x[1];
  y[1] = A[2]*x[0] + A[3]*x[1];
}

static void mat_app_2c(real y[2], const real A[4], const real x[2])
{
  y[0] = A[0]*x[0] + A[2]*x[1];
  y[1] = A[1]*x[0] + A[3]*x[1];
}

static void mat_app_3r(real y[3], const real A[9], const real x[3])
{
  y[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2];
  y[1] = A[3]*x[0] + A[4]*x[1] + A[5]*x[2];
  y[2] = A[6]*x[0] + A[7]*x[1] + A[8]*x[2];
}

static void mat_app_3c(real y[3], const real A[9], const real x[3])
{
  y[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
  y[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
  y[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2];
}

static void tinyla_solve_2(real x[2], const real A[4], const real b[2])
{
  real inv[4];
  mat_inv_2(A,inv);
  mat_app_2r(x,inv,b);
}

static void tinyla_solve_3(real x[3], const real A[9], const real b[3])
{
  real inv[9];
  mat_inv_3(A,inv);
  mat_app_3r(x,inv,b);
}

/* solve
   A[0] x0 + A[2] x1 = b0,
   A[2] x0 + A[1] x1 = b1
*/
static void tinyla_solve_sym_2(real *x0, real *x1, const real A[3],
                               real b0, real b1)
{
  const real idet = 1/(A[0]*A[1] - A[2]*A[2]);
  *x0 = idet * (A[1]*b0 - A[2]*b1);
  *x1 = idet * (A[0]*b1 - A[2]*b0);
}

/*--------------------------------------------------------------------------
   Oriented Bounding Box

   Suffixes on names are _2 for 2-d and _3 for 3-d

   Needed inputs are the Gauss-Lobatto quadrature nodes and weights:   
     unsigned nr = ..., ns = ...;
     real zr[nr], wr[nr];
     real zs[ns], ws[ns];
     
     lobatto_nodes(zr,nr); lobatto_weights(zr,wr,nr);
     lobatto_nodes(zs,ns); lobatto_weights(zs,ws,ns);

   The number of points in the constructed piecewise (bi-)linear bounds
   for the boundaries is a parameter; more points give tighter bounds
   
     unsigned mr = 2*nr, ms = 2*ns;
  
   Bounding boxes are increased by a relative amount as a parameter
     
     real tol = 0.01;

   Setup is accomplished via:

     const real *z[2] = {zr,zs}, *w[2] = {wr,ws};
     const unsigned n[2] = {nr,ns}, m[2] = {mr,ms};
     obbox_data_2 *data = obbox_setup_2(z,w,n,m);
   
   Bounding box data may then be computed:
     
     obbox_2 box;                 // will store bounding box information
     real xm[ns][nr], ym[ns][nr]; // x, y coordinates of the element nodes
     
     obbox_calc_2(data, tol, (const real *)&xm[0][0],
                             (const real *)&ym[0][0], &box);
   
   A point may be tested:
   
     const real x[2]; // point to test
     real r[2];
     
     if( obbox_axis_test_2(&box, x) )
       ... // x failed axis-aligned bounding box test
       
     if( obbox_test_2(&box, x, r) )
       ... // x failed oriented bounding box test
     else
       ... // r suitable as initial guess for parametric coords
  
   Once all bounding box information has been computed
  
     obbox_free_2(data);
   
   to free the memory allocated with obbox_setup_2.

  --------------------------------------------------------------------------*/

typedef struct {
  lob_bnd_base dr, ds;
  real *Jr0, *Dr0, *Js0, *Ds0, *work;
} obbox_data_2;

typedef struct {
  lob_bnd_base dr;
  lob_bnd_ext ds, dt;
  real *Jr0, *Dr0, *Js0, *Ds0, *Jt0, *Dt0, *work;
} obbox_data_3;

static void obbox_data_alloc_2(obbox_data_2 *p,
                               const unsigned n[2], const unsigned m[2])
{
  const unsigned max_npm = umax_2(n[0]+m[0],n[1]+m[1]);
  lob_bnd_base_alloc(&p->dr, n[0], m[0]);
  lob_bnd_base_alloc(&p->ds, n[1], m[1]);
  p->Jr0 = tmalloc(real,2*n[0]+2*n[1]+2*max_npm);
  p->Dr0 = p->Jr0 + n[0];
  p->Js0 = p->Dr0 + n[0];
  p->Ds0 = p->Js0 + n[1];
  p->work = p->Ds0 + n[1];
}

static void obbox_data_free_2(obbox_data_2 *p)
{
  lob_bnd_base_free(&p->dr);
  lob_bnd_base_free(&p->ds);
  free(p->Jr0);
}

static void obbox_data_alloc_3(obbox_data_3 *p,
                               const unsigned n[3], const unsigned m[3])
{
  const unsigned wk1 = 3*n[0]*n[1] + 2*m[0] + 2*m[0]*n[1] + 2*m[0]*m[1];
  const unsigned wk2 = 3*n[0]*n[2] + 2*m[0] + 2*m[0]*n[2] + 2*m[0]*m[2];
  const unsigned wk3 = 3*n[1]*n[2] + 2*m[1] + 2*m[1]*n[2] + 2*m[1]*m[2];
  const unsigned wk_max = umax_3(wk1,wk2,wk3);
  lob_bnd_base_alloc(&p->dr, n[0], m[0]);
  lob_bnd_ext_alloc(&p->ds, n[1], m[1]);
  lob_bnd_ext_alloc(&p->dt, n[2], m[2]);
  p->Jr0 = tmalloc(real,2*n[0]+2*n[1]+2*n[2] + wk_max);
  p->Dr0 = p->Jr0 + n[0];
  p->Js0 = p->Dr0 + n[0];
  p->Ds0 = p->Js0 + n[1];
  p->Jt0 = p->Ds0 + n[1];
  p->Dt0 = p->Jt0 + n[2];
  p->work = p->Dt0 + n[2];
}

static void obbox_data_free_3(obbox_data_3 *p)
{
  lob_bnd_base_free(&p->dr);
  lob_bnd_ext_free(&p->ds);
  lob_bnd_ext_free(&p->dt);
  free(p->Jr0);
}

static obbox_data_2 *obbox_setup_2(const real *const z[2],
                                   const real *const w[2],
                                   const unsigned n[2], const unsigned m[2])
{
  const real zero = 0;
  real *work;
  obbox_data_2 *p = tmalloc(obbox_data_2,1);
  obbox_data_alloc_2(p,n,m);
  lob_bnd_base_setup(&p->dr,z[0],w[0]);
  lob_bnd_base_setup(&p->ds,z[1],w[1]);
  work = tmalloc(real,6*umax_2(n[0],n[1]));
  lagrange_weights_deriv(z[0],n[0],&zero,1,p->Jr0,p->Dr0,work);
  lagrange_weights_deriv(z[1],n[1],&zero,1,p->Js0,p->Ds0,work);
  free(work);
  return p;
}

static obbox_data_3 *obbox_setup_3(const real *const z[3],
                                   const real *const w[3],
                                   const unsigned n[3], const unsigned m[3])
{
  const real zero = 0;
  real *work;
  obbox_data_3 *p = tmalloc(obbox_data_3,1);
  obbox_data_alloc_3(p,n,m);
  lob_bnd_base_setup(&p->dr,z[0],w[0]);
  lob_bnd_ext_setup(&p->ds,z[1],w[1]);
  lob_bnd_ext_setup(&p->dt,z[2],w[2]);
  work = tmalloc(real,6*umax_3(n[0],n[1],n[2]));
  lagrange_weights_deriv(z[0],n[0],&zero,1,p->Jr0,p->Dr0,work);
  lagrange_weights_deriv(z[1],n[1],&zero,1,p->Js0,p->Ds0,work);
  lagrange_weights_deriv(z[2],n[2],&zero,1,p->Jt0,p->Dt0,work);
  free(work);
  return p;
}

static void obbox_free_2(obbox_data_2 *p)
{
  obbox_data_free_2(p);
  free(p);
}

static void obbox_free_3(obbox_data_3 *p)
{
  obbox_data_free_3(p);
  free(p);
}

typedef struct {
  real x[2], A[4], axis_bnd[4];
} obbox_2;

typedef struct {
  real x[3], A[9], axis_bnd[6];
} obbox_3;

static int obbox_axis_test_2(const obbox_2 *p, const real x[2])
{
  return (x[0]<p->axis_bnd[0] || x[0]>p->axis_bnd[1] ||
          x[1]<p->axis_bnd[2] || x[1]>p->axis_bnd[3]);
}

static int obbox_axis_test_3(const obbox_3 *p, const real x[3])
{
  return (x[0]<p->axis_bnd[0] || x[0]>p->axis_bnd[1] ||
          x[1]<p->axis_bnd[2] || x[1]>p->axis_bnd[3] ||
          x[2]<p->axis_bnd[4] || x[2]>p->axis_bnd[5]);
}

static int obbox_test_2(const obbox_2 *p, const real x[2], real r[2])
{
  const real xt[2] = {x[0]-p->x[0],x[1]-p->x[1]};
  r[0] = p->A[0]*xt[0] + p->A[1]*xt[1];
  if(fabsr(r[0])>1) return 1;
  r[1] = p->A[2]*xt[0] + p->A[3]*xt[1];
  return fabsr(r[1])>1;
}

static int obbox_test_3(const obbox_3 *p, const real x[3], real r[3])
{
  const real xt[3] = {x[0]-p->x[0],x[1]-p->x[1],x[2]-p->x[2]};
  r[0] = p->A[0]*xt[0] + p->A[1]*xt[1] + p->A[2]*xt[2];
  if(fabsr(r[0])>1) return 1;
  r[1] = p->A[3]*xt[0] + p->A[4]*xt[1] + p->A[5]*xt[2];
  if(fabsr(r[1])>1) return 1;
  r[2] = p->A[6]*xt[0] + p->A[7]*xt[1] + p->A[8]*xt[2];
  return fabsr(r[2])>1;
}

static void obbox_calc_tfm_2(const real *x, const real *y,
                             unsigned n, unsigned s,
                             const real c0[2], const real A[4], real *u)
{
  unsigned i;
  real *v = u+n;
  for(i=0; i<n; ++i,x+=s,y+=s) {
    const real xt = *x-c0[0], yt = *y-c0[1];
    u[i] = A[0]*xt + A[1]*yt;
    v[i] = A[2]*xt + A[3]*yt;
  }
}

static void obbox_calc_tfm_3(const real *x, const real *y, const real *z,
                             unsigned nr, unsigned sr, unsigned ns, unsigned ss,
                             const real c0[3], const real A[9], real *u)
{
  unsigned i,j;
  real *v = u+nr*ns, *w=v+nr*ns;
  for(j=0; j<ns; ++j,x+=ss,y+=ss,z+=ss) {
    for(i=0; i<nr; ++i,x+=sr,y+=sr,z+=sr) {
      const real xt = *x-c0[0], yt = *y-c0[1], zt = *z-c0[2];
      *u++ = A[0]*xt + A[1]*yt + A[2]*zt;
      *v++ = A[3]*xt + A[4]*yt + A[5]*zt;
      *w++ = A[6]*xt + A[7]*yt + A[8]*zt;
    }
  }
}

static void obbox_merge_2(real *b, const real *ob)
{
  if(ob[0]<b[0]) b[0]=ob[0];
  if(ob[1]>b[1]) b[1]=ob[1];
  if(ob[2]<b[2]) b[2]=ob[2];
  if(ob[3]>b[3]) b[3]=ob[3];
}

static void obbox_merge_3(real *b, const real *ob)
{
  if(ob[0]<b[0]) b[0]=ob[0];
  if(ob[1]>b[1]) b[1]=ob[1];
  if(ob[2]<b[2]) b[2]=ob[2];
  if(ob[3]>b[3]) b[3]=ob[3];
  if(ob[4]<b[4]) b[4]=ob[4];
  if(ob[5]>b[5]) b[5]=ob[5];
}

/* work holds 2*n + 2*m reals */
static void obbox_side_2(const real *x, const real *y,
                         unsigned n, unsigned s,
                         const real c0[2], const real A[4], real *work,
                         const lob_bnd_base *lbd, real bnd[4])
{
  obbox_calc_tfm_2(x,y,n,s,c0,A,work);
  lob_bnd_1(lbd,work  ,bnd  ,work+2*n);
  lob_bnd_1(lbd,work+n,bnd+2,work+2*n);
}

/* work holds 3*nr*ns + 2*mr + 2*mr*ns + 2*mr*ms reals */
static void obbox_side_3(const real *x, const real *y, const real *z,
                         unsigned nr, unsigned sr, unsigned ns, unsigned ss,
                         const real c0[3], const real A[9], real *work,
                         const lob_bnd_base *dr, const lob_bnd_ext *ds,
                         real bnd[6])
{
  obbox_calc_tfm_3(x,y,z,nr,sr,ns,ss,c0,A,work);
  lob_bnd_2(dr,ds,work        ,bnd  ,work+3*nr*ns);
  lob_bnd_2(dr,ds,work+  nr*ns,bnd+2,work+3*nr*ns);
  lob_bnd_2(dr,ds,work+2*nr*ns,bnd+4,work+3*nr*ns);
}

/* return bounds on u = A (x - c0)
   bnd[0] <= u_0 <= bnd[1]
   bnd[2] <= u_1 <= bnd[3] */
static void obbox_bnd_2(const obbox_data_2 *p,
                        const real *x, const real *y,
                        const real c0[2], const real A[4],
                        real bnd[4])
{
  unsigned i, nr = p->dr.n, ns = p->ds.n;
  real obnd[4];

  i = nr*(ns-1);
  obbox_side_2(x  ,y  , nr, 1, c0,A,p->work, &p->dr, bnd);
  obbox_side_2(x+i,y+i, nr, 1, c0,A,p->work, &p->dr, obnd);
  obbox_merge_2(bnd,obnd);

  i = nr-1;
  obbox_side_2(x  ,y  , ns,nr, c0,A,p->work, &p->ds, obnd);
  obbox_merge_2(bnd,obnd);
  obbox_side_2(x+i,y+i, nr,nr, c0,A,p->work, &p->ds, obnd);
  obbox_merge_2(bnd,obnd);
}

/* return bounds on u = A (x - c0)
   bnd[0] <= u_0 <= bnd[1]
   bnd[2] <= u_1 <= bnd[3]
   bnd[4] <= u_2 <= bnd[5] */
static void obbox_bnd_3(const obbox_data_3 *p,
                        const real *x, const real *y, const real *z,
                        const real c0[3], const real A[9],
                        real bnd[6])
{
  unsigned i, nr = p->dr.n, ns = p->ds.b.n, nt = p->dt.b.n;
  real obnd[6];

  i = nr*ns*(nt-1);
  obbox_side_3(x  ,y  ,z  , nr, 1,ns,0, c0,A,p->work, &p->dr  ,&p->ds, bnd);
  obbox_side_3(x+i,y+i,z+i, nr, 1,ns,0, c0,A,p->work, &p->dr  ,&p->ds, obnd);
  obbox_merge_3(bnd,obnd);

  i = nr*(ns-1);
  obbox_side_3(x  ,y  ,z  , nr, 1,nt,i, c0,A,p->work, &p->dr  ,&p->dt, obnd);
  obbox_merge_3(bnd,obnd);
  obbox_side_3(x+i,y+i,z+i, nr, 1,nt,i, c0,A,p->work, &p->dr  ,&p->dt, obnd);
  obbox_merge_3(bnd,obnd);

  i = nr-1;
  obbox_side_3(x  ,y  ,z  , ns,nr,nt,0, c0,A,p->work, &p->ds.b,&p->dt, obnd);
  obbox_merge_3(bnd,obnd);
  obbox_side_3(x+i,y+i,z+i, ns,nr,nt,0, c0,A,p->work, &p->ds.b,&p->dt, obnd);
  obbox_merge_3(bnd,obnd);
}

static void obbox_calc_2(const obbox_data_2 *p, real tol,
                         const real *x, const real *y, obbox_2 *b)
{
  const real zero[2] = {0,0}, id[4] = {1,0,0,1};
  real c0[2], jac[4], inv[4], bnd[4], u0[2], d[2];
  
  obbox_bnd_2(p,x,y,zero,id,b->axis_bnd);
  d[0] = b->axis_bnd[1]-b->axis_bnd[0];
  d[1] = b->axis_bnd[3]-b->axis_bnd[2];
  b->axis_bnd[0] -= tol*d[0], b->axis_bnd[1] += tol*d[0];
  b->axis_bnd[2] -= tol*d[1], b->axis_bnd[3] += tol*d[1];
  
  c0[0] = tensor_ig2(p->Jr0,p->Dr0,p->dr.n,
                     p->Js0,p->Ds0,p->ds.n,
                     x, jac  , p->work);
  c0[1] = tensor_ig2(p->Jr0,p->Dr0,p->dr.n,
                     p->Js0,p->Ds0,p->ds.n,
                     y, jac+2, p->work);
  mat_inv_2(jac,inv);
  
  obbox_bnd_2(p,x,y,c0,inv,bnd);
  
  u0[0] = (bnd[0]+bnd[1])/2;
  u0[1] = (bnd[2]+bnd[3])/2;
  d[0] = 2/((1+tol)*(bnd[1]-bnd[0]));
  d[1] = 2/((1+tol)*(bnd[3]-bnd[2]));
  b->x[0] = c0[0] + jac[0]*u0[0] + jac[1]*u0[1];
  b->x[1] = c0[1] + jac[2]*u0[0] + jac[3]*u0[1];
  b->A[0] = d[0]*inv[0], b->A[1] = d[0]*inv[1];
  b->A[2] = d[1]*inv[2], b->A[3] = d[1]*inv[3];
}

static void obbox_calc_3(const obbox_data_3 *p, real tol,
                         const real *x, const real *y, const real *z,
                         obbox_3 *b)
{
  const real zero[3] = {0,0}, id[9] = {1,0,0,0,1,0,0,0,1};
  real c0[3], jac[9], inv[9], bnd[6], u0[3], d[3];
  
  obbox_bnd_3(p,x,y,z,zero,id,b->axis_bnd);
  d[0] = b->axis_bnd[1]-b->axis_bnd[0];
  d[1] = b->axis_bnd[3]-b->axis_bnd[2];
  d[2] = b->axis_bnd[5]-b->axis_bnd[4];
  b->axis_bnd[0] -= tol*d[0], b->axis_bnd[1] += tol*d[0];
  b->axis_bnd[2] -= tol*d[1], b->axis_bnd[3] += tol*d[1];
  b->axis_bnd[4] -= tol*d[2], b->axis_bnd[5] += tol*d[2];
  
  c0[0] = tensor_ig3(p->Jr0,p->Dr0,p->dr.n,
                     p->Js0,p->Ds0,p->ds.b.n,
                     p->Jt0,p->Dt0,p->dt.b.n,
                     x, jac  , p->work);
  c0[1] = tensor_ig3(p->Jr0,p->Dr0,p->dr.n,
                     p->Js0,p->Ds0,p->ds.b.n,
                     p->Jt0,p->Dt0,p->dt.b.n,
                     y, jac+3, p->work);
  c0[2] = tensor_ig3(p->Jr0,p->Dr0,p->dr.n,
                     p->Js0,p->Ds0,p->ds.b.n,
                     p->Jt0,p->Dt0,p->dt.b.n,
                     z, jac+6, p->work);
  mat_inv_3(jac,inv);

  obbox_bnd_3(p,x,y,z,c0,inv,bnd);  
  
  u0[0] = (bnd[0]+bnd[1])/2;
  u0[1] = (bnd[2]+bnd[3])/2;
  u0[2] = (bnd[4]+bnd[5])/2;
  d[0] = 2/((1+tol)*(bnd[1]-bnd[0]));
  d[1] = 2/((1+tol)*(bnd[3]-bnd[2]));
  d[2] = 2/((1+tol)*(bnd[5]-bnd[4]));
  b->x[0] = c0[0] + jac[0]*u0[0] + jac[1]*u0[1] + jac[2]*u0[2];
  b->x[1] = c0[1] + jac[3]*u0[0] + jac[4]*u0[1] + jac[5]*u0[2];
  b->x[2] = c0[2] + jac[6]*u0[0] + jac[7]*u0[1] + jac[8]*u0[2];
  b->A[0] = d[0]*inv[0], b->A[1] = d[0]*inv[1], b->A[2] = d[0]*inv[2];
  b->A[3] = d[1]*inv[3], b->A[4] = d[1]*inv[4], b->A[5] = d[1]*inv[5];
  b->A[6] = d[2]*inv[6], b->A[7] = d[2]*inv[7], b->A[8] = d[2]*inv[8];
}

/*--------------------------------------------------------------------------
   Point to Possible Elements Hashing
   
   Initializing the data:
     unsigned nel;        // number of elements
     const unsigned n[3]; // number of nodes in r, s, t directions
     const real *xm[3];   // n[0]*n[1]*n[2]*nel x,y,z coordinates
     real tol = 0.01;     // how far point is allowed to be outside element
                          //   relative to element size
     unsigned max_size = n[0]*n[1]*n[2]*nel; // maximum size of hash table
     
     hash_data_3 data;
     hash_build_3(&data, xm, n, nel, max_size, tol);
     
   Using the data:
     real x[3];   // point to find
     
     unsigned index = hash_index_3(&data, x);
     unsigned i, b = data.offset[index], e = data.offset[index+1];
     
     // point may be in elements
     //   data.offset[b], data.offset[b+1], ... , data.offset[e-1]
     //
     // list has maximum size data.max (e.g., e-b <= data.max)
   
     for(i=b; i!=e; ++i) {
       unsigned el = data.offset[i];
       const obbox_3 *obb = &data.obb[el]; // bounding box data for element el
       ...
     }
   
   When done:
     hash_free_3(&data);
     
  --------------------------------------------------------------------------*/

typedef struct {
  unsigned hash_n;
  real bnd[4]; /* bounds for all elements */
  real fac[2]; /* fac[i] = hash_n / (bnd[2*i+1]-bnd[2*i]) */
  obbox_2 *obb; /* obb[nel] -- bounding box info for each element */
  uint *offset; /* hash table -- for cell i,j:
                         uint index = j*hash_n+i,
                                  b = offset[index  ],
                                  e = offset[index+1];
                         elements in cell are
                           offset[b], offset[b+1], ..., offset[e-1] */
  unsigned max; /* maximum # of elements in any cell */
} hash_data_2;

typedef struct {
  unsigned hash_n;
  real bnd[6]; /* bounds for all elements */
  real fac[3]; /* fac[i] = hash_n / (bnd[2*i+1]-bnd[2*i]) */
  obbox_3 *obb; /* obb[nel] -- bounding box info for each element */
  uint *offset; /* hash table -- for cell i,j,k:
                         uint index = (k*hash_n+j)*hash_n+i,
                                  b = offset[index  ],
                                  e = offset[index+1];
                         elements in cell are
                           offset[b], offset[b+1], ..., offset[e-1] */
  unsigned max; /* maximum # of elements in any cell */
} hash_data_3;

static int ifloor(real x)
{
  /*
  int y = x;
  return (double)y > x ? y-1 : y;
  */
  return floorr(x);
}

static int iceil(real x)
{
  /*
  int y = x;
  return (double)y < x ? y+1 : y;
  */
  return ceilr(x);
}

static unsigned hash_index_helper(real low, real fac, unsigned n, real x)
{
  const int i = ifloor((x-low)*fac);
  if(i<0) return 0;
  return umin_2(i,n-1);
}

static uint hash_index_2(const hash_data_2 *p, const real x[2])
{
  const unsigned n = p->hash_n;
  return (uint)hash_index_helper(p->bnd[2],p->fac[1],n,x[1])*n
              +hash_index_helper(p->bnd[0],p->fac[0],n,x[0]);
}

static uint hash_index_3(const hash_data_3 *p, const real x[3])
{
  const unsigned n = p->hash_n;
  return ( (uint)hash_index_helper(p->bnd[4],p->fac[2],n,x[2])  *n
                +hash_index_helper(p->bnd[2],p->fac[1],n,x[1]) )*n
                +hash_index_helper(p->bnd[0],p->fac[0],n,x[0]);
}

static void hash_setfac_2(hash_data_2 *p, unsigned n)
{
  p->hash_n = n;
  p->fac[0] = n/(p->bnd[1] - p->bnd[0]);
  p->fac[1] = n/(p->bnd[3] - p->bnd[2]);
}

static void hash_setfac_3(hash_data_3 *p, unsigned n)
{
  p->hash_n = n;
  p->fac[0] = n/(p->bnd[1] - p->bnd[0]);
  p->fac[1] = n/(p->bnd[3] - p->bnd[2]);
  p->fac[2] = n/(p->bnd[5] - p->bnd[4]);
}

static void hash_range_2(const hash_data_2 *p, uint i, unsigned d,
                         unsigned *ia, unsigned *ib)
{
  const real a = p->obb[i].axis_bnd[d*2  ];
  const real b = p->obb[i].axis_bnd[d*2+1];
  const int      i0 = ifloor( (a - p->bnd[d*2]) * p->fac[d] );
  const unsigned i1 = iceil(  (b - p->bnd[d*2]) * p->fac[d] );
  *ia = imax_2(0,i0);
  *ib = imin_2(i1,p->hash_n);
  if(*ib == *ia) ++(*ib);
}

static void hash_range_3(const hash_data_3 *p, uint i, unsigned d,
                         unsigned *ia, unsigned *ib)
{
  const real a = p->obb[i].axis_bnd[d*2  ];
  const real b = p->obb[i].axis_bnd[d*2+1];
  const int      i0 = ifloor( (a - p->bnd[d*2]) * p->fac[d] );
  const unsigned i1 = iceil(  (b - p->bnd[d*2]) * p->fac[d] );
  *ia = imax_2(0,i0);
  *ib = imin_2(i1,p->hash_n);
  if(*ib == *ia) ++(*ib);
}

static uint hash_count_2(hash_data_2 *p, uint nel, unsigned n)
{
  uint i,count=0;
  hash_setfac_2(p,n);
  for(i=0;i<nel;++i) {
    unsigned ia, ib; uint ci;
    hash_range_2(p,i,0,&ia,&ib); ci  = ib-ia;
    hash_range_2(p,i,1,&ia,&ib); ci *= ib-ia;
    count+=ci;
  }
  return count;
}

static uint hash_count_3(hash_data_3 *p, uint nel, unsigned n)
{
  uint i,count=0;
  hash_setfac_3(p,n);
  for(i=0;i<nel;++i) {
    unsigned ia, ib; uint ci;
    hash_range_3(p,i,0,&ia,&ib); ci  = ib-ia;
    hash_range_3(p,i,1,&ia,&ib); ci *= ib-ia;
    hash_range_3(p,i,2,&ia,&ib); ci *= ib-ia;
    count+=ci;
  }
  return count;
}

static uint hash_opt_size_2(hash_data_2 *p, uint nel, uint max_size)
{
  unsigned nl=1, nu = ceil(sqrt(max_size-nel));
  uint size_low = 2+nel;
  while(nu-nl>1) {
    unsigned nm = nl+(nu-nl)/2;
    uint size = (uint)nm*nm+1+hash_count_2(p,nel,nm);
    if(size<=max_size) nl=nm,size_low=size; else nu=nm;
  }
  hash_setfac_2(p,nl);
  return size_low;
}

static uint hash_opt_size_3(hash_data_3 *p, uint nel, uint max_size)
{
  unsigned nl=1, nu = ceil(pow(max_size-nel,1.0/3));
  uint size_low = 2+nel;
  while(nu-nl>1) {
    unsigned nm = nl+(nu-nl)/2;
    uint size = (uint)nm*nm*nm+1+hash_count_3(p,nel,nm);
    if(size<=max_size) nl=nm,size_low=size; else nu=nm;
  }
  hash_setfac_3(p,nl);
  return size_low;
}

static void hash_getbb_2(hash_data_2 *p, const real *const elx[2],
                         const unsigned n[2], uint nel, real tol)
{
  obbox_data_2 *data;
  const real *x[2]={elx[0],elx[1]};
  real *z[2], *w[2];
  uint i; unsigned d;
  const unsigned nn = n[0]*n[1], m[2] = {2*n[0],2*n[1]};
  
  z[0] = tmalloc(real,2*(n[0]+n[1]));
  w[0] = z[0] + n[0];
  z[1] = w[0] + n[0], w[1] = z[1] + n[1];
  for(d=0;d<2;++d)
    lobatto_nodes(z[d],n[d]), lobatto_weights(z[d],w[d],n[d]);
  data = obbox_setup_2((const real *const*)z,(const real *const*)w,n,m);
  obbox_calc_2(data,tol,x[0],x[1],&p->obb[0]);
  memcpy(&p->bnd[0],(const real*)&p->obb[0].axis_bnd[0],4*sizeof(real));
  for(i=0;i<nel;++i,x[0]+=nn,x[1]+=nn) {
    obbox_calc_2(data,tol,x[0],x[1],&p->obb[i]);
    obbox_merge_2(&p->bnd[0],(const real*)&p->obb[i].axis_bnd[0]);
  }
  obbox_free_2(data);
  free(z[0]);
}

static void hash_getbb_3(hash_data_3 *p, const real *const elx[3],
                         const unsigned n[3], uint nel, real tol)
{
  obbox_data_3 *data;
  const real *x[3]={elx[0],elx[1],elx[2]};
  real *z[3], *w[3];
  uint i; unsigned d;
  const unsigned nn = n[0]*n[1]*n[2], m[3] = {2*n[0],2*n[1],2*n[2]};
  
  z[0] = tmalloc(real,2*(n[0]+n[1]+n[2]));
  w[0] = z[0] + n[0];
  for(d=1;d<3;++d) z[d]=w[d-1]+n[d-1], w[d]=z[d]+n[d];
  for(d=0;d<3;++d)
    lobatto_nodes(z[d],n[d]), lobatto_weights(z[d],w[d],n[d]);
  data = obbox_setup_3((const real *const*)z,(const real *const*)w,n,m);
  obbox_calc_3(data,tol,x[0],x[1],x[2],&p->obb[0]);
  memcpy(&p->bnd[0],(const real*)&p->obb[0].axis_bnd[0],6*sizeof(real));
  for(i=0;i<nel;++i,x[0]+=nn,x[1]+=nn,x[2]+=nn) {
    obbox_calc_3(data,tol,x[0],x[1],x[2],&p->obb[i]);
    obbox_merge_3(&p->bnd[0],(const real*)&p->obb[i].axis_bnd[0]);
  }
  obbox_free_3(data);
  free(z[0]);
}

static void hash_build_2(hash_data_2 *p, const real *const x[2],
                         const unsigned n[2], uint nel,
                         uint max_hash_size, real tol)
{
  uint i,el,size,hn2,sum; unsigned hn;
  unsigned *count;
  p->obb = tmalloc(obbox_2,nel);
  hash_getbb_2(p,x,n,nel,tol);
  size = hash_opt_size_2(p,nel,max_hash_size);
  p->offset = tmalloc(uint,size);
  hn = p->hash_n;
  hn2 = (uint)hn*hn;
  count = tcalloc(unsigned,hn2);
  for(el=0;el<nel;++el) {
    unsigned i,ia,ib, j,ja,jb;
    hash_range_2(p,el,0,&ia,&ib);
    hash_range_2(p,el,1,&ja,&jb);
    for(j=ja;j<jb;++j) for(i=ia;i<ib;++i)
      ++count[(uint)j*hn+i];
  }
  sum=hn2+1, p->max=count[0];
  p->offset[0]=sum;
  for(i=0;i<hn2;++i) {
    if(count[i]>p->max) p->max=count[i];
    sum+=count[i];
    p->offset[i+1]=sum;
  }
  for(el=0;el<nel;++el) {
    unsigned i,ia,ib, j,ja,jb;
    hash_range_2(p,el,0,&ia,&ib);
    hash_range_2(p,el,1,&ja,&jb);
    for(j=ja;j<jb;++j) for(i=ia;i<ib;++i) {
      uint index = (uint)j*hn+i;
      p->offset[p->offset[index+1] - count[index]] = el;
      --count[index];
    }
  }
  free(count);
}

static void hash_build_3(hash_data_3 *p, const real *const x[3],
                         const unsigned n[3], uint nel,
                         uint max_hash_size, real tol)
{
  uint i,el,size,hn3,sum; unsigned hn;
  unsigned *count;
  p->obb = tmalloc(obbox_3,nel);
  hash_getbb_3(p,x,n,nel,tol);
  size = hash_opt_size_3(p,nel,max_hash_size);
  p->offset = tmalloc(uint,size);
  hn = p->hash_n;
  hn3 = (uint)hn*hn*hn;
  count = tcalloc(unsigned,hn3);
  for(el=0;el<nel;++el) {
    unsigned i,ia,ib, j,ja,jb, k,ka,kb;
    hash_range_3(p,el,0,&ia,&ib);
    hash_range_3(p,el,1,&ja,&jb);
    hash_range_3(p,el,2,&ka,&kb);
    for(k=ka;k<kb;++k) for(j=ja;j<jb;++j) for(i=ia;i<ib;++i)
      ++count[((uint)k*hn+j)*hn+i];
  }
  sum=hn3+1, p->max=count[0];
  p->offset[0]=sum;
  for(i=0;i<hn3;++i) {
    if(count[i]>p->max) p->max=count[i];
    sum+=count[i];
    p->offset[i+1]=sum;
  }
  for(el=0;el<nel;++el) {
    unsigned i,ia,ib, j,ja,jb, k,ka,kb;
    hash_range_3(p,el,0,&ia,&ib);
    hash_range_3(p,el,1,&ja,&jb);
    hash_range_3(p,el,2,&ka,&kb);
    for(k=ka;k<kb;++k) for(j=ja;j<jb;++j) for(i=ia;i<ib;++i) {
      uint index = ((uint)k*hn+j)*hn+i;
      p->offset[p->offset[index+1] - count[index]] = el;
      --count[index];
    }
  }
  free(count);
}

static void hash_free_2(hash_data_2 *p)
{
  free(p->obb);
  free(p->offset);
}

static void hash_free_3(hash_data_3 *p)
{
  free(p->obb);
  free(p->offset);
}

/*--------------------------------------------------------------------------
   Optimization algorithm to find a point within an element
          
   Given x(r)  (as values of x,y,z at all Lobatto nodes) and x_star,
   find the r that minimizes || x_star - x(r) ||_2
   
   As a minimization problem, the Newton step is
   
               __ 3
     [ J^T J - >_ d=1  resid_d H_d ] dr = J^t resid
     
   where resid = x_star - x(r), J = [ dx_i/dr_j ],
   and H_d = [ d^2 x_d/dr_i dr_j ].
  
   This is the appropriate step to take whenever constraints are active,
   and the current iterate is on a boundary of the element. When the current
   iterate is inside, J is square ( dim r = dim x ), resid will become small,
   and the step
  
     J dr = resid
  
   may be used instead, still giving quadratic convergence.

   
   Names use a _3 suffix for 3-d and _2 for 2-d.
   The routines require an initialized lagrange_data array as input:
     unsigned d, n[3] = { ... };
     real *z[3] = { tmalloc(real, n[0]), ... };
     for(d=0;d<3;++d) lobatto_nodes(z[d],n[d]);
     
     lagrange_data ld[3];
     for(d=0;d<3;++d) lagrange_setup(&ld[d],z[d],n[d]);
   
   Initialization:
     opt_data_3 data;
     opt_alloc_3(&data, ld);
  
   Use:
     const real *xm[3];  // 3 pointers, each to n[0]*n[1]*n[2] reals
                         //   giving the nodal x, y, or z coordinates
    
     const real x_star[3] = { ... };  // point to find
     real r[3] = { 0,0,0 };           // initial guess with
     unsigned c = opt_no_constraints_3;   //   these constraints active
     
     real dist = opt_findpt_3(&data,xm,x_star,r,&c);
     // minimizer is r with constraints c; 2-norm of resid is dist
   
   Clean-up:
     opt_free_3(&data);
     
     for(d=0;d<3;++d) lagrange_free(&ld[d]);
     for(d=0;d<3;++d) free(z[d]);
  
   The constraint number works as follows. Let cr be the constraints
   on the r variable:
      cr = 0       r fixed at -1
      cr = 1       r not fixed
      cr = 2       r fixed at 1
   Then the constraint number is (ct*3+cs)*3+cr
     
  --------------------------------------------------------------------------*/

static const unsigned opt_no_constraints_2 = 3+1;
static const unsigned opt_no_constraints_3 = 9+3+1;

/* how many directions are constrained? */
static const char opt_constr_num_2[9] = {2,1,2, 1,0,1, 2,1,2};
static const char opt_constr_num_3[27] = {
  3,2,3, 2,1,2, 3,2,3, 
  2,1,2, 1,0,1, 2,1,2, 
  3,2,3, 2,1,2, 3,2,3
};

/* which direction is constrained? */
static const char opt_constr_dir_2[9] = {-1, 1,-1,  0,-1, 0, -1, 1,-1};
static const char opt_constr_dir_3[27] = {
  -1,-1,-1, -1, 2,-1, -1,-1,-1, 
  -1, 1,-1,  0,-1, 0, -1, 1,-1, 
  -1,-1,-1, -1, 2,-1, -1,-1,-1
};

/* which direction is not constrained? */
static const char opt_constr_not[27] = {
  -1, 0,-1,  1,-1, 1, -1, 0,-1, 
   2,-1, 2, -1,-1,-1,  2,-1, 2, 
  -1, 0,-1,  1,-1, 1, -1, 0,-1
};

static const char opt_constr_wide[27] = {
  0x00,0x01,0x02, 0x04,0x05,0x06, 0x08,0x09,0x0a,
  0x10,0x11,0x12, 0x14,0x15,0x16, 0x18,0x19,0x1a,
  0x20,0x21,0x22, 0x24,0x25,0x26, 0x28,0x29,0x2a
};

static const unsigned opt_other1_3[3] = {1,0,0},
                      opt_other2_3[3] = {2,2,1};

static unsigned opt_constr(unsigned constraints, unsigned d)
{
  return (opt_constr_wide[constraints]>>(d*2))&3;
}

static void opt_constr_unpack_2(unsigned constraints, unsigned *c)
{
  const char cw = opt_constr_wide[constraints];
  c[0] = cw & 3;
  c[1] = cw >> 2;
}

static void opt_constr_unpack_3(unsigned constraints, unsigned *c)
{
  const char cw = opt_constr_wide[constraints];
  c[0] = cw & 3;
  c[1] = (cw >> 2) & 3;
  c[2] = cw >> 4;
}

static unsigned opt_constr_pack_2(const unsigned *c)
{
  return c[1]*3+c[0];
}

static unsigned opt_constr_pack_3(const unsigned *c)
{
  return (c[2]*3+c[1])*3+c[0];
}

/*--------------------------------------------------------------------------
   
   3 - D
     
  --------------------------------------------------------------------------*/

typedef struct {
  unsigned constraints;
  unsigned dn, d1, d2;
  real *x[3], *fdn[3];
} opt_face_data_3;

typedef struct {
  unsigned constraints;
  unsigned de, d1, d2;
  real *x[3], *fd1[3], *fd2[3];
} opt_edge_data_3;

typedef struct {
  unsigned constraints;
  real x[3], jac[9];
} opt_point_data_3;

typedef struct {
  lagrange_data *ld;
  unsigned size[4];
  const real *elx[3];
  opt_face_data_3 fd;
  opt_edge_data_3 ed;
  opt_point_data_3 pd;
  real *work;
  real x[3], jac[9];
} opt_data_3;

static void opt_alloc_3(opt_data_3 *p, lagrange_data *ld)
{
  const unsigned nr = ld[0].n, ns = ld[1].n, nt = ld[2].n,
                 nf = umax_3(nr*ns,nr*nt,ns*nt),
                 ne = umax_3(nr,ns,nt),
                 nw = 2*ns*nt + 3*ns;
  p->size[0] = 1;
  p->size[1] = nr;
  p->size[2] = nr*ns;
  p->size[3] = p->size[2]*nt;
  p->ld = ld;
  p->work = tmalloc(real, 6*nf + 9*ne + nw);
  p->fd.x[0] = p->work + nw;
  p->fd.x[1] = p->fd.x[0] + nf; 
  p->fd.x[2] = p->fd.x[1] + nf; 
  p->fd.fdn[0] = p->fd.x[2] + nf; 
  p->fd.fdn[1] = p->fd.fdn[0] + nf; 
  p->fd.fdn[2] = p->fd.fdn[1] + nf; 
  p->ed.x[0] = p->fd.fdn[2] + nf; 
  p->ed.x[1] = p->ed.x[0] + ne; 
  p->ed.x[2] = p->ed.x[1] + ne; 
  p->ed.fd1[0] = p->ed.x[2] + ne; 
  p->ed.fd1[1] = p->ed.fd1[0] + ne; 
  p->ed.fd1[2] = p->ed.fd1[1] + ne; 
  p->ed.fd2[0] = p->ed.fd1[2] + ne; 
  p->ed.fd2[1] = p->ed.fd2[0] + ne; 
  p->ed.fd2[2] = p->ed.fd2[1] + ne; 
}

static void opt_free_3(opt_data_3 *p)
{
  free(p->work);
}

static void opt_vol_set_3(opt_data_3 *p, const real r[3])
{
  lagrange_1(&p->ld[0],r[0]);
  lagrange_1(&p->ld[1],r[1]);
  lagrange_1(&p->ld[2],r[2]);
}

/* work holds 2*ns*nt + 3*ns reals */
static void opt_vol_intp_3(opt_data_3 *p)
{
  unsigned d;
  const lagrange_data *ld = p->ld;
  
  for(d=0;d<3;++d) 
    p->x[d] = tensor_ig3(ld[0].J,ld[0].D,ld[0].n,
                         ld[1].J,ld[1].D,ld[1].n,
                         ld[2].J,ld[2].D,ld[2].n,
                         p->elx[d], &p->jac[d*3], p->work);
}

static void opt_vol_set_intp_3(opt_data_3 *p, const real r[3])
{
  opt_vol_set_3(p,r);
  opt_vol_intp_3(p);
}

static void opt_face_proj_3(opt_data_3 *p)
{
  unsigned d, off=0;
  const unsigned dn = p->fd.dn, d1 = p->fd.d1, d2 = p->fd.d2,
                 so = p->size[d2]-p->size[d1+1],
                 s1 = p->size[d1], sn = p->size[dn],
                 n1 = p->ld[d1].n, n2 = p->ld[d2].n, nn = p->ld[dn].n;
  const real *D = p->ld[dn].D_z0;
  if(opt_constr(p->fd.constraints,dn)==2)
    off = p->size[dn+1]-p->size[dn],
    D   = p->ld[dn].D_zn;
  for(d=0;d<3;++d) {
    unsigned i,j,k,index=0;
    const real *in = p->elx[d]+off;
    for(j=n2;j;--j,in+=so)
      for(i=n1;i;--i,++index,in+=s1) {
        const real *ind = in-off;
        real *fdn = &p->fd.fdn[d][index];
        p->fd.x[d][index] = *in;
        *fdn = 0;
        for(k=0;k<nn;++k,ind+=sn)
          *fdn += *ind * D[k];
      }
  }
}

static void opt_face_set_3(opt_data_3 *p, const real r[3], unsigned constr)
{
  if(p->fd.constraints!=constr) {
    p->fd.constraints=constr;
    p->fd.dn = opt_constr_dir_3[constr];
    p->fd.d1 = opt_other1_3[p->fd.dn];
    p->fd.d2 = opt_other2_3[p->fd.dn];
    opt_face_proj_3(p);
  }
  lagrange_1(&p->ld[p->fd.d1],r[p->fd.d1]);
  lagrange_1(&p->ld[p->fd.d2],r[p->fd.d2]);
}

/* work holds 2*ld[d2].n reals */
static void opt_face_intp_3(opt_data_3 *p)
{
  unsigned d;
  const unsigned dn = p->fd.dn, d1 = p->fd.d1, d2 = p->fd.d2,
                 n1 = p->ld[d1].n, n2 = p->ld[d2].n;
  const real *J1 = p->ld[d1].J, *J2 = p->ld[d2].J,
             *D1 = p->ld[d1].D, *D2 = p->ld[d2].D;
  
  for(d=0;d<3;++d) {
    real g[2];
    p->x[d] = tensor_ig2(J1,D1,n1, J2,D2,n2, p->fd.x[d], &g[0], p->work);
    p->jac[d*3+d1] = g[0];
    p->jac[d*3+d2] = g[1];
    p->jac[d*3+dn] = tensor_i2(J1,n1, J2,n2, p->fd.fdn[d], p->work);
  }
}

static void opt_face_set_intp_3(opt_data_3 *p, const real r[3], unsigned constr)
{
  opt_face_set_3(p,r,constr);
  opt_face_intp_3(p);
}

static void opt_face_hess_3(opt_data_3 *p, real hess[9])
{
  unsigned d;
  const unsigned d1 = p->fd.d1, d2 = p->fd.d2,
                 n1 = p->ld[d1].n, n2 = p->ld[d2].n;
  const real *J1 = p->ld[d1].J , *J2 = p->ld[d2].J,
             *D1 = p->ld[d1].D , *D2 = p->ld[d2].D,
             *S1 = p->ld[d1].D2, *S2 = p->ld[d2].D2;
  
  lagrange_2u(&p->ld[d1]);
  lagrange_2u(&p->ld[d2]);
  
  for(d=0;d<3;++d) {
    (void) tensor_ig2(J1,S1,n1, J2,S2,n2,  p->fd.x[d], hess+d*3, p->work);
    hess[d*3+0] = tensor_i2(S1,n1, J2,n2, p->fd.x[d], p->work);
    hess[d*3+1] = tensor_i2(J1,n1, S2,n2, p->fd.x[d], p->work);
    hess[d*3+2] = tensor_i2(D1,n1, D2,n2, p->fd.x[d], p->work);
  }
}

static void opt_edge_proj_3(opt_data_3 *p)
{
  unsigned d, off, off1=0, off2=0;
  const unsigned de=p->ed.de,    d1=p->ed.d1,    d2=p->ed.d2,
                 se=p->size[de], s1=p->size[d1], s2=p->size[d2],
                 ne=p->ld[de].n, n1=p->ld[d1].n, n2=p->ld[d2].n;
  const real *fD1, *fD2;
  if(opt_constr(p->ed.constraints,d1)==0)
    fD1=p->ld[d1].D_z0;
  else
    fD1=p->ld[d1].D_zn, off1 = p->size[d1+1]-p->size[d1];
  if(opt_constr(p->ed.constraints,d2)==0)
    fD2=p->ld[d2].D_z0;
  else
    fD2=p->ld[d2].D_zn, off2 = p->size[d2+1]-p->size[d2];
  off = off1+off2;
  for(d=0;d<3;++d) {
    unsigned i,j;
    const real *in = p->elx[d]+off;
    for(i=0;i<ne;++i,in+=se) {
      const real *in1 = in - off1, *in2 = in - off2;
      real *fd1 = &p->ed.fd1[d][i], *fd2 = &p->ed.fd2[d][i];
      p->ed.x[d][i] = *in;
      *fd1 = *fd2 = 0;
      for(j=0;j<n1;++j,in1+=s1) *fd1 += *in1 * fD1[j];
      for(j=0;j<n2;++j,in2+=s2) *fd2 += *in2 * fD2[j];
    }
  }
}

static void opt_edge_set_3(opt_data_3 *p, const real r[3], unsigned constr)
{
  if(p->ed.constraints!=constr) {
    p->ed.constraints=constr;
    p->ed.de = opt_constr_not[constr];
    p->ed.d1 = opt_other1_3[p->ed.de];
    p->ed.d2 = opt_other2_3[p->ed.de];
    opt_edge_proj_3(p);
  }
  lagrange_1(&p->ld[p->ed.de],r[p->ed.de]);
}

static void opt_edge_intp_3(opt_data_3 *p)
{
  unsigned d;
  const unsigned de = p->ed.de, d1 = p->ed.d1, d2 = p->ed.d2,
                 n = p->ld[de].n;
  const real *J = p->ld[de].J, *D = p->ld[de].D;
  
  for(d=0;d<3;++d) {
    p->x[d] = tensor_ig1(J,D,n, p->ed.x[d], &p->jac[d*3+de]);
    p->jac[d*3+d1] = tensor_i1(J,n, p->ed.fd1[d]);
    p->jac[d*3+d2] = tensor_i1(J,n, p->ed.fd2[d]);
  }
}

static void opt_edge_set_intp_3(opt_data_3 *p, const real r[3], unsigned constr)
{
  opt_edge_set_3(p,r,constr);
  opt_edge_intp_3(p);
}

static void opt_edge_hess_3(opt_data_3 *p, real hess[3])
{
  unsigned d;
  const unsigned de = p->ed.de, n = p->ld[de].n;
  const real *D2 = p->ld[de].D2;
  lagrange_2u(&p->ld[de]);
  for(d=0;d<3;++d) hess[d] = tensor_i1(D2,n, p->ed.x[d]);
}

static void opt_point_proj_3(opt_data_3 *p)
{
  unsigned off[3], offt, d, c[3];
  const real *fD[3];
  opt_constr_unpack_3(p->pd.constraints,c);
  for(d=0;d<3;++d)
    if(c[d]==0)
      fD[d]=p->ld[d].D_z0,off[d]=0;
    else
      fD[d]=p->ld[d].D_zn,off[d]=p->size[d+1]-p->size[d];
  offt = off[0]+off[1]+off[2];
  for(d=0;d<9;++d) p->pd.jac[d]=0;
  for(d=0;d<3;++d) {
    unsigned i,j;
    p->pd.x[d] = p->elx[d][offt];
    for(i=0;i<3;++i) {
      const real *in = p->elx[d]+offt-off[i];
      for(j=0;j<p->ld[i].n;++j,in+=p->size[i])
        p->pd.jac[d*3+i] += *in * fD[i][j];
    }
  }
}

static void opt_point_set_3(opt_data_3 *p, unsigned constr)
{
  if(p->pd.constraints!=constr) {
    p->pd.constraints=constr;
    opt_point_proj_3(p);
  }
}

static void opt_point_intp_3(opt_data_3 *p)
{
  memcpy(p->x,p->pd.x,3*sizeof(real));
  memcpy(p->jac,p->pd.jac,9*sizeof(real));
}

static void opt_point_set_intp_3(opt_data_3 *p, unsigned constr)
{
  opt_point_set_3(p,constr);
  opt_point_intp_3(p);
}

#define DIAGNOSTICS 0

static double opt_findpt_3(opt_data_3 *p, const real *const elx[3],
                           const real xstar[3], real r[3], unsigned *constr)
{
  real dr[3], resid[3], steep[3];

  unsigned c=*constr,ac,d,cc[3],step=0;
  
  p->elx[0]=elx[0], p->elx[1]=elx[1], p->elx[2]=elx[2];
  
  p->fd.constraints = opt_no_constraints_3;
  p->ed.constraints = opt_no_constraints_3;
  p->pd.constraints = opt_no_constraints_3;
  
# if DIAGNOSTICS
  printf("opt_findpt: xstar = %g, %g, %g\n", xstar[0], xstar[1], xstar[2]);
# endif
  
  do {
    ++step;
    if(step==50) fail("%s: opt_findpt_3 did not converge\n",__FILE__);
#   if DIAGNOSTICS
    printf("  iteration %u\n", step);
    printf("    %d constraint(s) active\n", (int)opt_constr_num_3[c]);
#   endif
    /* update face/edge/point data if necessary,
       and evaluate x(r) as well as the jacobian */
    switch(opt_constr_num_3[c]) {
      case 0: opt_vol_set_intp_3(p,r); break;
      case 1: opt_face_set_intp_3(p,r,c); break;
      case 2: opt_edge_set_intp_3(p,r,c); break;
      case 3: opt_point_set_intp_3(p,c); break;
    }
#   if DIAGNOSTICS
    printf("    r = %g, %g, %g\n", r[0], r[1], r[2]);
    printf("    x = %g, %g, %g\n", p->x[0], p->x[1], p->x[2]);
#   endif
    /* compute residual */
    for(d=0;d<3;++d) resid[d]=xstar[d]-p->x[d];
#   if DIAGNOSTICS
    printf("    resid = %g, %g, %g\n", resid[0], resid[1], resid[2]);
    printf("    2-norm = %g\n", r2norm_3(resid[0],resid[1],resid[2]));
#   endif
    /* check constraints against steepest descent direction */
    ac = c;
    if(opt_constr_num_3[c]) {
      opt_constr_unpack_3(c,cc);
      mat_app_3c(steep,p->jac,resid); /* steepest descent = J^T r */
#     if DIAGNOSTICS
      printf("    steepest descent = %g, %g, %g\n", steep[0],steep[1],steep[2]);
#     endif
      for(d=0;d<3;++d)
        if((cc[d]==0 && steep[d]>0) || (cc[d]==2 && steep[d]<0)) cc[d]=1;
      ac = opt_constr_pack_3(cc);
    }
    /* update face/edge/point data if necessary */
    if(ac!=c) {
      c=ac;
#     if DIAGNOSTICS
      printf("    relaxed to %d constraints\n", (int)opt_constr_num_3[c]);
#     endif
      switch(opt_constr_num_3[c]) {
        case 1: opt_face_set_3(p,r,c); break;
        case 2: opt_edge_set_3(p,r,c); break;
        case 3: opt_point_set_3(p,c); break;
      }
    }
    /* compute Newton step */
    switch(opt_constr_num_3[c]) {
      case 0: tinyla_solve_3(dr,p->jac,resid); break;
      case 1: {
        const unsigned dn = p->fd.dn, d1 = p->fd.d1, d2 = p->fd.d2;
        real A[4], H[9];
        const real *J = p->jac;
        opt_face_hess_3(p,H);
        A[0] = J[d1]*J[d1] + J[3+d1]*J[3+d1] + J[6+d1]*J[6+d1];
        A[1] = J[d2]*J[d2] + J[3+d2]*J[3+d2] + J[6+d2]*J[6+d2];
        A[2] = J[d1]*J[d2] + J[3+d1]*J[3+d2] + J[6+d1]*J[6+d2];
        A[0] -= resid[0]*H[0] + resid[1]*H[3] + resid[2]*H[6];
        A[1] -= resid[0]*H[1] + resid[1]*H[4] + resid[2]*H[7];
        A[2] -= resid[0]*H[2] + resid[1]*H[5] + resid[2]*H[8];
        tinyla_solve_sym_2(&dr[d1],&dr[d2],A,steep[d1],steep[d2]);
        dr[dn]=0;
      } break;
      case 2: {
        const unsigned de = p->ed.de, d1 = p->ed.d1, d2 = p->ed.d2;
        real fac, H[3];
        const real *J = p->jac+de;
        opt_edge_hess_3(p,H);
        fac = J[0]*J[0]+J[3]*J[3]+J[6]*J[6]
             -(resid[0]*H[0]+resid[1]*H[1]+resid[2]*H[2]);
        dr[de] = steep[de] / fac;
        dr[d1] = 0, dr[d2] = 0;
      } break;
      case 3:
        dr[0] = dr[1] = dr[2] = 0;
        break;
    }
#   if DIAGNOSTICS
    printf("    dr = %g, %g, %g\n", dr[0], dr[1], dr[2]);
#   endif
    /* project new iteration onto [-1,1]^3 */
    opt_constr_unpack_3(c,cc);
    for(d=0;d<3;++d) {
      if(cc[d]!=1) continue;
      r[d] += dr[d];
      if(r[d] <= -1)
        dr[d] -= r[d]+1, r[d] = -1, cc[d]=0;
      else if(r[d] >= 1)
        dr[d] -= r[d]-1, r[d] = 1, cc[d]=2;
    }
    c = opt_constr_pack_3(cc);
  } while(r1norm_3(dr[0],dr[1],dr[2]) > 30*EPS);
  *constr = c;
# if 0
  printf("opt_findpt_3 converged in %u iterations\n", step);
# endif
  return r2norm_3(resid[0],resid[1],resid[2]);
}

#undef DIAGNOSTICS

/*--------------------------------------------------------------------------
   
   2 - D
     
  --------------------------------------------------------------------------*/

typedef struct {
  unsigned constraints;
  unsigned de, d1;
  real *x[2], *fd1[2];
} opt_edge_data_2;

typedef struct {
  unsigned constraints;
  real x[2], jac[4];
} opt_point_data_2;

typedef struct {
  lagrange_data *ld;
  unsigned size[3];
  const real *elx[2];
  opt_edge_data_2 ed;
  opt_point_data_2 pd;
  real *work;
  real x[2], jac[4];
} opt_data_2;

static void opt_alloc_2(opt_data_2 *p, lagrange_data *ld)
{
  const unsigned nr = ld[0].n, ns = ld[1].n,
                 ne = umax_2(nr,ns),
                 nw = 2*ns;
  p->size[0] = 1;
  p->size[1] = nr;
  p->size[2] = nr*ns;
  p->ld = ld;
  p->work = tmalloc(real, 4*ne + nw);
  p->ed.x[0] = p->work + nw; 
  p->ed.x[1] = p->ed.x[0] + ne; 
  p->ed.fd1[0] = p->ed.x[1] + ne; 
  p->ed.fd1[1] = p->ed.fd1[0] + ne; 
}

static void opt_free_2(opt_data_2 *p)
{
  free(p->work);
}

static void opt_area_set_2(opt_data_2 *p, const real r[2])
{
  lagrange_1(&p->ld[0],r[0]);
  lagrange_1(&p->ld[1],r[1]);
}

/* work holds 2*ns reals */
static void opt_area_intp_2(opt_data_2 *p)
{
  unsigned d;
  const lagrange_data *ld = p->ld;
  
  for(d=0;d<2;++d) 
    p->x[d] = tensor_ig2(ld[0].J,ld[0].D,ld[0].n,
                         ld[1].J,ld[1].D,ld[1].n,
                         p->elx[d], &p->jac[d*2], p->work);
}

static void opt_area_set_intp_2(opt_data_2 *p, const real r[2])
{
  opt_area_set_2(p,r);
  opt_area_intp_2(p);
}

static void opt_edge_proj_2(opt_data_2 *p)
{
  unsigned d, off=0;
  const unsigned de = p->ed.de, d1 = p->ed.d1,
                 se=p->size[de], s1=p->size[d1],
                 ne=p->ld[de].n, n1=p->ld[d1].n;
  const real *fD1;
  if(opt_constr(p->ed.constraints,d1)==0)
    fD1=p->ld[d1].D_z0;
  else
    fD1=p->ld[d1].D_zn, off = p->size[d1+1]-p->size[d1];
  for(d=0;d<2;++d) {
    unsigned i,j;
    const real *in = p->elx[d]+off;
    for(i=0;i<ne;++i,in+=se) {
      const real *in1 = in - off;
      real *fd1 = &p->ed.fd1[d][i];
      p->ed.x[d][i] = *in;
      *fd1 = 0;
      for(j=0;j<n1;++j,in1+=s1) *fd1 += *in1 * fD1[j];
    }
  }
}

static void opt_edge_set_2(opt_data_2 *p, const real r[2], unsigned constr)
{
  if(p->ed.constraints!=constr) {
    p->ed.constraints=constr;
    p->ed.de = opt_constr_not[constr];
    p->ed.d1 = 1 - p->ed.de;
    opt_edge_proj_2(p);
  }
  lagrange_1(&p->ld[p->ed.de],r[p->ed.de]);
}

static void opt_edge_intp_2(opt_data_2 *p)
{
  unsigned d;
  const unsigned de = p->ed.de, d1 = p->ed.d1, n = p->ld[de].n;
  const real *J = p->ld[de].J, *D = p->ld[de].D;
  for(d=0;d<2;++d) {
    p->x[d] = tensor_ig1(J,D,n, p->ed.x[d], &p->jac[d*2+de]);
    p->jac[d*2+d1] = tensor_i1(J,n, p->ed.fd1[d]);
  }
}

static void opt_edge_set_intp_2(opt_data_2 *p, const real r[2], unsigned constr)
{
  opt_edge_set_2(p,r,constr);
  opt_edge_intp_2(p);
}

static void opt_edge_hess_2(opt_data_2 *p, real hess[2])
{
  unsigned d;
  const unsigned de = p->ed.de, n = p->ld[de].n;
  const real *D2 = p->ld[de].D2;
  lagrange_2u(&p->ld[de]);
  for(d=0;d<2;++d) hess[d] = tensor_i1(D2,n, p->ed.x[d]);
}

static void opt_point_proj_2(opt_data_2 *p)
{
  unsigned off[2], offt, d, c[2];
  const real *fD[2];
  opt_constr_unpack_2(p->pd.constraints,c);
  for(d=0;d<2;++d)
    if(c[d]==0)
      fD[d]=p->ld[d].D_z0,off[d]=0;
    else
      fD[d]=p->ld[d].D_zn,off[d]=p->size[d+1]-p->size[d];
  offt = off[0]+off[1];
  for(d=0;d<4;++d) p->pd.jac[d]=0;
  for(d=0;d<2;++d) {
    unsigned i,j;
    p->pd.x[d] = p->elx[d][offt];
    for(i=0;i<2;++i) {
      const real *in = p->elx[d]+offt-off[i];
      for(j=0;j<p->ld[i].n;++j,in+=p->size[i])
        p->pd.jac[d*2+i] += *in * fD[i][j];
    }
  }
}

static void opt_point_set_2(opt_data_2 *p, unsigned constr)
{
  if(p->pd.constraints!=constr) {
    p->pd.constraints=constr;
    opt_point_proj_2(p);
  }
}

static void opt_point_intp_2(opt_data_2 *p)
{
  memcpy(p->x,p->pd.x,2*sizeof(real));
  memcpy(p->jac,p->pd.jac,4*sizeof(real));
}

static void opt_point_set_intp_2(opt_data_2 *p, unsigned constr)
{
  opt_point_set_2(p,constr);
  opt_point_intp_2(p);
}

#define DIAGNOSTICS 0

static double opt_findpt_2(opt_data_2 *p, const real *const elx[2],
                           const real xstar[2], real r[2], unsigned *constr)
{
  real dr[2], resid[2], steep[2];

  unsigned c=*constr,ac,d,cc[2],step=0;
  
  p->elx[0]=elx[0], p->elx[1]=elx[1];
  
  p->ed.constraints = opt_no_constraints_2;
  p->pd.constraints = opt_no_constraints_2;
  
# if DIAGNOSTICS
  printf("opt_findpt: xstar = %g, %g\n", xstar[0], xstar[1]);
# endif
  
  do {
    ++step;
    if(step==150) fail("%s: opt_findpt_2 did not converge\n",__FILE__);
#   if DIAGNOSTICS
    printf("  iteration %u\n", step);
    printf("    %d constraint(s) active\n", (int)opt_constr_num_2[c]);
#   endif
    /* update face/edge/point data if necessary,
       and evaluate x(r) as well as the jacobian */
    switch(opt_constr_num_2[c]) {
      case 0: opt_area_set_intp_2(p,r); break;
      case 1: opt_edge_set_intp_2(p,r,c); break;
      case 2: opt_point_set_intp_2(p,c); break;
    }
#   if DIAGNOSTICS
    printf("    r = %g, %g\n", r[0], r[1]);
    printf("    x = %g, %g\n", p->x[0], p->x[1]);
#   endif
    /* compute residual */
    for(d=0;d<2;++d) resid[d]=xstar[d]-p->x[d];
#   if DIAGNOSTICS
    printf("    resid = %g, %g\n", resid[0], resid[1]);
    printf("    2-norm = %g\n", r2norm_2(resid[0],resid[1]));
#   endif
    /* check constraints against steepest descent direction */
    ac = c;
    if(opt_constr_num_2[c]) {
      opt_constr_unpack_2(c,cc);
      mat_app_2c(steep,p->jac,resid); /* steepest descent = J^T r */
#     if DIAGNOSTICS
      printf("    steepest descent = %g, %g\n", steep[0], steep[1]);
#     endif
      for(d=0;d<2;++d)
        if((cc[d]==0 && steep[d]>0) || (cc[d]==2 && steep[d]<0)) cc[d]=1;
      ac = opt_constr_pack_2(cc);
    }
    /* update face/edge/point data if necessary */
    if(ac!=c) {
      c=ac;
#     if DIAGNOSTICS
      printf("    relaxed to %d constraints\n", (int)opt_constr_num_2[c]);
#     endif
      switch(opt_constr_num_2[c]) {
        case 1: opt_edge_set_2(p,r,c); break;
        case 2: opt_point_set_2(p,c); break;
      }
    }
    /* compute Newton step */
    switch(opt_constr_num_2[c]) {
      case 0: tinyla_solve_2(dr,p->jac,resid); break;
      case 1: {
        const unsigned de = p->ed.de, d1 = p->ed.d1;
        real fac, H[2];
        const real *J = p->jac+de;
        opt_edge_hess_2(p,H);
        fac = J[0]*J[0]+J[2]*J[2]-(resid[0]*H[0]+resid[1]*H[1]);
        dr[de] = steep[de] / fac;
        dr[d1] = 0;
      } break;
      case 2:
        dr[0] = dr[1] = 0;
        break;
    }
#   if DIAGNOSTICS
    printf("    dr = %g, %g\n", dr[0], dr[1]);
#   endif
    /* project new iteration onto [-1,1]^2 */
    opt_constr_unpack_2(c,cc);
    for(d=0;d<2;++d) {
      if(cc[d]!=1) continue;
      r[d] += dr[d];
      if(r[d] <= -1)
        dr[d] -= r[d]+1, r[d] = -1, cc[d]=0;
      else if(r[d] >= 1)
        dr[d] -= r[d]-1, r[d] = 1, cc[d]=2;
    }
    c = opt_constr_pack_2(cc);
  } while(r1norm_2(dr[0],dr[1]) > 30*EPS);
  *constr = c;
  return r2norm_2(resid[0],resid[1]);
}

#undef DIAGNOSTICS

/*--------------------------------------------------------------------------
   Point Finding  (interface/top-level)
          
   Initializing the data:
     unsigned nel;        // number of elements
     const unsigned n[3]; // number of nodes in r, s, t directions
     const real *xm[3];   // n[0]*n[1]*n[2]*nel x,y,z coordinates
     real tol = 0.01;     // how far point is allowed to be outside element
                          //   relative to element size
     unsigned max_size = n[0]*n[1]*n[2]*nel; // maximum size of hash table
                          
     findpt_data_3 *data = findpt_setup_3(xm,n,nel,max_size,tol);
     
   Using the data:
     real x[3] = { ... };   // point to find
     int el; // element number
     real r[3]; // parametric coordinates
     int guess = 0; // do we use (el,r,s,t) as an initial guess?
     int code; // 0 : normal, -1 : outside all elements,
               // 1 : border, or outside but within tolerance
     real dist; // distance in xyz space from returned (el,r,s,t) to given
                // (x,y,z)

     code = findpt_3(data, x, guess, &el, r, &dist);
   
   When done:
     findpt_free_3(&data);
     
  --------------------------------------------------------------------------*/

typedef struct {
  uint el;
  real r[3];
  real dist;
} findpt_listel;

/* heap sort on A[0:n-1] with key A[i]->dist
   precondition: n!=0 */
static void findpt_list_sort(findpt_listel **A, unsigned n)
{
  unsigned i;
  --A; /* make A have a base index of 1 */
  /* build heap */
  for(i=2;i<=n;++i) {
    findpt_listel *item = A[i];
    unsigned hole = i, parent = hole>>1;
    if(A[parent]->dist >= item->dist) continue;
    do {
      A[hole] = A[parent];
      hole = parent;
      parent>>=1;
    } while(parent && A[parent]->dist < item->dist);
    A[hole] = item;
  }
  /* extract */
  for(i=n-1;i;--i) {
    findpt_listel *item = A[i+1];
    unsigned hole = 1;
    A[i+1] = A[1];
    for(;;) {
      unsigned ch = hole<<1, r = ch+1;
      if(r<=i && A[ch]->dist < A[r]->dist) ch=r;
      if(ch>i || item->dist >= A[ch]->dist) break;
      A[hole]=A[ch];
      hole=ch;
    }
    A[hole] = item;
  }
}

typedef struct {
  const real *xw[2];   /* geometry data */
  real *z[2];          /* lobatto nodes */
  lagrange_data ld[2]; /* interpolation, derivative weights & data */
  unsigned nptel;      /* nodes per element */
  hash_data_2 *hash;   /* geometric hashing data */
  findpt_listel *list, **sorted, **end; /* pre-allocated list of elements to
                                           check (found by hashing), and
                                           pre-allocated list of pointers into
                                           the first list for sorting */
  opt_data_2 *od; /* data for the optimization algorithm */
  real *od_work;
} findpt_data_2;

typedef struct {
  const real *xw[3];   /* geometry data */
  real *z[3];          /* lobatto nodes */
  lagrange_data ld[3]; /* interpolation, derivative weights & data */
  unsigned nptel;      /* nodes per element */
  hash_data_3 *hash;   /* geometric hashing data */
  findpt_listel *list, **sorted, **end; /* pre-allocated list of elements to
                                           check (found by hashing), and
                                           pre-allocated list of pointers into
                                           the first list for sorting */
  opt_data_3 *od; /* data for the optimization algorithm */
  real *od_work;
} findpt_data_3;

findpt_data_2 *findpt_setup_2(
          const real *const xw[2], const unsigned n[2], uint nel,
          uint max_hash_size, real bbox_tol)
{
  unsigned d;
  findpt_data_2 *p = tmalloc(findpt_data_2,1);

  p->hash = tmalloc(hash_data_2,1);
  p->od = tmalloc(opt_data_2,1);

  for(d=0;d<2;++d) p->xw[d]=xw[d];
  p->nptel = n[0]*n[1];

  hash_build_2(p->hash,xw,n,nel,max_hash_size,bbox_tol);

  for(d=0;d<2;++d) {
    p->z[d] = tmalloc(real,n[d]);
    lobatto_nodes(p->z[d],n[d]);
    lagrange_setup(&p->ld[d],p->z[d],n[d]);
  }

  p->list   = tmalloc(findpt_listel , p->hash->max);
  p->sorted = tmalloc(findpt_listel*, p->hash->max);
  
  opt_alloc_2(p->od,p->ld);
  p->od_work = p->od->work;
  
  return p;
}

findpt_data_3 *findpt_setup_3(
          const real *const xw[3], const unsigned n[3], uint nel,
          uint max_hash_size, real bbox_tol)
{
  unsigned d;
  findpt_data_3 *p = tmalloc(findpt_data_3,1);

  p->hash = tmalloc(hash_data_3,1);
  p->od = tmalloc(opt_data_3,1);

  for(d=0;d<3;++d) p->xw[d]=xw[d];
  p->nptel = n[0]*n[1]*n[2];

  hash_build_3(p->hash,xw,n,nel,max_hash_size,bbox_tol);

  for(d=0;d<3;++d) {
    p->z[d] = tmalloc(real,n[d]);
    lobatto_nodes(p->z[d],n[d]);
    lagrange_setup(&p->ld[d],p->z[d],n[d]);
  }

  p->list   = tmalloc(findpt_listel , p->hash->max);
  p->sorted = tmalloc(findpt_listel*, p->hash->max);
  
  opt_alloc_3(p->od,p->ld);
  p->od_work = p->od->work;
  
  return p;
}

void findpt_free_2(findpt_data_2 *p)
{
  unsigned d;
  opt_free_2(p->od); free(p->od);
  hash_free_2(p->hash); free(p->hash);
  free(p->list);
  free(p->sorted);
  for(d=0;d<2;++d) free(p->z[d]);
  free(p);
}

void findpt_free_3(findpt_data_3 *p)
{
  unsigned d;
  opt_free_3(p->od); free(p->od);
  hash_free_3(p->hash); free(p->hash);
  free(p->list);
  free(p->sorted);
  for(d=0;d<3;++d) free(p->z[d]);
  free(p);
}

const real *findpt_allbnd_2(const findpt_data_2 *p)
{
  return p->hash->bnd;
}

const real *findpt_allbnd_3(const findpt_data_3 *p)
{
  return p->hash->bnd;
}

static void findpt_hash_2(findpt_data_2 *p, const real x[2])
{
  findpt_listel *list = p->list, **sorted = p->sorted;
  const uint hi = hash_index_2(p->hash, x);
  const uint *offset = p->hash->offset;
  uint i; const uint b = offset[hi], e = offset[hi+1];
  for(i=b;i!=e;++i) {
    const uint el = offset[i];
    real *r = &list->r[0];
    const obbox_2 *obb = &p->hash->obb[el];
    if(obbox_axis_test_2(obb,x)) continue;
    if(obbox_test_2(obb,x,r)) continue;
    list->el = el;
    list->dist = r1norm_2(r[0],r[1]);
    *sorted++ = list++;
  }
  p->end = sorted;
  if(p->end!=p->sorted)
    findpt_list_sort(p->sorted,p->end - p->sorted);
}

static void findpt_hash_3(findpt_data_3 *p, const real x[3])
{
  findpt_listel *list = p->list, **sorted = p->sorted;
  const uint hi = hash_index_3(p->hash, x);
  const uint *offset = p->hash->offset;
  uint i; const uint b = offset[hi], e = offset[hi+1];
  for(i=b;i!=e;++i) {
    const uint el = offset[i];
    real *r = &list->r[0];
    const obbox_3 *obb = &p->hash->obb[el];
    if(obbox_axis_test_3(obb,x)) continue;
    if(obbox_test_3(obb,x,r)) continue;
    list->el = el;
    list->dist = r1norm_3(r[0],r[1],r[2]);
    *sorted++ = list++;
  }
  p->end = sorted;
  if(p->end!=p->sorted)
    findpt_list_sort(p->sorted,p->end - p->sorted);
}

static int findpt_guess_2(findpt_data_2 *p, const real x[2],
                          uint el, real r[2], real *dist)
{
  const uint index = p->nptel*el;
  const real *elx[2] = {p->xw[0]+index,p->xw[1]+index};
  real g[2];
  unsigned c = opt_no_constraints_2;
  const obbox_2 *obb = &p->hash->obb[el];
  if(obbox_axis_test_2(obb,x) || obbox_test_2(obb,x,g)) return 0;
  *dist = opt_findpt_2(p->od,elx,x,r,&c);
  return c==opt_no_constraints_2;
}

static int findpt_guess_3(findpt_data_3 *p, const real x[3],
                          uint el, real r[3], real *dist)
{
  const uint index = p->nptel*el;
  const real *elx[3] = {p->xw[0]+index,p->xw[1]+index,p->xw[2]+index};
  real g[3];
  unsigned c = opt_no_constraints_3;
  const obbox_3 *obb = &p->hash->obb[el];
  if(obbox_axis_test_3(obb,x) || obbox_test_3(obb,x,g)) return 0;
  *dist = opt_findpt_3(p->od,elx,x,r,&c);
  return c==opt_no_constraints_3;
}

#define DIAGNOSTICS 0

static int findpt_pass_2(findpt_data_2 *p, const real x[2],
                         uint *el, real r[2], real *dist_min)
{
  findpt_listel **qq = p->sorted;
  const real *bnd;
  do {
    findpt_listel *q = *qq;
    const uint index = p->nptel*q->el;
    const real *elx[2] = {p->xw[0]+index,p->xw[1]+index};
    unsigned c = opt_no_constraints_2;
    const real dist = opt_findpt_2(p->od,elx,x,q->r,&c);
    if(qq==p->sorted || dist<*dist_min || c==opt_no_constraints_2) {
      *dist_min = dist;
      *el = q->el;
      memcpy(r, q->r, 2*sizeof(real));
      if(c==opt_no_constraints_2) return 0;
    }
  } while(++qq != p->end);
  bnd = p->hash->obb[*el].axis_bnd;
  return *dist_min>r2norm_2(bnd[1]-bnd[0],bnd[3]-bnd[2]) ? -1 : 1;
}

static int findpt_pass_3(findpt_data_3 *p, const real x[3],
                         uint *el, real r[3], real *dist_min)
{
  findpt_listel **qq = p->sorted;
  const real *bnd;
  do {
    findpt_listel *q = *qq;
    const uint index = p->nptel*q->el;
    const real *elx[3] = {p->xw[0]+index,p->xw[1]+index,p->xw[2]+index};
    unsigned c = opt_no_constraints_3;
    const real dist = opt_findpt_3(p->od,elx,x,q->r,&c);
    if(qq==p->sorted || dist<*dist_min || c==opt_no_constraints_3) {
      *dist_min = dist;
      *el = q->el;
      memcpy(r, q->r, 3*sizeof(real));
      if(c==opt_no_constraints_3) {
#       if DIAGNOSTICS
        printf("point found in element #%d\n", qq-p->sorted);
#       endif
        return 0;
      }
    }
  } while(++qq != p->end);
  bnd = p->hash->obb[*el].axis_bnd;
  return *dist_min>r2norm_3(bnd[1]-bnd[0],bnd[3]-bnd[2],bnd[5]-bnd[4]) ? -1 : 1;
}

int findpt_2(findpt_data_2 *p, const real x[2], int guess,
             uint *el, real r[2], real *dist)
{
  if(guess && findpt_guess_2(p,x,*el,r,dist)) return 0;
  findpt_hash_2(p,x);
  if(p->sorted==p->end) return -1;
  return findpt_pass_2(p,x,el,r,dist);
}

int findpt_3(findpt_data_3 *p, const real x[3], int guess,
             uint *el, real r[3], real *dist)
{
  if(guess && findpt_guess_3(p,x,*el,r,dist)) return 0;
  findpt_hash_3(p,x);
# if DIAGNOSTICS
  printf("hashing leaves %d elements to consider\n",p->end-p->sorted);
# endif
  if(p->sorted==p->end) return -1;
  return findpt_pass_3(p,x,el,r,dist);
}

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

