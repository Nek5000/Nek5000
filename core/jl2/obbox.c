#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "tensor.h"
#include "poly.h"
#include "lob_bnd.h"

#define obbox_2_calc  PREFIXED_NAME(obbox_2_calc)
#define obbox_3_calc  PREFIXED_NAME(obbox_3_calc)

typedef struct { double c0[2], A[4];
                 struct dbl_range x,y; } obbox_2;

typedef struct { double c0[3], A[9];
                 struct dbl_range x,y,z; } obbox_3;


static void copy_strided(double *out, const double *in,
                         unsigned g, unsigned s, unsigned n)
{
  if(g==1) for(;n;--n,in+=s) *out++ = *in;
  else {
    s *= g;
    for(;n;--n,in+=s) memcpy(out,in,g*sizeof(double)), out+=g;
  }
}

static void mat_inv_2(double inv[4], const double A[4])
{
  const double idet = 1/(A[0]*A[3]-A[1]*A[2]);
  inv[0] =   idet*A[3];
  inv[1] = -(idet*A[1]);
  inv[2] = -(idet*A[2]);
  inv[3] =   idet*A[0];
}

static void mat_inv_3(double inv[9], const double A[9])
{
  const double a = A[4]*A[8]-A[5]*A[7],
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

static struct dbl_range dbl_range_merge(struct dbl_range a, struct dbl_range b)
{
  struct dbl_range m = { b.min<a.min?b.min:a.min,
                         a.max>b.max?a.max:b.max };
  return m;
}

static struct dbl_range dbl_range_expand(struct dbl_range b, double tol)
{
  double a = (b.min+b.max)/2, l = (b.max-b.min)*(1+tol)/2;
  struct dbl_range m = { a-l,a+l };
  return m;
}

static void bbox_2_tfm(double *out, const double x0[2], const double Ji[4],
                       const double *x, const double *y, unsigned n)
{
  unsigned i;
  for(i=0;i<n;++i) {
    const double dx = x[i]-x0[0], dy = y[i]-x0[1];
    out[  i] = Ji[0]*dx + Ji[1]*dy;
    out[n+i] = Ji[2]*dx + Ji[3]*dy;
  }
}

static void bbox_3_tfm(double *out, const double x0[3], const double Ji[9],
                       const double *x, const double *y, const double *z,
                       unsigned n)
{
  unsigned i;
  for(i=0;i<n;++i) {
    const double dx = x[i]-x0[0], dy = y[i]-x0[1], dz = z[i]-x0[2];
    out[    i] = Ji[0]*dx + Ji[1]*dy + Ji[2]*dz;
    out[  n+i] = Ji[3]*dx + Ji[4]*dy + Ji[5]*dz;
    out[2*n+i] = Ji[6]*dx + Ji[7]*dy + Ji[8]*dz;
  }
}

#if 0

/* positive when possibly inside */
double obbox_2_axis_test(const obbox_2 *b, double x, double y)
{
  double bx = (x-b->x.min)*(b->x.max-x);
  return bx<0 ? bx : (y-b->y.min)*(b->y.max-y);
}

/* positive when possibly inside */
double obbox_2_test(const obbox_2 *b, double x, double y)
{
  double bxy = obbox_2_axis_test(b,x,y);
  if(bxy<0) return bxy; else {
    double dx = x-b->c0[0], dy = y-b->c0[1];
    double r = b->A[0]*dx + b->A[1]*dy,
           s = b->A[2]*dx + b->A[3]*dy;
    double br = (r+1)*(1-r);
    return br<0 ? br : (s+1)*(1-s);
  }
}

#endif

#define DO_MAX(a,b) do { unsigned temp = b; if(temp>a) a=temp; } while(0)

void obbox_2_calc(obbox_2 *out,
  const double *x, const double *y,
  unsigned nr, unsigned ns, uint n,
  unsigned mr, unsigned ms, double tol)
{
  unsigned nrs = nr*ns;
  double *data;
  unsigned lbsize[2] = {lob_bnd_size(nr,mr), lob_bnd_size(ns,ms)};
  unsigned wsize = 4*ns+2*ms;
  DO_MAX(wsize,2*nr+2*mr);
  DO_MAX(wsize,gll_lag_size(nr));
  DO_MAX(wsize,gll_lag_size(ns));
  data = tmalloc(double, 2*(nr+ns)+lbsize[0]+lbsize[1]+wsize);

  {
    double *const I0[2] = { data, data+2*nr };
    double *const lob_bnd_data_r = data+2*(nr+ns),
           *const lob_bnd_data_s = data+2*(nr+ns)+lbsize[0];
    double *const work = data+2*(nr+ns)+lbsize[0]+lbsize[1];

    #define SETUP_DIR(i,r) do { \
      lagrange_fun *const lag = gll_lag_setup(work, n##r); \
      lag(I0[i], work,n##r,1, 0); \
      lob_bnd_setup(lob_bnd_data_##r, n##r,m##r); \
    } while(0)
    
    SETUP_DIR(0,r); SETUP_DIR(1,s);
    
    #undef SETUP_DIR
    
    for(;n;--n,x+=nrs,y+=nrs,++out) {
      double x0[2], J[4], Ji[4];
      struct dbl_range ab[2], tb[2];
  
      /* double work[2*nr] */
      x0[0] = tensor_ig2(J  , I0[0],nr, I0[1],ns, x, work);
      x0[1] = tensor_ig2(J+2, I0[0],nr, I0[1],ns, y, work);
      mat_inv_2(Ji, J);

      /* double work[2*m##r] */
      #define DO_BOUND(bnd,merge,r,x,work) do { \
        struct dbl_range b = \
        lob_bnd_1(lob_bnd_data_##r,n##r,m##r, x, work); \
        if(merge) bnd=dbl_range_merge(bnd,b); else bnd=b; \
      } while(0)

      /* double work[2*n##r + 2*m##r] */
      #define DO_EDGE(merge,r,x,y,work) do { \
        DO_BOUND(ab[0],merge,r,x,work); \
        DO_BOUND(ab[1],merge,r,y,work); \
        bbox_2_tfm(work, x0,Ji, x,y,n##r); \
        DO_BOUND(tb[0],merge,r,(work)     ,(work)+2*n##r); \
        DO_BOUND(tb[1],merge,r,(work)+n##r,(work)+2*n##r); \
      } while(0)

      DO_EDGE(0,r,x,y,work);
      DO_EDGE(1,r,&x[nrs-nr],&y[nrs-nr],work);

      /* double work[4*ns + 2*ms] */
      #define GET_EDGE(off) do { \
        copy_strided(work   , x+off,1,nr,ns); \
        copy_strided(work+ns, y+off,1,nr,ns); \
        DO_EDGE(1,s,work,work+ns,work+2*ns); \
      } while(0)
  
      GET_EDGE(0);
      GET_EDGE(nr-1);
  
      #undef GET_EDGE
      #undef DO_EDGE
      #undef DO_BOUND

      out->x = dbl_range_expand(ab[0],tol),
      out->y = dbl_range_expand(ab[1],tol);
  
      {
        double av[2] = {(tb[0].min+tb[0].max)/2,(tb[1].min+tb[1].max)/2};
        out->c0[0] = x0[0] + J[0]*av[0] + J[1]*av[1];
        out->c0[1] = x0[1] + J[2]*av[0] + J[3]*av[1];
      }
      {
        double di[2] = {2/((1+tol)*(tb[0].max-tb[0].min)),
                        2/((1+tol)*(tb[1].max-tb[1].min))};
        out->A[0]=di[0]*Ji[0], out->A[1]=di[0]*Ji[1];
        out->A[2]=di[1]*Ji[2], out->A[3]=di[1]*Ji[3];
      }

    }
  }
  
  free(data);  
}

void obbox_3_calc(obbox_3 *out,
  const double *x, const double *y, const double *z,
  unsigned nr, unsigned ns, unsigned nt, uint n,
  unsigned mr, unsigned ms, unsigned mt, double tol)
{
  unsigned nrs = nr*ns, nrst = nr*ns*nt;
  double *data;
  unsigned lbsize[3] = 
    {lob_bnd_size(nr,mr), lob_bnd_size(ns,ms), lob_bnd_size(nt,mt)};
  unsigned wsize = 3*nr*ns+2*mr*(ns+ms+1);
  DO_MAX(wsize,6*nr*nt+2*mr*(nt+mt+1));
  DO_MAX(wsize,6*ns*nt+2*ms*(nt+mt+1));
  DO_MAX(wsize,2*nr*ns+3*nr);
  DO_MAX(wsize,gll_lag_size(nr));
  DO_MAX(wsize,gll_lag_size(ns));
  DO_MAX(wsize,gll_lag_size(nt));
  data = tmalloc(double, 2*(nr+ns+nt)+lbsize[0]+lbsize[1]+lbsize[2]+wsize);

  {
    double *const I0[3] = { data, data+2*nr, data+2*(nr+ns) };
    double *const lob_bnd_data_r = data+2*(nr+ns+nt),
           *const lob_bnd_data_s = data+2*(nr+ns+nt)+lbsize[0],
           *const lob_bnd_data_t = data+2*(nr+ns+nt)+lbsize[0]+lbsize[1];
    double *const work = data+2*(nr+ns+nt)+lbsize[0]+lbsize[1]+lbsize[2];
    
    #define SETUP_DIR(i,r) do { \
      lagrange_fun *const lag = gll_lag_setup(work, n##r); \
      lag(I0[i], work,n##r,1, 0); \
      lob_bnd_setup(lob_bnd_data_##r, n##r,m##r); \
    } while(0)
    
    SETUP_DIR(0,r); SETUP_DIR(1,s); SETUP_DIR(2,t);
    
    #undef SETUP_DIR
    
    for(;n;--n,x+=nrst,y+=nrst,z+=nrst,++out) {
      double x0[3], J[9], Ji[9];
      struct dbl_range ab[3], tb[3];
  
      /* double work[2*nrs+3*nr] */
      #define EVAL_AT_0(d,x) \
        x0[d] = tensor_ig3(J+3*d, I0[0],nr, I0[1],ns, I0[2],nt, x, work)
      EVAL_AT_0(0,x); EVAL_AT_0(1,y); EVAL_AT_0(2,z);                          
      mat_inv_3(Ji, J);
      #undef EVAL_AT_0
 
      /* double work[2*m##r*(n##s+m##s+1)] */
      #define DO_BOUND(bnd,merge,r,s,x,work) do { \
        struct dbl_range b = \
        lob_bnd_2(lob_bnd_data_##r,n##r,m##r, \
                  lob_bnd_data_##s,n##s,m##s, x, work); \
        if(merge) bnd=dbl_range_merge(bnd,b); else bnd=b; \
      } while(0)

      /* double work[3*n##r*n##s+2*m##r*(n##s+m##s+1)] */
      #define DO_FACE(merge,r,s,x,y,z,work) do { \
        DO_BOUND(ab[0],merge,r,s,x,work); \
        DO_BOUND(ab[1],merge,r,s,y,work); \
        DO_BOUND(ab[2],merge,r,s,z,work); \
        bbox_3_tfm(work, x0,Ji, x,y,z,n##r*n##s); \
        DO_BOUND(tb[0],merge,r,s,(work)            ,(work)+3*n##r*n##s); \
        DO_BOUND(tb[1],merge,r,s,(work)+  n##r*n##s,(work)+3*n##r*n##s); \
        DO_BOUND(tb[2],merge,r,s,(work)+2*n##r*n##s,(work)+3*n##r*n##s); \
      } while(0)

      DO_FACE(0,r,s,x,y,z,work);
      DO_FACE(1,r,s,&x[nrst-nrs],&y[nrst-nrs],&z[nrst-nrs],work);

      /* double work[6*n##r*n##s+2*m##r*(n##s+m##s+1)] */
      #define GET_FACE(r,s,off,n1,n2,n3) do { \
        copy_strided(work            , x+off,n1,n2,n3); \
        copy_strided(work+  n##r*n##s, y+off,n1,n2,n3); \
        copy_strided(work+2*n##r*n##s, z+off,n1,n2,n3); \
        DO_FACE(1,r,s,work,work+n##r*n##s,work+2*n##r*n##s,work+3*n##r*n##s); \
      } while(0)
  
      GET_FACE(r,t,0     ,nr,ns,nt);
      GET_FACE(r,t,nrs-nr,nr,ns,nt);
      GET_FACE(s,t,0     , 1,nr,ns*nt);
      GET_FACE(s,t,nr-1  , 1,nr,ns*nt);
      
      #undef GET_FACE
      #undef DO_FACE
      #undef DO_BOUND

      out->x = dbl_range_expand(ab[0],tol),
      out->y = dbl_range_expand(ab[1],tol);
      out->z = dbl_range_expand(ab[2],tol);
  
      {
        double av[3] = {(tb[0].min+tb[0].max)/2,
                        (tb[1].min+tb[1].max)/2,
                        (tb[2].min+tb[2].max)/2};
        out->c0[0] = x0[0] + J[0]*av[0] + J[1]*av[1] + J[2]*av[2];
        out->c0[1] = x0[1] + J[3]*av[0] + J[4]*av[1] + J[5]*av[2];
        out->c0[2] = x0[2] + J[6]*av[0] + J[7]*av[1] + J[8]*av[2];
      }
      {
        double di[3] = {2/((1+tol)*(tb[0].max-tb[0].min)),
                        2/((1+tol)*(tb[1].max-tb[1].min)),
                        2/((1+tol)*(tb[2].max-tb[2].min))};
        out->A[0]=di[0]*Ji[0], out->A[1]=di[0]*Ji[1], out->A[2]=di[0]*Ji[2];
        out->A[3]=di[1]*Ji[3], out->A[4]=di[1]*Ji[4], out->A[5]=di[1]*Ji[5];
        out->A[6]=di[2]*Ji[6], out->A[7]=di[2]*Ji[7], out->A[8]=di[2]*Ji[8];
      }

    }
  }
  
  free(data);  
}



#if 0

/*--------------------------------------------------------------------------
  Borders
  --------------------------------------------------------------------------*/

static void get_edge(double *out, const double *u, unsigned nr, unsigned ns)
{
  for(;ns;--ns,u+=nr) *out++ = *u;
}

unsigned borders_2_size(unsigned nr, unsigned ns)
{
  return 2*ns;
}

void borders_2(const double *border[9], double *extra,
               const double *u, unsigned nr, unsigned ns)
{
  unsigned j, nrs=nr*ns;
  for(j=0;j<ns;++j) extra[j] = u+j*nr;
  for(j=0;j<ns;++j) extra[ns+j] = u+j*nr+(nr-1);
  border[0] = &u[0];
  border[1] = &u[0];
  border[2] = &u[nr-1];
  border[3] = &extra[0];
  border[4] = &u[0];
  border[5] = &extra[ns];
  border[6] = &u[nrs-nr];
  border[7] = &u[nrs-nr];
  border[8] = &u[nrs-1];
}

unsigned borders_3_size(unsigned nr, unsigned ns, unsigned nt)
{
  return 2*ns*nt + 2*nt*nr;
}

void borders_3(const double *border[27], double *extra,
               const double *u, unsigned nr, unsigned ns, unsigned nt)
{
  unsigned i,j,k, nrs=nr*ns, nst=ns*nt, ntr = nt*nr, nrst = nrs*nt;
  double *out=extra; const double *in;
  for(j=ns*nt,in=u       ;j;--j,in+=nr) *out++ = *in;
  for(j=ns*nt,in=u+(nr-1);j;--j,in+=nr) *out++ = *in;
  for(i=0;i<nr;++i)
    for(j=nt,in=u+i;j;--j,in+=nrs) *out++ = *in;
  for(i=0;i<nr;++i)
    for(j=nt,in=u+nsr-nr+i;j;--j,in+=nrs) *out++ = *in;
    
  border[0] = &u[0];
  border[1] = &u[0];
  border[2] = &u[nr-1];
  border[3] = &extra[0];
  border[4] = &u[0];
  border[5] = &extra[nst];
  border[6] = &u[nrs-nr];
  border[7] = &u[nrs-nr];
  border[8] = &u[nrs-1];

  border[ 9] = &extra[2*nst];
  border[10] = &extra[2*nst];
  border[11] = &extra[2*nst + ntr-nt];
  border[12] = &extra[0];
  border[13] = &u[0];
  border[14] = &extra[nst];
  border[15] = &extra[2*nst + ntr];
  border[16] = &extra[2*nst + ntr];
  border[17] = &extra[2*nst + 2*ntr-nt];

  border[18] = &u[nrst-nrs];
  border[19] = &u[nrst-nrs];
  border[20] = &u[nrst-nrs+(nr-1)];
  border[21] = &extra[nst-ns];
  border[22] = &u[nrst-nrs];
  border[23] = &extra[2*nst-ns];
  border[24] = &u[nrst-nr];
  border[25] = &u[nrst-nr];
  border[26] = &u[nrst-1];

}

#endif
