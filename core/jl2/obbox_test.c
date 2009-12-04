#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "types.h"
#include "name.h"
#include "fail.h"
#include "mem.h"
#include "tensor.h"
#include "poly.h"
#include "lob_bnd.h"
#include "obbox.h"

static double det_2(const double A[4]) { return A[0]*A[3]-A[1]*A[2]; }

static double quad_2(const double x0, const double g[2], const double H[3],
                     const double r[2])
{
  return x0 + (g[0]*r[0]+g[1]*r[1])
            + (  r[0] * (H[0]*r[0]+H[1]*r[1])
               + r[1] * (H[1]*r[0]+H[2]*r[1]) )/2;
}

static void quad_2_grad(double grad[2], const double g[2], const double H[3],
                        const double r[2])
{
  grad[0] = g[0] + (H[0]*r[0]+H[1]*r[1]);
  grad[1] = g[1] + (H[1]*r[0]+H[2]*r[1]);
}

static double quad_2_jac(const double g[4], const double H[6],
                         const double r[2])
{
  double J[4];
  quad_2_grad(J  ,g  ,H  ,r);
  quad_2_grad(J+2,g+2,H+3,r);
  return det_2(J);
}

static double det_3(const double A[9])
{
  const double a = A[4]*A[8]-A[5]*A[7],
               b = A[5]*A[6]-A[3]*A[8],
               c = A[3]*A[7]-A[4]*A[6];
  return A[0]*a+A[1]*b+A[2]*c;
}

static double quad_3(const double x0, const double g[3], const double H[6],
                     const double r[3])
{
  return x0 + (g[0]*r[0]+g[1]*r[1]+g[2]*r[2])
            + (  r[0] * (H[0]*r[0]+H[1]*r[1]+H[2]*r[2])
               + r[1] * (H[1]*r[0]+H[3]*r[1]+H[4]*r[2])
               + r[2] * (H[2]*r[0]+H[4]*r[1]+H[5]*r[2]) )/2;
}

static void quad_3_grad(double grad[3], const double g[3], const double H[6],
                        const double r[3])
{
  grad[0] = g[0] + (H[0]*r[0]+H[1]*r[1]+H[2]*r[2]);
  grad[1] = g[1] + (H[1]*r[0]+H[3]*r[1]+H[4]*r[2]);
  grad[2] = g[2] + (H[2]*r[0]+H[4]*r[1]+H[5]*r[2]);
}

static double quad_3_jac(const double g[9], const double H[18],
                         const double r[3])
{
  double J[9];
  quad_3_grad(J  ,g  ,H   ,r);
  quad_3_grad(J+3,g+3,H+ 6,r);
  quad_3_grad(J+6,g+6,H+12,r);
  return det_3(J);
}

static void rand_elt_2(double *x, double *y,
                       const double *zr, unsigned nr,
                       const double *zs, unsigned ns)
{
  static int init=0;
  static double z4[4], lob_bnd_data[16+3*4*(2*16+1)],
                work[2*16*(4+16+1)];
  unsigned i,j;
  double x0[2], g[4], H[6], jac[4*4], r[2];
  struct dbl_range jr;
  if(!init) {
    init=1;
    lobatto_nodes(z4,4);
    lob_bnd_setup(lob_bnd_data,4,16);
  }
  do {
    for(i=0;i<4;++i) g[i] = -1+2*(rand()/(double)RAND_MAX);
    for(i=0;i<6;++i) H[i] =.5*(-1+2*(rand()/(double)RAND_MAX));
    for(j=0;j<4;++j) { r[1] = z4[j];
      for(i=0;i<4;++i) { r[0] = z4[i];
        jac[j*4+i] = quad_2_jac(g,H,r);
      }
    }
    jr = lob_bnd_2(lob_bnd_data,4,16, lob_bnd_data,4,16, jac, work);
    /*printf("Jacobian range %g, %g\n", jr.min, jr.max);*/
  } while(jr.max*jr.min<=0);
  for(i=0;i< 2;++i) x0[i] = -1+2*(rand()/(double)RAND_MAX);
  for(j=0;j<ns;++j) {   r[1] = zs[j];
    for(i=0;i<nr;++i) { r[0] = zr[i];
      x[j*nr+i] = quad_2(x0[0],g  ,H  ,r);
      y[j*nr+i] = quad_2(x0[1],g+2,H+3,r);
    }
  }
}

static void rand_elt_3(double *x, double *y, double *z,
                       const double *zr, unsigned nr,
                       const double *zs, unsigned ns,
                       const double *zt, unsigned nt)
{
  static int init=0;
  static double z4[4], lob_bnd_data[16+3*4*(2*16+1)],
                work[2*16*16*(4+16+1)];
  unsigned i,j,k;
  double x0[3], g[9], H[18], jac[4*4*4], r[3];
  struct dbl_range jr;
  if(!init) {
    init=1;
    lobatto_nodes(z4,4);
    lob_bnd_setup(lob_bnd_data,4,16);
  }
  do {
    for(i=0;i< 9;++i) g[i] = -1+2*(rand()/(double)RAND_MAX);
    for(i=0;i<18;++i) H[i] =.5*(-1+2*(rand()/(double)RAND_MAX));
    for(k=0;k<4;++k) { r[2] = z4[k];
      for(j=0;j<4;++j) { r[1] = z4[j];
        for(i=0;i<4;++i) { r[0] = z4[i];
          jac[(k*4+j)*4+i] = quad_3_jac(g,H,r);
        }
      }
    }
    jr = lob_bnd_3(lob_bnd_data,4,16, lob_bnd_data,4,16, lob_bnd_data,4,16,
                   jac, work);
    /*printf("Jacobian range %g, %g\n", jr.min, jr.max);*/
  } while(jr.max*jr.min<=0);
  for(i=0;i< 3;++i) x0[i] = -1+2*(rand()/(double)RAND_MAX);
  for(k=0;k<nt;++k) {     r[2] = zt[k];
    for(j=0;j<ns;++j) {   r[1] = zs[j];
      for(i=0;i<nr;++i) { r[0] = zr[i];
        x[(k*ns+j)*nr+i] = quad_3(x0[0],g  ,H   ,r);
        y[(k*ns+j)*nr+i] = quad_3(x0[1],g+3,H+ 6,r);
        z[(k*ns+j)*nr+i] = quad_3(x0[2],g+6,H+12,r);
      }
    }
  }
}

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

static void bubble_elt(double *x, double *y, double *z,
                       const double *zr, unsigned nr,
                       const double *zs, unsigned ns,
                       const double *zt, unsigned nt, int type)
{
  unsigned i,j,k;
  for(k=0;k<nt;++k) for(j=0;j<ns;++j) for(i=0;i<nr;++i) {
    double dx=0,dy=0,dz=0;
    switch(type) {
      case 0: dx =  cos(PI*zs[j]/2)*cos(PI*zt[k]/2); break;
      case 1: dx = -cos(PI*zs[j]/2)*cos(PI*zt[k]/2); break;
      case 2: dy =  cos(PI*zt[k]/2)*cos(PI*zr[i]/2); break;
      case 3: dy = -cos(PI*zt[k]/2)*cos(PI*zr[i]/2); break;
      case 4: dz =  cos(PI*zr[i]/2)*cos(PI*zs[j]/2); break;
      case 5: dz = -cos(PI*zr[i]/2)*cos(PI*zs[j]/2); break;
    }
    x[(k*ns+j)*nr+i] = zr[i] + dx;
    y[(k*ns+j)*nr+i] = zs[j] + dy;
    z[(k*ns+j)*nr+i] = zt[k] + dz;
  }
}

#define REPEAT 20

#define N 100
#define NR 7
#define MR (4*NR)
#define NS 8
#define MS (4*NS)
#define NT 9
#define MT (4*NT)

#define TOL 0.00001

static double zr[NR], zs[NS], zt[NT];
static double x[NR*NS*NT*N], y[NR*NS*NT*N], z[NR*NS*NT*N];
static double tx[3][NR*NS*NT];

static obbox_2 ob2[N*NT];
static obbox_3 ob3[N];

static struct dbl_range dbl_range_expand(struct dbl_range b, double tol)
{
  double a = (b.min+b.max)/2, l = (b.max-b.min)*(1+tol)/2;
  struct dbl_range m = { a-l,a+l };
  return m;
}

int main()
{
  int failure=0;
  unsigned i;
  
  double *lob_bnd_data_r = tmalloc(double,
      lob_bnd_size(NR,MR)+lob_bnd_size(NS,MS)+lob_bnd_size(NT,MT)),
    *lob_bnd_data_s = lob_bnd_data_r + lob_bnd_size(NR,MR),
    *lob_bnd_data_t = lob_bnd_data_s + lob_bnd_size(NS,MS);

  lobatto_nodes(zr,NR); lob_bnd_setup(lob_bnd_data_r,NR,MR);
  lobatto_nodes(zs,NS); lob_bnd_setup(lob_bnd_data_s,NS,MS);
  lobatto_nodes(zt,NT); lob_bnd_setup(lob_bnd_data_t,NT,MT);

  /* 2-D */
  for(i=0;i<REPEAT;++i) {
    unsigned n; double *x_ = x, *y_ = y;
    for(n=0;n<N && n<6;++n, x_+=NR*NS*NT, y_+=NR*NS*NT)
      bubble_elt(x_,y_,z, zr,NR, zs,NS, zt,NT, n);
    for(n=N-6;n;--n, x_+=NR*NS*NT, y_+=NR*NS*NT)
      rand_elt_2(x_,y_, zr,NR, zs,NS);
    obbox_2_calc(ob2, x,y, NR,NS,NT*N, MR,MS, TOL);
    x_=x, y_=y;
    for(n=0;n<N*NT;++n, x_+=NR*NS, y_+=NR*NS) {
      const obbox_2 *ob = &ob2[n];
      struct dbl_range xr,yr, tr[2];
      static double work[2*MR*(NS+MS+1)];
      unsigned j;
      for(j=0;j<NR*NS;++j) {
        const double dx=x_[j]-ob->c0[0], dy=y_[j]-ob->c0[1];
        tx[0][j] = ob->A[0]*dx+ob->A[1]*dy;
        tx[1][j] = ob->A[2]*dx+ob->A[3]*dy;
        if(   (x_[j]-ob->x.min)*(ob->x.max-x_[j]) < 0
           || (y_[j]-ob->y.min)*(ob->y.max-y_[j]) < 0 )
          failure=1,
          printf("%d %d (%g,%g) not in [%g,%g] x [%g,%g]\n", n, j,
            x_[j],y_[j], ob->x.min,ob->x.max, ob->y.min,ob->y.max);
        if(   (tx[0][j]+1)*(1-tx[0][j]) < 0
           || (tx[1][j]+1)*(1-tx[1][j]) < 0 )
          failure=1,
          printf("%d %d (%g,%g) not in [-1,1]^2\n", n, j,
            tx[0][j],tx[1][j]);
        if(failure) break;
      }
      
      xr = dbl_range_expand(lob_bnd_2(lob_bnd_data_r,NR,MR,
                                      lob_bnd_data_s,NS,MS, x_, work), TOL);
      yr = dbl_range_expand(lob_bnd_2(lob_bnd_data_r,NR,MR,
                                      lob_bnd_data_s,NS,MS, y_, work), TOL);

      for(j=0;j<2;++j) tr[j] = dbl_range_expand(
          lob_bnd_2(lob_bnd_data_r,NR,MR, lob_bnd_data_s,NS,MS, tx[j], work)
        , TOL);
        
      if(   ob->x.min < xr.min - DBL_EPSILON*128
         || ob->x.max > xr.max + DBL_EPSILON*128 ) failure = 1;
      if(   ob->y.min < yr.min - DBL_EPSILON*128
         || ob->y.max > yr.max + DBL_EPSILON*128 ) failure = 1;
      
      for(j=0;j<2;++j)
        if(   tr[j].min > -1 + DBL_EPSILON*128
           || tr[j].max <  1 - DBL_EPSILON*128 ) failure = 1;
           
      if((i==0&&n==0) || failure) {
        printf("x: [%g,%g] in [%g,%g]\n", ob->x.min, ob->x.max,
                                          xr.min, xr.max);
        printf("y: [%g,%g] in [%g,%g]\n", ob->y.min, ob->y.max,
                                          yr.min, yr.max);
        for(j=0;j<2;++j)
          printf("r %d: [%g,%g]\n", j, tr[j].min, tr[j].max);
      }
      if(failure) break;
    }
    if(failure) break;
    printf("."); fflush(stdout);
  }
  printf("\n");

  /* 3-D */
  for(i=0;!failure && i<REPEAT;++i) {
    unsigned n; double *x_ = x, *y_ = y, *z_ = z;
    for(n=0;n<N && n<6;++n, x_+=NR*NS*NT, y_+=NR*NS*NT, z_+=NR*NS*NT)
      bubble_elt(x_,y_,z_, zr,NR, zs,NS, zt,NT, n);
    for(n=N-6;n;--n, x_+=NR*NS*NT, y_+=NR*NS*NT, z_+=NR*NS*NT)
      rand_elt_3(x_,y_,z_, zr,NR, zs,NS, zt,NT);
    obbox_3_calc(ob3, x,y,z, NR,NS,NT,N, MR,MS,MT, TOL);
    x_=x, y_=y, z_=z;
    for(n=0;n<N;++n, x_+=NR*NS*NT, y_+=NR*NS*NT, z_+=NR*NS*NT) {
      const obbox_3 *ob = &ob3[n];
      struct dbl_range xr,yr,zr, tr[3];
      static double work[2*MR*MS*(NT+MT+1)];
      unsigned j;
      for(j=0;j<NR*NS*NT;++j) {
        const double dx=x_[j]-ob->c0[0], dy=y_[j]-ob->c0[1], dz=z_[j]-ob->c0[2];
        tx[0][j] = ob->A[0]*dx+ob->A[1]*dy+ob->A[2]*dz;
        tx[1][j] = ob->A[3]*dx+ob->A[4]*dy+ob->A[5]*dz;
        tx[2][j] = ob->A[6]*dx+ob->A[7]*dy+ob->A[8]*dz;
        if(   (x_[j]-ob->x.min)*(ob->x.max-x_[j]) < 0
           || (y_[j]-ob->y.min)*(ob->y.max-y_[j]) < 0
           || (z_[j]-ob->z.min)*(ob->z.max-z_[j]) < 0 )
          failure=1,
          printf("%d %d (%g,%g,%g) not in [%g,%g] x [%g,%g] x [%g,%g]\n", n, j,
            x_[j],y_[j],z_[j], ob->x.min,ob->x.max, ob->y.min,ob->y.max,
            ob->z.min,ob->z.max);
        if(   (tx[0][j]+1)*(1-tx[0][j]) < 0
           || (tx[1][j]+1)*(1-tx[1][j]) < 0
           || (tx[2][j]+1)*(1-tx[2][j]) < 0 )
          failure=1,
          printf("%d %d (%g,%g,%g) not in [-1,1]^3\n", n, j,
            tx[0][j],tx[1][j],tx[2][j]);
        if(failure) break;
      }
      
      xr = dbl_range_expand(lob_bnd_3(lob_bnd_data_r,NR,MR,
                                      lob_bnd_data_s,NS,MS,
                                      lob_bnd_data_t,NT,MT, x_, work), TOL);
      yr = dbl_range_expand(lob_bnd_3(lob_bnd_data_r,NR,MR,
                                      lob_bnd_data_s,NS,MS,
                                      lob_bnd_data_t,NT,MT, y_, work), TOL);
      zr = dbl_range_expand(lob_bnd_3(lob_bnd_data_r,NR,MR,
                                      lob_bnd_data_s,NS,MS,
                                      lob_bnd_data_t,NT,MT, z_, work), TOL);

      for(j=0;j<3;++j) tr[j] = dbl_range_expand(
          lob_bnd_3(lob_bnd_data_r,NR,MR, lob_bnd_data_s,NS,MS,
                    lob_bnd_data_t,NT,MT, tx[j], work)
        , TOL);
        
      if(   ob->x.min < xr.min - DBL_EPSILON*128
         || ob->x.max > xr.max + DBL_EPSILON*128 ) failure = 1;
      if(   ob->y.min < yr.min - DBL_EPSILON*128
         || ob->y.max > yr.max + DBL_EPSILON*128 ) failure = 1;
      if(   ob->z.min < zr.min - DBL_EPSILON*128
         || ob->z.max > zr.max + DBL_EPSILON*128 ) failure = 1;
      
      for(j=0;j<3;++j)
        if(   tr[j].min > -1 + DBL_EPSILON*128
           || tr[j].max <  1 - DBL_EPSILON*128 ) failure = 1;
           
      if((i==0&&n==0) || failure) {
        printf("x: [%g,%g] in [%g,%g]\n", ob->x.min, ob->x.max,
                                          xr.min, xr.max);
        printf("y: [%g,%g] in [%g,%g]\n", ob->y.min, ob->y.max,
                                          yr.min, yr.max);
        printf("z: [%g,%g] in [%g,%g]\n", ob->z.min, ob->z.max,
                                          zr.min, zr.max);
        for(j=0;j<3;++j)
          printf("r %d: [%g,%g]\n", j, tr[j].min, tr[j].max);
      }
      if(failure) break;
    }
    if(failure) break;
    printf("."); fflush(stdout);
  }
  printf("\n");
        
  free(lob_bnd_data_r);
 
  printf("Tests %s\n", failure?"failed":"successful");

  return failure;
}
