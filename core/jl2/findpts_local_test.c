#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "name.h"
#include "fail.h"
#include "mem.h"
#include "types.h"
#include "poly.h"
#include "obbox.h"
#include "findpts_el.h"
#include "findpts_local.h"
#include "rand_elt_test.h"

#define NR 5
#define NS 8
#define NT 6
#define K 4
#define NEL (K*K*K)
#define TN 4

#define NPT_MAX 256
#define BBOX_TOL 0.01
#define NEWT_TOL 1024*DBL_EPSILON
#define MAX_HASH_SIZE NEL*NR*NS*NT

/*
#define NPT_MAX 256
#define BBOX_TOL 1.00
#define NEWT_TOL 1024*DBL_EPSILON
#define MAX_HASH_SIZE NEL*3
*/

static const unsigned nr[3] = {NR,NS,NT};
static const unsigned mr[3] = {4*NR,4*NS,4*NT};
static double zr[NR], zs[NS], zt[NT];
static double x3[3][27];
static double mesh[3][NEL*NR*NS*NT];
static const double *const elx[3] = {mesh[0],mesh[1],mesh[2]};

static double testx[NEL*TN*TN*TN*3];
struct pt_data { double r[3], dist2; uint code, el; };
static struct pt_data testp[NEL*TN*TN*TN];

static double quad_eval(const double coef[27],
                        const double r, const double s, const double t)
{
  const double lr0=r*(r-1)/2, lr1=(1+r)*(1-r), lr2=r*(r+1)/2;
  const double ls0=s*(s-1)/2, ls1=(1+s)*(1-s), ls2=s*(s+1)/2;
  const double lt0=t*(t-1)/2, lt1=(1+t)*(1-t), lt2=t*(t+1)/2;
  return ( (coef[ 0]*lr0+coef[ 1]*lr1+coef[ 2]*lr2)*ls0
          +(coef[ 3]*lr0+coef[ 4]*lr1+coef[ 5]*lr2)*ls1
          +(coef[ 6]*lr0+coef[ 7]*lr1+coef[ 8]*lr2)*ls2)*lt0
        +( (coef[ 9]*lr0+coef[10]*lr1+coef[11]*lr2)*ls0
          +(coef[12]*lr0+coef[13]*lr1+coef[14]*lr2)*ls1
          +(coef[15]*lr0+coef[16]*lr1+coef[17]*lr2)*ls2)*lt1
        +( (coef[18]*lr0+coef[19]*lr1+coef[20]*lr2)*ls0
          +(coef[21]*lr0+coef[22]*lr1+coef[23]*lr2)*ls1
          +(coef[24]*lr0+coef[25]*lr1+coef[26]*lr2)*ls2)*lt2;
}

static void rand_mesh(void)
{
  const double fac = 1.0/K;
  const double z3[3] = {-1,0,1};
  unsigned ki,kj,kk;
  rand_elt_3(x3[0],x3[1],x3[2], z3,3,z3,3,z3,3);
  for(kk=0;kk<K;++kk) for(kj=0;kj<K;++kj) for(ki=0;ki<K;++ki) {
    unsigned off = ((kk*K+kj)*K+ki)*NR*NS*NT;
    unsigned i,j,k;
    double base[3] = {-1+2*fac*ki,-1+2*fac*kj,-1+2*fac*kk};
    for(k=0;k<NT;++k) { const double t = base[2]+fac*(1+zt[k]);
    for(j=0;j<NS;++j) { const double s = base[1]+fac*(1+zs[j]);
    for(i=0;i<NR;++i) { const double r = base[0]+fac*(1+zr[i]);
      mesh[0][off+(k*NS+j)*NR+i] = quad_eval(x3[0],r,s,t);
      mesh[1][off+(k*NS+j)*NR+i] = quad_eval(x3[1],r,s,t);
      mesh[2][off+(k*NS+j)*NR+i] = quad_eval(x3[2],r,s,t);
    }}}
  }
}

static void test_mesh(void)
{
  const double fac = 1.0/K, step = 2.0/(K*(TN-1));
  unsigned ki,kj,kk;
  for(kk=0;kk<K;++kk) for(kj=0;kj<K;++kj) for(ki=0;ki<K;++ki) {
    unsigned off = ((kk*K+kj)*K+ki)*TN*TN*TN;
    unsigned i,j,k;
    double base[3] = {-1+2*fac*ki,-1+2*fac*kj,-1+2*fac*kk};
    for(k=0;k<TN;++k) { const double t = base[2]+step*k;
    for(j=0;j<TN;++j) { const double s = base[1]+step*j;
    for(i=0;i<TN;++i) { const double r = base[0]+step*i;
      testx[(off+(k*TN+j)*TN+i)*3+0] = quad_eval(x3[0],r,s,t);
      testx[(off+(k*TN+j)*TN+i)*3+1] = quad_eval(x3[1],r,s,t);
      testx[(off+(k*TN+j)*TN+i)*3+2] = quad_eval(x3[2],r,s,t);
    }}}
  }
}

static void print_ptdata(void)
{
  uint i,notfound=0;
  double dist2max=0;
  for(i=0;i<NEL*TN*TN*TN;++i) {
    printf("code=%u, el=%u, dist2=%g, r=(%.17g,%.17g,%.17g)\n",
      testp[i].code,testp[i].el,testp[i].dist2,
      testp[i].r[0],testp[i].r[1],testp[i].r[2]);
    dist2max=testp[i].dist2>dist2max?testp[i].dist2:dist2max;
    if(testp[i].code==2) ++notfound;
  }
  printf("Maximum distance = %g\n%u points not found\n",
         sqrt(dist2max), (unsigned)notfound);
}

static void test(buffer *buf)
{
  const double *const x_base[3]={testx,testx+1,testx+2};
  const unsigned x_stride[3]=
    {3*sizeof(double),3*sizeof(double),3*sizeof(double)};
  struct findpts_local_data_3 fld;
  rand_mesh();
  test_mesh();
  findpts_local_setup_3(&fld,elx,nr,NEL,mr,BBOX_TOL,MAX_HASH_SIZE,
                        NPT_MAX,NEWT_TOL);
  findpts_local_3(&testp[0].code , sizeof(struct pt_data),
                  &testp[0].el   , sizeof(struct pt_data),
                   testp[0].r    , sizeof(struct pt_data),
                  &testp[0].dist2, sizeof(struct pt_data),
                  x_base, x_stride,
                  NEL*TN*TN*TN, &fld, buf);
  findpts_local_free_3(&fld);
  print_ptdata();
}

int main()
{
  buffer buf = null_buffer;
  lobatto_nodes(zr,NR),lobatto_nodes(zs,NS),lobatto_nodes(zt,NT);
  test(&buf);
  buffer_free(&buf);
  return 0;
}
