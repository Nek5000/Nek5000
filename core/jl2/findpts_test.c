#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "poly.h"
#include "gs_defs.h"
#include "comm.h"
#include "rand_elt_test.h"
#include "findpts.h"
#include "crystal.h"
#include "sarray_transfer.h"

#define NR 5
#define NS 7
#define NT 6
#define K 4
#define NEL (K*K*K)
#define TN 4

#define NPT_MAX 256
#define BBOX_TOL 0.01
#define NEWT_TOL 1024*DBL_EPSILON
#define LOC_HASH_SIZE NEL*NR*NS*NT
#define GBL_HASH_SIZE NEL*NR*NS*NT

/*
#define NPT_MAX 256
#define BBOX_TOL 1.00
#define NEWT_TOL 1024*DBL_EPSILON
#define LOC_HASH_SIZE NEL*3
#define GBL_HASH_SIZE NEL*3
*/

static uint np, id;

static const unsigned nr[3] = {NR,NS,NT};
static const unsigned mr[3] = {2*NR,2*NS,2*NT};
static double zr[NR], zs[NS], zt[NT];
static double x3[3][27];
static double mesh[3][NEL*NR*NS*NT];
static const double *const elx[3] = {mesh[0],mesh[1],mesh[2]};

struct pt_data { double x[3], r[3], dist2, ex[3]; uint code, proc, el; };
static array testp;

static crystal_data cr;

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
  const uint pn = ceil(pow(np,1/3.0));
  const uint pi=id%pn, pj=(id/pn)%pn, pk=(id/pn)/pn;
  const double pfac = 1.0/pn;
  const double pbase[3] = {-1+2*pfac*pi, -1+2*pfac*pj, -1+2*pfac*pk};
  const double fac = 1.0/K;
  const double z3[3] = {-1,0,1};
  unsigned ki,kj,kk;
  if(id==0) printf("Global division: %u^3\n",(unsigned)pn);
  if(0)
  printf("  proc %u at (%u,%u,%u) - (%g,%g,%g)\n",
         id,pi,pj,pk,pbase[0],pbase[1],pbase[2]);
  rand_elt_3(x3[0],x3[1],x3[2], z3,3,z3,3,z3,3);
  for(kk=0;kk<K;++kk) for(kj=0;kj<K;++kj) for(ki=0;ki<K;++ki) {
    unsigned off = ((kk*K+kj)*K+ki)*NR*NS*NT;
    unsigned i,j,k;
    double base[3] = {-1+2*fac*ki,-1+2*fac*kj,-1+2*fac*kk};
    for(k=0;k<NT;++k) { const double t=pbase[2]+pfac*(1+base[2]+fac*(1+zt[k]));
    for(j=0;j<NS;++j) { const double s=pbase[1]+pfac*(1+base[1]+fac*(1+zs[j]));
    for(i=0;i<NR;++i) { const double r=pbase[0]+pfac*(1+base[0]+fac*(1+zr[i]));
      mesh[0][off+(k*NS+j)*NR+i] = quad_eval(x3[0],r,s,t);
      mesh[1][off+(k*NS+j)*NR+i] = quad_eval(x3[1],r,s,t);
      mesh[2][off+(k*NS+j)*NR+i] = quad_eval(x3[2],r,s,t);
    }}}
  }
}

static void test_mesh(void)
{
  const uint pn = ceil(pow(np,1/3.0));
  const uint pi=id%pn, pj=(id/pn)%pn, pk=(id/pn)/pn;
  const double pfac = 1.0/pn;
  const double pbase[3] = {-1+2*pfac*pi, -1+2*pfac*pj, -1+2*pfac*pk};
  const double fac = 1.0/K, step = 2.0/(K*(TN-1));
  unsigned ki,kj,kk;
  struct pt_data *out = testp.ptr;
  testp.n = NEL*TN*TN*TN;
  memset(testp.ptr,0,testp.n*sizeof(struct pt_data));
  for(kk=0;kk<K;++kk) for(kj=0;kj<K;++kj) for(ki=0;ki<K;++ki) {
    unsigned i,j,k;
    double base[3] = {-1+2*fac*ki,-1+2*fac*kj,-1+2*fac*kk};
    for(k=0;k<TN;++k) { const double t = pbase[2]+pfac*(1+base[2]+step*k);
    for(j=0;j<TN;++j) { const double s = pbase[1]+pfac*(1+base[1]+step*j);
    for(i=0;i<TN;++i) { const double r = pbase[0]+pfac*(1+base[0]+step*i);
      out->x[0] = quad_eval(x3[0],r,s,t);
      out->x[1] = quad_eval(x3[1],r,s,t);
      out->x[2] = quad_eval(x3[2],r,s,t);
      out->proc = rand()%np;
      ++out;
    }}}
  }
  sarray_transfer(struct pt_data,&testp,proc,&cr);
  if(0)
    printf("%u: %u shuffled points\n",id,(unsigned)testp.n);
}

static void print_ptdata(const struct comm *const comm)
{
  uint notfound=0;
  double dist2max=0, ed2max=0;
  const struct pt_data *pt = testp.ptr, *const end = pt+testp.n;
  for(;pt!=end;++pt) {
    if(0&&id==0)
    printf("code=%u, proc=%u, el=%u, dist2=%g, r=(%.17g,%.17g,%.17g), "
           "x=(%.17g,%.17g,%.17g), ex=(%.17g,%.17g,%.17g)\n",
      pt->code,pt->proc,pt->el,pt->dist2,
      pt->r[0],pt->r[1],pt->r[2],
      pt->x[0],pt->x[1],pt->x[2],
      pt->ex[0],pt->ex[1],pt->ex[2]);
    if(pt->code==2) ++notfound;
    else {
      double ed2=0, dx;
      unsigned d; for(d=0;d<3;++d) dx=pt->x[d]-pt->ex[d], ed2+=dx*dx;
      dist2max=pt->dist2>dist2max?pt->dist2:dist2max;
      ed2max=ed2>ed2max?ed2:ed2max;
    }
  }
  {
    double distmax=sqrt(dist2max), edmax=sqrt(ed2max);
    slong total=testp.n;
    if(0)
    printf("%u: maximum distance = %g (adv), %g (eval);"
           " %u/%u points not found\n",
         (unsigned)id, distmax, edmax,
         (unsigned)notfound, (unsigned)testp.n);
    distmax = comm_reduce_double(comm,gs_max,&distmax,1);
    edmax   = comm_reduce_double(comm,gs_max,&edmax  ,1);
    notfound = comm_reduce_sint(comm,gs_add,(sint*)&notfound,1);
    total    = comm_reduce_slong(comm,gs_add,&total,1);
    if(id==0)
      printf("maximum distance = %g (adv), %g (eval);"
           " %u/%lu points not found\n",
           distmax, edmax,
           (unsigned)notfound, (unsigned long)total);
  }
}

static void test(const struct comm *const comm)
{
  const double *x_base[3];
  const unsigned x_stride[3] = {sizeof(struct pt_data),
                                sizeof(struct pt_data),
                                sizeof(struct pt_data)};
  struct findpts_data_3 *fd;
  struct pt_data *pt;
  unsigned d;
  if(id==0) printf("Initializing mesh\n");
  rand_mesh();
  test_mesh();
  pt = testp.ptr;
  if(id==0) printf("calling findpts_setup\n");
  fd=findpts_setup_3(comm,elx,nr,NEL,mr,BBOX_TOL,
                     LOC_HASH_SIZE,GBL_HASH_SIZE,
                     NPT_MAX,NEWT_TOL);
  if(id==0) printf("calling findpts\n");
  x_base[0]=pt->x, x_base[1]=pt->x+1, x_base[2]=pt->x+2;
  findpts_3(&pt->code , sizeof(struct pt_data),
            &pt->proc , sizeof(struct pt_data),
            &pt->el   , sizeof(struct pt_data),
             pt->r    , sizeof(struct pt_data),
            &pt->dist2, sizeof(struct pt_data),
             x_base   , x_stride, testp.n, fd);
  for(d=0;d<3;++d) {
    if(id==0) printf("calling findpts_eval (%u)\n",d);
    findpts_eval_3(&pt->ex[d], sizeof(struct pt_data),
                   &pt->code , sizeof(struct pt_data),
                   &pt->proc , sizeof(struct pt_data),
                   &pt->el   , sizeof(struct pt_data),
                    pt->r    , sizeof(struct pt_data),
                    testp.n, mesh[d], fd);
  }
  findpts_free_3(fd);
  print_ptdata(comm);
}

int main(int narg, char *arg[])
{
  comm_ext world;
  struct comm comm;
  
#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
#else
  world=0;
#endif

  comm_init(&comm,world);
  id=comm.id, np=comm.np;

  lobatto_nodes(zr,NR),lobatto_nodes(zs,NS),lobatto_nodes(zt,NT);
  array_init(struct pt_data,&testp,NEL*TN*TN*TN);
  crystal_init(&cr,&comm);
  test(&comm);
  crystal_free(&cr);
  array_free(&testp);
  
  comm_free(&comm);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
