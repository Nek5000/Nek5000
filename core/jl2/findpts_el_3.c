#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "tensor.h"
#include "poly.h"

#define findpts_el_3_setup   PREFIXED_NAME(findpts_el_3_setup )
#define findpts_el_3_free    PREFIXED_NAME(findpts_el_3_free  )
#define findpts_el_3         PREFIXED_NAME(findpts_el_3       )

/* A is row-major */
static void lin_solve_3(double x[3], const double A[9], const double y[3])
{
  const double a = A[4]*A[8]-A[5]*A[7],
               b = A[5]*A[6]-A[3]*A[8],
               c = A[3]*A[7]-A[4]*A[6],
            idet = 1/(A[0]*a+A[1]*b+A[2]*c);
  const double
    inv0 = a,
    inv1 = A[2]*A[7]-A[1]*A[8],
    inv2 = A[1]*A[5]-A[2]*A[4],
    inv3 = b,
    inv4 = A[0]*A[8]-A[2]*A[6],
    inv5 = A[2]*A[3]-A[0]*A[5],
    inv6 = c,
    inv7 = A[1]*A[6]-A[0]*A[7],
    inv8 = A[0]*A[4]-A[1]*A[3];
  x[0] = idet*(inv0*y[0] + inv1*y[1] + inv2*y[2]);
  x[1] = idet*(inv3*y[0] + inv4*y[1] + inv5*y[2]);
  x[2] = idet*(inv6*y[0] + inv7*y[1] + inv8*y[2]);
}

static void lin_solve_sym_2(double x[2], const double A[3], const double y[2])
{
  const double idet = 1/(A[0]*A[2] - A[1]*A[1]);
  x[0] = idet*(A[2]*y[0] - A[1]*y[1]);
  x[1] = idet*(A[0]*y[1] - A[1]*y[0]);
}


struct findpts_el_3_pt { double x[3], r[3]; double dist2; unsigned flags; };

/* the bit structure of flags is CTTSSRR
   the C bit --- 1<<6 --- is set when the point is converged
   RR is 0 = 00b if r is unconstrained,
         1 = 01b if r is constrained at -1
         2 = 10b if r is constrained at +1
   SS, TT are similarly for s and t constraints
*/

static unsigned num_constrained(const unsigned flags)
{
  const unsigned y = flags | flags>>1;
  return (y&1u) + (y>>2 & 1u) + (y>>4 & 1u);
}

static unsigned pt_flags_to_bin_noC(const unsigned flags)
{
  return ((flags>>4 & 3u)*3 + (flags>>2 & 3u))*3 + (flags & 3u);
}

/* map flags to 27 if the C bit is set,
   else to [0,26] --- the 27 valid configs of TTSSRR */
static unsigned pt_flags_to_bin(const unsigned flags)
{
  const unsigned mask = 0u - (flags>>6); /* 0 or 0xfff... when converged */
  return (mask & 27u) | (~mask & pt_flags_to_bin_noC(flags));
}

/* assumes x = 0, 1, or 2 */
static unsigned plus_1_mod_3(const unsigned x) { return ((x | x>>1)+1) & 3u; } 
static unsigned plus_2_mod_3(const unsigned x)
{
  const unsigned y = (x-1) & 3u;
  return y ^ (y>>1);
}

/* assumes x = 1 << i, with i < 6, returns i+1 */
static unsigned which_bit(const unsigned x)
{
  const unsigned y = x&7u;
  return (y-(y>>2)) | ((x-1)&4u) | (x>>4);
}

static unsigned face_index(const unsigned x) { return which_bit(x)-1; }

static unsigned edge_index(const unsigned x)
{
  const unsigned y = ~((x>>1) | x);
  const unsigned RTSR = ((x>>1)&1u) | ((x>>2)&2u) | ((x>>3)&4u) | ((x<<2)&8u);
  const unsigned re = RTSR>>1;
  const unsigned se = 4u | RTSR>>2;
  const unsigned te = 8u | (RTSR&3u);
  return   ( (0u - ( y    &1u)) & re )
         | ( (0u - ((y>>2)&1u)) & se )
         | ( (0u - ((y>>4)&1u)) & te );
}

static unsigned point_index(const unsigned x)
{
  return ((x>>1)&1u) | ((x>>2)&2u) | ((x>>3)&4u);
}

/* extra data

  we need x, dx/dn for each face
    rs: x at 0, nrst - nrs,
      6*nrs extra for dx/dn
    st: 12*nst extra
    tr: 12*ntr extra
      (transposed order for embedded t-edges)

  for each edge,
    have x, dx/dn2 already as part of face data
    need dx/dn1 (strided in face data)
    need d^2x/dn1^2, d^2x/dn2^2 possibly, if constraints relax
      thats 3*4*(nr+ns+nt) extra

*/

struct findpts_el_3_gface { const double *x[3], *dxdn[3]; };
struct findpts_el_3_gedge { const double *x[3], *dxdn1[3], *dxdn2[3],
                                         *d2xdn1[3], *d2xdn2[3]; };
struct findpts_el_3_gpt   { double x[3], jac[9], hes[18]; };

struct findpts_el_3_data {
  struct findpts_el_3_pt *p;

  unsigned n[3];
  lagrange_fun *lag[3];
  double *lag_data[3];
  double *wtend[3];
  
  const double *x[3];
  
  unsigned side_init;
  double *sides;
  struct findpts_el_3_gface face[6]; /* ST R=-1,R=+1; TR S=-1,S=+1; ... */
  struct findpts_el_3_gedge edge[12]; /* R S=-1,T=-1; R S=1,T=-1; ... */
  struct findpts_el_3_gpt pt[8];

  double *work;
};

/* work[2*nt+2*nrs] */
/* work[4*nr+4*nst] */
/* work[4*ns+4*nr] */
/* work[4*n1+4*n], work[2*n2+2*n] */
/* work[4*nr+4], work[2*nt+2] */
/* work[(3+9+2*(nr+ns+nt+nrs))*pn + max(2*nr,ns) ] */
/* work[(3+9+3+3*(n1+n2+n1))*pn ] */
/* work[ 3*n ] */
static uint work_size(
  const unsigned nr, const unsigned ns, const unsigned nt,
  const uint npt_max)
{
  unsigned n1, n2;
  uint wsize;
  if(nr>ns) {
    if(nr>nt) n1=nr, n2 = (ns>nt ? ns : nt);
    else      n1=nt, n2 = nr;
  } else {
    if(ns>nt) n1=ns, n2 = (nr>nt ? nr : nt);
    else      n1=nt, n2 = ns;
  }
  #define DO_MAX(x) do { const uint temp=(x); \
                         wsize=temp>wsize?temp:wsize; } while(0)
  wsize = (12 + 2*(nr+ns+nt+nr*ns)) * npt_max + (2*nr>ns?2*nr:ns);
  DO_MAX(2*(nt+nr*ns));
  DO_MAX(4*(nr+ns*nt));
  DO_MAX(4*(n1+n2));
  DO_MAX(npt_max*(15+3*(2*n1+n2)));
  #undef DO_MAX
  return wsize;
}

struct findpts_el_3_data *findpts_el_3_setup(
  const unsigned nr, const unsigned ns, const unsigned nt,
  const uint npt_max)
{
  struct findpts_el_3_data *fd = tmalloc(struct findpts_el_3_data,1);
  const unsigned nrs = nr*ns, nst=ns*nt, ntr=nt*nr;
  const unsigned face_size = 12*nst + 12*ntr + 6*nrs;
  const unsigned off_es = face_size + 36*nr, off_et = off_es + 36*ns,
                 tot = off_et + 36*nt;
  unsigned d,i, lag_size[3];

  fd->p = tmalloc(struct findpts_el_3_pt, npt_max*2);

  fd->n[0]=nr, fd->n[1]=ns, fd->n[2]=nt;
  for(d=0;d<3;++d) lag_size[d] = gll_lag_size(fd->n[d]);
  
  fd->lag_data[0] = tmalloc(double,lag_size[0]+lag_size[1]+lag_size[2]
                                   +6*(nr+ns+nt) + tot +
                                   work_size(nr,ns,nt,npt_max));
  fd->lag_data[1] = fd->lag_data[0]+lag_size[0];
  fd->lag_data[2] = fd->lag_data[1]+lag_size[1];
  fd->wtend[0]    = fd->lag_data[2]+lag_size[2];
  fd->wtend[1]    = fd->wtend[0]+6*nr;
  fd->wtend[2]    = fd->wtend[1]+6*ns;
  fd->sides       = fd->wtend[2]+6*nt;
  fd->work        = fd->sides + tot;

  fd->side_init = 0;
  
  for(d=0;d<3;++d) {
    double *wt=fd->wtend[d]; unsigned n=fd->n[d];
    fd->lag[d] = gll_lag_setup(fd->lag_data[d],n);
    fd->lag[d](wt    , fd->lag_data[d],n,2,-1);
    fd->lag[d](wt+3*n, fd->lag_data[d],n,2, 1);
    
    wt[0]=1; for(i=1;i<n;++i) wt[i]=0;
    wt+=3*n; { for(i=0;i<n-1;++i) wt[i]=0; } wt[i]=1;
  }

  #define SET_FACE(i,base,n) do { for(d=0;d<3;++d) \
    fd->face[2*i  ].x[d]    = fd->sides + base +    d *n, \
    fd->face[2*i  ].dxdn[d] = fd->sides + base + (3+d)*n, \
    fd->face[2*i+1].x[d]    = fd->sides + base + (6+d)*n, \
    fd->face[2*i+1].dxdn[d] = fd->sides + base + (9+d)*n; \
  } while(0)
  SET_FACE(0,0,nst);
  SET_FACE(1,12*nst,ntr);
  #undef SET_FACE

  for(d=0;d<3;++d)
    fd->face[4].x[d] = 0, /* will point to user data */
    fd->face[4].dxdn[d] = fd->sides + 12*(nst+ntr) + d*nrs,
    fd->face[5].x[d] = 0, /* will point to user data */
    fd->face[5].dxdn[d] = fd->sides + 12*(nst+ntr) + (3+d)*nrs;

  #define SET_EDGE1(j,d,rd,rn,base) \
    for(i=0;i<2;++i) \
      fd->edge[4*j+2*i+0].dxdn2[d] = fd->face[2*j+i].dxdn[d], \
      fd->edge[4*j+2*i+1].dxdn2[d] = fd->face[2*j+i].dxdn[d]+n##rd##rn-n##rn;
  #define SET_EDGE2(j,d,rd,rn,base) \
    for(i=0;i<4;++i) \
      fd->edge[4*j+i].dxdn1 [d] = fd->sides + base + (9*i  +d)*n##rd, \
      fd->edge[4*j+i].d2xdn1[d] = fd->sides + base + (9*i+3+d)*n##rd, \
      fd->edge[4*j+i].d2xdn2[d] = fd->sides + base + (9*i+6+d)*n##rd;
  #define SET_EDGE(j,rd,rn,base) do { \
    for(d=0;d<3;++d) { SET_EDGE1(j,d,rd,rn,base); \
                       SET_EDGE2(j,d,rd,rn,base); } \
  } while(0)
  SET_EDGE(0,r,s,face_size);
  SET_EDGE(1,s,t,off_es);
  SET_EDGE(2,t,r,off_et);
  #undef SET_EDGE
  #undef SET_EDGE2
  #undef SET_EDGE1
  
  return fd;
}

void findpts_el_3_free(struct findpts_el_3_data *fd)
{
  free(fd->p);
  free(fd->lag_data[0]);
  free(fd);
}

/* work[2*nt+2*nrs] */
static void compute_face_data_rs(struct findpts_el_3_data *fd)
{
  const unsigned nr = fd->n[0], ns=fd->n[1], nt=fd->n[2],
                 nrs = nr*ns, nst=ns*nt, ntr = nt*nr, nrstm1 = nrs*(nt-1);
  unsigned d;
  double *work = fd->work, *out = fd->sides + 12*(nst+ntr);
  memcpy(work   , fd->wtend[2]+  nt, nt*sizeof(double));
  memcpy(work+nt, fd->wtend[2]+4*nt, nt*sizeof(double));
  for(d=0;d<3;++d) {
    tensor_mxm(work+2*nt,nrs, fd->x[d],nt, work,2);
    memcpy(out+   d *nrs, work+2*nt       , nrs*sizeof(double));
    memcpy(out+(3+d)*nrs, work+2*nt+nrs   , nrs*sizeof(double));
    fd->face[4].x[d] = fd->x[d];
    fd->face[5].x[d] = fd->x[d] + nrstm1;
  }
}

/* work[4*nr+4*nst] */
static void compute_face_data_st(struct findpts_el_3_data *fd)
{
  const unsigned nr = fd->n[0], ns=fd->n[1], nt=fd->n[2], nst=ns*nt;
  unsigned i;
  double *work = fd->work, *out = fd->sides;
  memcpy(work     , fd->wtend[0]     , 2*nr*sizeof(double));
  memcpy(work+2*nr, fd->wtend[0]+3*nr, 2*nr*sizeof(double));
  for(i=0;i<3;++i) {
    tensor_mtxm(work+4*nr,nst, fd->x[i],nr, work,4);
    memcpy(out+   i *nst, work+4*nr      , nst*sizeof(double));
    memcpy(out+(3+i)*nst, work+4*nr+  nst, nst*sizeof(double));
    memcpy(out+(6+i)*nst, work+4*nr+2*nst, nst*sizeof(double));
    memcpy(out+(9+i)*nst, work+4*nr+3*nst, nst*sizeof(double));
  }
}

/* work[4*ns+4*nr] */
static void compute_face_data_tr(struct findpts_el_3_data *fd)
{
  const unsigned nr = fd->n[0], ns=fd->n[1], nt=fd->n[2],
                 nrs = nr*ns, nst=ns*nt, ntr=nt*nr;
  unsigned i,k,d;
  double *work = fd->work, *out = fd->sides + 12*nst;
  memcpy(work     , fd->wtend[1]     , 2*ns*sizeof(double));
  memcpy(work+2*ns, fd->wtend[1]+3*ns, 2*ns*sizeof(double));
  for(d=0;d<3;++d) {
    for(k=0;k<nt;++k) {
      double *outk; double *in = work+4*ns;
      tensor_mxm(in,nr, fd->x[d]+k*nrs,ns, work,4);
      for(outk=out+   d *ntr+k,i=0;i<nr;++i,outk+=nt) *outk=*in++;
      for(outk=out+(3+d)*ntr+k,i=0;i<nr;++i,outk+=nt) *outk=*in++;
      for(outk=out+(6+d)*ntr+k,i=0;i<nr;++i,outk+=nt) *outk=*in++;
      for(outk=out+(9+d)*ntr+k,i=0;i<nr;++i,outk+=nt) *outk=*in++;
    }
  }
}

static const struct findpts_el_3_gface *get_face(
  struct findpts_el_3_data *fd, unsigned fi)
{
  const unsigned mask = 1<<(fi/2);
  if((fd->side_init&mask)==0) {
    switch(fi/2) {
      case 1: compute_face_data_rs(fd); break;
      case 2: compute_face_data_st(fd); break;
      case 3: compute_face_data_tr(fd); break;
    }
    fd->side_init |= mask;
  }
  return &fd->face[fi];
}

/* work[4*n1+4*n], work[2*n2+2*n] */
static void compute_edge_data(struct findpts_el_3_data *fd, unsigned d)
{
  const unsigned dn1 = d  +1>=3?0:d  +1;
  const unsigned dn2 = dn1+1>=3?0:dn1+1;
  const unsigned n = fd->n[d], n1 = fd->n[dn1], n2 = fd->n[dn2];
  const unsigned nr=fd->n[0],ns=fd->n[1],nt=fd->n[2],
                 nrs=nr*ns,nst=ns*nt,ntr=nt*nr;
  const unsigned base = 6*nrs + 12*nst + 12*ntr
                        + (d>0 ? 36*nr : 0) + (d>1 ? 36*ns : 0);
  #define DXDN1(i,d)  (fd->sides+base+(9*(i)  +(d))*n)
  #define D2XDN1(i,d) (fd->sides+base+(9*(i)+3+(d))*n)
  #define D2XDN2(i,d) (fd->sides+base+(9*(i)+6+(d))*n)
  const struct findpts_el_3_gface *face_d_n1 = get_face(fd,2*d),
                                  *face_n2_d = get_face(fd,2*dn2);
  struct findpts_el_3_gedge *e = fd->edge + 4*d;
  unsigned i,xd;
  double *work = fd->work;
  for(xd=0;xd<3;++xd) for(i=0;i<2;++i)
    e[2*i  ].x[xd] = face_d_n1[i].x[xd],
    e[2*i+1].x[xd] = face_d_n1[i].x[xd]+n*(dn1-1);
  memcpy(work     , fd->wtend[dn1]+  n1,2*n1*sizeof(double));
  memcpy(work+2*n1, fd->wtend[dn1]+4*n1,2*n1*sizeof(double));
  for(i=0;i<2;++i) for(xd=0;xd<3;++xd) {
    tensor_mxm(work+4*n1,n, face_d_n1[i].x[xd],n1, work,4);
    memcpy( DXDN1(2*i+0,xd), work+4*n1    , n*sizeof(double));
    memcpy(D2XDN1(2*i+0,xd), work+4*n1+  n, n*sizeof(double));
    memcpy( DXDN1(2*i+1,xd), work+4*n1+2*n, n*sizeof(double));
    memcpy(D2XDN1(2*i+1,xd), work+4*n1+3*n, n*sizeof(double));
  }
  memcpy(work   , fd->wtend[dn2]+2*n2,n2*sizeof(double));
  memcpy(work+n2, fd->wtend[dn2]+5*n2,n2*sizeof(double));
  for(i=0;i<2;++i) for(xd=0;xd<3;++xd) {
    tensor_mtxm(work+2*n2,n, face_n2_d[i].x[xd],n2, work,2);
    memcpy(D2XDN2(  i,xd), work+2*n2  , n*sizeof(double));
    memcpy(D2XDN2(2+i,xd), work+2*n2+n, n*sizeof(double));
  }
  #undef D2XDN2
  #undef D2XDN1
  #undef DXDN1
}

static const struct findpts_el_3_gedge *get_edge(
  struct findpts_el_3_data *fd, unsigned ei)
{
  const unsigned mask = 8<<(ei/4);
  if((fd->side_init&mask)==0)
    compute_edge_data(fd,ei/4), fd->side_init |= mask;
  return &fd->edge[ei];
}

/* work[4*nr+4], work[2*nt+2] */
static void compute_pt_data(struct findpts_el_3_data *fd)
{
  const unsigned nr = fd->n[0], nt = fd->n[2];
  const struct findpts_el_3_gedge *e = get_edge(fd,0);
  unsigned d,i;
  double *work = fd->work;
  for(i=0;i<4;++i) for(d=0;d<3;++d)
    fd->pt[2*i  ].x[d] = e[i].x[d][0],
    fd->pt[2*i  ].jac[3*d+1] = e[i].dxdn1[d][0],
    fd->pt[2*i  ].jac[3*d+2] = e[i].dxdn2[d][0],
    fd->pt[2*i  ].hes[6*d+3] = e[i].d2xdn1[d][0],
    fd->pt[2*i  ].hes[6*d+5] = e[i].d2xdn2[d][0],
    fd->pt[2*i+1].x[d] = e[i].x[d][nr-1],
    fd->pt[2*i+1].jac[3*d+1] = e[i].dxdn1[d][nr-1],
    fd->pt[2*i+1].jac[3*d+2] = e[i].dxdn2[d][nr-1],
    fd->pt[2*i+1].hes[6*d+3] = e[i].d2xdn1[d][nr-1],
    fd->pt[2*i+1].hes[6*d+5] = e[i].d2xdn2[d][nr-1];
  memcpy(work     , fd->wtend[0]+  nr, 2*nr*sizeof(double));
  memcpy(work+2*nr, fd->wtend[0]+4*nr, 2*nr*sizeof(double));
  for(i=0;i<4;++i) for(d=0;d<3;++d) {
    tensor_mtxv(work+4*nr,4, work, e[i].x[d],nr);
    fd->pt[2*i  ].jac[3*d  ] = work[4*nr  ];
    fd->pt[2*i  ].hes[6*d  ] = work[4*nr+1];
    fd->pt[2*i+1].jac[3*d  ] = work[4*nr+2];
    fd->pt[2*i+1].hes[6*d  ] = work[4*nr+3];
  }
  memcpy(work+nr,work+2*nr,nr*sizeof(double));
  for(i=0;i<4;++i) for(d=0;d<3;++d) {
    tensor_mtxv(work+2*nr,2, work, e[i].dxdn1[d],nr);
    fd->pt[2*i  ].hes[6*d+1] = work[2*nr  ];
    fd->pt[2*i+1].hes[6*d+1] = work[2*nr+1];
    tensor_mtxv(work+2*nr,2, work, e[i].dxdn2[d],nr);
    fd->pt[2*i  ].hes[6*d+2] = work[2*nr  ];
    fd->pt[2*i+1].hes[6*d+2] = work[2*nr+1];
  }
  e = get_edge(fd,8);
  memcpy(work   , fd->wtend[2]+  nt, nt*sizeof(double));
  memcpy(work+nt, fd->wtend[2]+4*nt, nt*sizeof(double));
  for(i=0;i<4;++i) for(d=0;d<3;++d) {
    tensor_mtxv(work+2*nt,2, work, e[i].dxdn2[d],nt);
    fd->pt[  i].hes[6*d+4] = work[2*nt  ];
    fd->pt[4+i].hes[6*d+4] = work[2*nt+1];
  }
}

static const struct findpts_el_3_gpt *get_pt(
  struct findpts_el_3_data *fd, unsigned pi)
{
  if((fd->side_init&0x40u)==0)
    compute_pt_data(fd), fd->side_init |= 0x40u;
  return &fd->pt[pi];
}

static void newton_vol(struct findpts_el_3_pt *out,
                       const double *jac, const double *resid,
                       const struct findpts_el_3_pt *p, double tol)
{
  double dr[3], r[3]; unsigned d,flags=0;
  lin_solve_3(dr, jac, resid);
  /* check constraints */
  for(d=0;d<3;++d) {
    const double oldr = p->r[d];
    r[d] = oldr + dr[d];
    if(r[d]<-1) r[d]=-1, dr[d]= -1-oldr, flags |= 1<<(2*d);
    else if(r[d]>1) r[d]=1, dr[d]=1-oldr, flags |= 2<<(2*d);
  }
  /* check convergence */
  if(fabs(dr[0])+fabs(dr[1])+fabs(dr[2]) < tol) flags |= 1<<6;
  for(d=0;d<3;++d) out->x[d]=p->x[d];
  for(d=0;d<3;++d) out->r[d]=r[d];
  out->dist2 = resid[0]*resid[0]+resid[1]*resid[1]+resid[2]*resid[2];
  out->flags = flags;
}

static void newton_face(struct findpts_el_3_pt *out,
                        const double *jac, const double *rhes,
                        const double *resid,
                        const unsigned d1, const unsigned d2, const unsigned dn,
                        unsigned flags,
                        const struct findpts_el_3_pt *p, double tol)
{
  /* A = J^T J - resid_d H_d */
  const double A[3] = { jac[  d1]*jac[  d1]
                       +jac[3+d1]*jac[3+d1]
                       +jac[6+d1]*jac[6+d1] - rhes[0],
                        jac[  d1]*jac[  d2]
                       +jac[3+d1]*jac[3+d2]
                       +jac[6+d1]*jac[6+d2] - rhes[1],
                        jac[  d2]*jac[  d2]
                       +jac[3+d2]*jac[3+d2]
                       +jac[6+d2]*jac[6+d2] - rhes[2] };
  /* y = J^T r */
  const double y[2] = { jac[  d1]*resid[0]
                       +jac[3+d1]*resid[1]
                       +jac[6+d1]*resid[2],
                        jac[  d2]*resid[0]
                       +jac[3+d2]*resid[1]
                       +jac[6+d2]*resid[2] };
  double dr[2], r[2];
  unsigned d;
  lin_solve_sym_2(dr,A,y);
  /* check constraints */
  #define UPDATE(d,d3) do { \
    r[d] = p->r[d3] + dr[d]; \
    if(r[d]<-1) r[d]=-1, dr[d]= -1-p->r[d3], flags |= 1<<(2*d3); \
    else if(r[d]>1) r[d]=1, dr[d]=1-p->r[d3], flags |= 2<<(2*d3); \
  } while(0)
  UPDATE(0,d1); UPDATE(1,d2);
  #undef UPDATE
  /* check convergence */
  if(fabs(dr[0])+fabs(dr[1]) < tol) flags |= 1<<6;
  for(d=0;d<3;++d) out->x[d]=p->x[d];
  out->r[dn]=p->r[dn], out->r[d1]=r[0], out->r[d2]=r[1];
  out->dist2 = resid[0]*resid[0]+resid[1]*resid[1]+resid[2]*resid[2];
  out->flags = flags;
}

static void newton_edge(struct findpts_el_3_pt *out,
  const double *jac, const double rhes, const double *resid,
  const unsigned de, const unsigned dn1, const unsigned dn2,
  unsigned flags,
  const struct findpts_el_3_pt *p, double tol)
{
  /* A = J^T J - resid_d H_d */
  const double A = jac[  de]*jac[  de]
                  +jac[3+de]*jac[3+de]
                  +jac[6+de]*jac[6+de] - rhes;
  /* y = J^T r */
  const double y = jac[  de]*resid[0]
                  +jac[3+de]*resid[1]
                  +jac[6+de]*resid[2];
  double dr = y/A;
  /* check constraints */
  double r = p->r[de] + dr;
  unsigned d;
  if     (r<-1) r=-1, dr= -1-p->r[de], flags |= 1<<(2*de);
  else if(r> 1) r= 1, dr=  1-p->r[de], flags |= 2<<(2*de);
  /* check convergence */
  if(fabs(dr) < tol) flags |= 1<<6;
  for(d=0;d<3;++d) out->x[d]=p->x[d];
  out->r[dn1]=p->r[dn1], out->r[dn2]=p->r[dn2], out->r[de]=r;
  out->dist2 = resid[0]*resid[0]+resid[1]*resid[1]+resid[2]*resid[2];
  out->flags = flags;
}

typedef void findpt_fun(struct findpts_el_3_pt *out,
                        struct findpts_el_3_data *fd,
                        const struct findpts_el_3_pt *p, uint pn, double tol);

/* work[(3+9+2*(nr+ns+nt+nrs))*pn + max(2*nr,ns) ] */
static void findpt_vol(struct findpts_el_3_pt *out,
                       struct findpts_el_3_data *fd,
                       const struct findpts_el_3_pt *p, uint pn, double tol)
{
  const unsigned nr=fd->n[0],ns=fd->n[1],nt=fd->n[2],
                 nrs=nr*ns;
  double *resid = fd->work, *jac = resid + 3*pn,
         *wtrs = jac+9*pn, *wtt = wtrs+2*(nr+ns)*pn,
         *slice = wtt+2*nt*pn, *temp = slice + 2*pn*nrs;
  uint i; unsigned d;
  /* evaluate x(r) and jacobian */
  for(i=0;i<pn;++i)
    fd->lag[0](wtrs+2*i*(nr+ns)     , fd->lag_data[0], nr, 1, p[i].r[0]);
  for(i=0;i<pn;++i)
    fd->lag[1](wtrs+2*i*(nr+ns)+2*nr, fd->lag_data[1], ns, 1, p[i].r[1]);
  for(i=0;i<pn;++i)
    fd->lag[2](wtt+2*i*nt           , fd->lag_data[2], nt, 1, p[i].r[2]);
  for(d=0;d<3;++d) {
    tensor_mxm(slice,nrs, fd->x[d],nt, wtt,2*pn);
    for(i=0;i<pn;++i) {
      const double *wtrs_i = wtrs+2*i*(nr+ns), *slice_i = slice+2*i*nrs;
      double *jac_i = jac+9*i+3*d;
      resid[3*i+d] = p[i].x[d] - tensor_ig2(jac_i,
        wtrs_i,nr, wtrs_i+2*nr,ns, slice_i, temp);
      jac_i[2] = tensor_i2(wtrs_i,nr, wtrs_i+2*nr,ns, slice_i+nrs, temp);
    }
  }
  /* perform Newton step */
  for(i=0;i<pn;++i) newton_vol(out+i, jac+9*i, resid+3*i, p+i, tol);
}

/* work[(3+9+3+3*(n1+n2+n1))*pn ] */
static void findpt_face(struct findpts_el_3_pt *out,
                        struct findpts_el_3_data *fd,
                        const struct findpts_el_3_pt *p, uint pn, double tol)
{
  const unsigned pflag = p->flags;
  const unsigned fi = face_index(pflag);
  const unsigned dn = fi>>1, d1 = plus_1_mod_3(dn), d2 = plus_2_mod_3(dn);
  const unsigned n1 = fd->n[d1], n2 = fd->n[d2];
  double *resid = fd->work, *jac = resid+3*pn, *hes = jac+9*pn,
         *wt1 = hes+3*pn, *wt2 = wt1+3*n1*pn,
         *slice = wt2+3*n2*pn;
  const struct findpts_el_3_gface *face = get_face(fd,fi);
  uint i; unsigned d;
  /* evaluate x(r), jacobian, hessian */
  for(i=0;i<pn;++i)
    fd->lag[d1](wt1+3*i*n1, fd->lag_data[d1], n1, 2, p[i].r[d1]);
  for(i=0;i<pn;++i)
    fd->lag[d2](wt2+3*i*n2, fd->lag_data[d2], n2, 2, p[i].r[d2]);
  for(i=0;i<3*pn;++i) hes[i]=0;
  for(d=0;d<3;++d) {
    tensor_mxm(slice,n1, face->x[d],n2, wt2,3*pn);
    for(i=0;i<pn;++i) {
      const double *wt1_i = wt1+3*i*n1, *slice_i = slice+3*i*n1;
      double v[9], r;
      tensor_mtxm(v,3, wt1_i,n1, slice_i,3);
      /* v[3*j + i] = d^i/dr1^i d^j/dr2^j x_d */
      resid[3*i+d] = r = p[i].x[d] - v[0];
      jac[9*i+3*d+d1] = v[1];
      jac[9*i+3*d+d2] = v[3];
      hes[3*i  ] += r * v[2];
      hes[3*i+1] += r * v[4];
      hes[3*i+2] += r * v[6];
    }
  }
  for(i=1;i<pn;++i) memcpy(wt2+i*n2, wt2+3*i*n2, n2*sizeof(double));
  for(d=0;d<3;++d) {
    tensor_mxm(slice,n1, face->dxdn[d],n2, wt2,pn);
    for(i=0;i<pn;++i)
      jac[9*i+3*d+dn] = tensor_dot(wt1+3*i*n1, slice+i*n1, n1);
  }
  /* perform Newton step */
  for(i=0;i<pn;++i) {
    double *resid_i = resid+3*i, *jac_i = jac+9*i, *hes_i = hes+3*i;
    /* check constraint */
    double steep = resid_i[0] * jac_i[  dn]
                  +resid_i[1] * jac_i[3+dn]
                  +resid_i[2] * jac_i[6+dn];
    if(steep * p[i].r[dn] < 0) /* relax constraint */
      newton_vol(out+i, jac_i, resid_i, p+i, tol);
    else
      newton_face(out+i, jac_i, hes_i, resid_i, d1,d2,dn,pflag, p+i, tol);
  }
}

/* work[ 3*n ] */
static void findpt_edge(struct findpts_el_3_pt *out,
  struct findpts_el_3_data *fd,
  const struct findpts_el_3_pt *p, uint pn, double tol)
{
  const unsigned pflag = p->flags;
  const unsigned ei = edge_index(pflag);
  const unsigned de = ei>>2, dn1 = plus_1_mod_3(de), dn2 = plus_2_mod_3(de);
  const unsigned n = fd->n[de];
  double *wt = fd->work;
  const struct findpts_el_3_gedge *edge = get_edge(fd,ei);
  uint i; unsigned d;
  for(i=0;i<pn;++i) {
    double dxi[3], resid[3], jac[9];
    double hes[5] = {0,0,0,0,0};
    /* evaluate x(r), jacobian, hessian */
    fd->lag[de](wt, fd->lag_data[de], n, 2, p[i].r[de]);
    for(d=0;d<3;++d) {
      double r;
      tensor_mtxv(dxi,3, wt, edge->x[d],n);
      resid[d] = r = p[i].x[d] - dxi[0];
      jac[3*d+de] = dxi[1];
      hes[0] += r * dxi[2];
      tensor_mtxv(dxi,2, wt, edge->dxdn1[d],n);
      jac[3*d+dn1] = dxi[0];
      hes[1] += r * dxi[1];
      tensor_mtxv(dxi,2, wt, edge->dxdn2[d],n);
      jac[3*d+dn2] = dxi[0];
      hes[2] += r * dxi[1];
      hes[3] += r * tensor_dot(wt, edge->d2xdn1[d], n);
      hes[4] += r * tensor_dot(wt, edge->d2xdn2[d], n);
    }
    /* check constraint */
    {
      const double steep[3] = {
        jac[0]*resid[0] + jac[3]*resid[1] + jac[6]*resid[2],
        jac[1]*resid[0] + jac[4]*resid[1] + jac[7]*resid[2],
        jac[2]*resid[0] + jac[5]*resid[1] + jac[8]*resid[2] };
      const double sr1 = steep[dn1]*p[i].r[dn1],
                   sr2 = steep[dn2]*p[i].r[dn2];
      if(sr1<0) {
        if(sr2<0)
          newton_vol(out+i, jac,resid, p+i, tol);
        else {
          const double rh[3] = { hes[0], hes[1], hes[3] };
          newton_face(out+i, jac,rh,resid, de,dn1,dn2,
                      pflag & (3u<<(dn2*2)), p+i, tol);
        }
      } else if(sr2<0) {
          const double rh[3] = { hes[4], hes[2], hes[0] };
          newton_face(out+i, jac,rh,resid, dn2,de,dn1,
                      pflag & (3u<<(dn1*2)), p+i, tol);
      } else
        newton_edge(out+i, jac,hes[0],resid, de,dn1,dn2, pflag, p+i, tol);
    }
  }
}

static void findpt_pt(struct findpts_el_3_pt *out,
  struct findpts_el_3_data *fd,
  const struct findpts_el_3_pt *p, uint pn, double tol)
{
  const unsigned pflag = p->flags;
  const unsigned pi = point_index(pflag);
  const struct findpts_el_3_gpt *gpt = get_pt(fd,pi);
  const double *x = gpt->x, *jac = gpt->jac, *hes = gpt->hes;
  uint i;
  for(i=0;i<pn;++i) {
    const double resid[3] = { p[i].x[0]-x[0],
                              p[i].x[1]-x[1],
                              p[i].x[2]-x[2] };
    const double steep[3] = {
      jac[0]*resid[0] + jac[3]*resid[1] + jac[6]*resid[2],
      jac[1]*resid[0] + jac[4]*resid[1] + jac[7]*resid[2],
      jac[2]*resid[0] + jac[5]*resid[1] + jac[8]*resid[2] };
    const double sr[3] = { steep[0]*p[i].r[0],
                           steep[1]*p[i].r[1], 
                           steep[2]*p[i].r[2] };
    unsigned d1,d2,dn, de,dn1,dn2, hi0,hi1,hi2;
    /* check constraints */
    if(sr[0]<0) {
      if(sr[1]<0) {
        if(sr[2]<0) goto findpt_pt_vol;
        else { d1=0,d2=1,dn=2, hi0=0,hi1=1,hi2=3; goto findpt_pt_face; }
      }
      else if(sr[2]<0) {d1=2,d2=0,dn=1, hi0=5,hi1=2,hi2=0; goto findpt_pt_face;}
      else { de=0,dn1=1,dn2=2, hi0=0; goto findpt_pt_edge; }
    }
    else if(sr[1]<0) {
      if(sr[2]<0) { d1=1,d2=2,dn=0, hi0=3,hi1=4,hi2=5; goto findpt_pt_face; }
      else { de=1,dn1=2,dn2=0, hi0=3; goto findpt_pt_edge; }
    }
    else if(sr[2]<0) { de=2,dn1=0,dn2=1, hi0=5; goto findpt_pt_edge; }
    out[i] = p[i], out[i].flags = pflag | (1<<6);
    out->dist2 = resid[0]*resid[0]+resid[1]*resid[1]+resid[2]*resid[2];
    continue;
    findpt_pt_vol:
      newton_vol(out+i, jac,resid, p+i, tol);
      continue;
    findpt_pt_face: {
      const double rh[3] = {
        resid[0]*hes[hi0]+resid[1]*hes[6+hi0]+resid[2]*hes[12+hi0],
        resid[0]*hes[hi1]+resid[1]*hes[6+hi1]+resid[2]*hes[12+hi1],
        resid[0]*hes[hi2]+resid[1]*hes[6+hi2]+resid[2]*hes[12+hi2] };
      newton_face(out+i, jac,rh,resid, d1,d2,dn,
                  pflag&(3u<<(2*dn)), p+i, tol);
    } continue;
    findpt_pt_edge: {
      const double rh =
        resid[0]*hes[hi0]+resid[1]*hes[6+hi0]+resid[2]*hes[12+hi0];
      newton_edge(out+i, jac,rh,resid, de,dn1,dn2,
                  pflag&~(3u<<(2*de)), p+i, tol);
    } continue;
  }
}

void findpts_el_3(struct findpts_el_3_data *fd, uint npt, const double tol)
{
  findpt_fun *const fun[4] = 
    { &findpt_vol, &findpt_face, &findpt_edge, &findpt_pt };
  struct findpts_el_3_pt *const pstart = fd->p, *const pbuf = fd->p + npt;
  unsigned step = 0;
  uint count[27] = { npt,0,0, 0,0,0, 0,0,0,
                       0,0,0, 0,0,0, 0,0,0,
                       0,0,0, 0,0,0, 0,0,0 } ;
  fd->p->flags = 0;
  while(npt && step++ < 50) {
    /* advance each group of points */
    struct findpts_el_3_pt *p, *const pe=pstart+npt, *pout; uint pn;
    for(p=pstart,pout=pbuf; p!=pe; p+=pn,pout+=pn) {
      pn = count[pt_flags_to_bin_noC(p->flags)];
      fun[num_constrained(p->flags)](pout, fd, p,pn, tol);
    }
    /* group points by contsraints */
    {
      uint offset[28] = { 0,0,0, 0,0,0, 0,0,0,
                          0,0,0, 0,0,0, 0,0,0,
                          0,0,0, 0,0,0, 0,0,0, 0 };
      struct findpts_el_3_pt *const pe = pbuf+npt;
      for(pout=pbuf; pout!=pe; ++pout)
        ++offset[pt_flags_to_bin(pout->flags)];
      {
        unsigned i; uint sum=0;
        for(i=0;i<27;++i) {
          uint ci=offset[i]; count[i]=ci, offset[i]=sum, sum+=ci;
        }
        npt = offset[27] = sum; /* last bin is converged; forget it */
      }
      for(pout=pbuf; pout!=pe; ++pout)
        pstart[offset[pt_flags_to_bin(pout->flags)]++] = *pout;
    }
  }
}

