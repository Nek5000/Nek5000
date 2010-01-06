#include <stdio.h>

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "tensor.h"
#include "poly.h"

#define findpts_el_setup_2   PREFIXED_NAME(findpts_el_setup_2)
#define findpts_el_free_2    PREFIXED_NAME(findpts_el_free_2 )
#define findpts_el_2         PREFIXED_NAME(findpts_el_2      )
#define findpts_el_eval_2    PREFIXED_NAME(findpts_el_eval_2 )

#define DIAGNOSTICS_ITERATIONS 0

struct findpts_el_pt_2 {
  double x[2],r[2],oldr[2],dist2,dist2p,tr;
  unsigned index,flags;
};

/* the bit structure of flags is CSSRR
   the C bit --- 1<<4 --- is set when the point is converged
   RR is 0 = 00b if r is unconstrained,
         1 = 01b if r is constrained at -1
         2 = 10b if r is constrained at +1
   SS is similarly for s constraints
*/

#define CONVERGED_FLAG (1u<<4)
#define FLAG_MASK 0x1fu

static unsigned num_constrained(const unsigned flags)
{
  const unsigned y = flags | flags>>1;
  return (y&1u) + (y>>2 & 1u);
}

static unsigned pt_flags_to_bin_noC(const unsigned flags)
{
  return (flags>>2 & 3u)*3 + (flags & 3u);
}

/* map flags to 9 if the C bit is set,
   else to [0,8] --- the 9 valid configs of SSRR */
static unsigned pt_flags_to_bin(const unsigned flags)
{
  const unsigned mask = 0u - (flags>>4); /* 0 or 0xfff... when converged */
  return (mask & 9u) | (~mask & pt_flags_to_bin_noC(flags));
}

/* assumes x = 0, or 1  */
static unsigned plus_1_mod_2(const unsigned x) { return x^1u; } 

/* assumes x = 1 << i, with i < 4, returns i+1 */
static unsigned which_bit(const unsigned x)
{
  const unsigned y = x&7u;
  return (y-(y>>2)) | ((x-1)&4u);
}

static unsigned edge_index(const unsigned x) { return which_bit(x)-1; }

static unsigned point_index(const unsigned x)
{
  return ((x>>1)&1u) | ((x>>2)&2u);
}

/* extra data

  we need x, dx/dn for each edge
    r: x at 0, nrs - nr,
      4*nr extra for dx/dn
    s: 8*ns extra

*/

struct findpts_el_gedge_2 { const double *x[2], *dxdn[2]; };
struct findpts_el_gpt_2   { double x[2], jac[4], hes[6]; };

struct findpts_el_data_2 {
  unsigned npt_max;
  struct findpts_el_pt_2 *p;

  unsigned n[2];
  double *z[2];
  lagrange_fun *lag[2];
  double *lag_data[2];
  double *wtend[2];
  
  const double *x[2];
  
  unsigned side_init;
  double *sides;
  struct findpts_el_gedge_2 edge[4]; /* R S=-1; R S=1; ... */
  struct findpts_el_gpt_2 pt[4];

  double *work;
};

static unsigned work_size(
  const unsigned nr, const unsigned ns, const unsigned npt_max)
{
  fail(1,"findpts_el_2 not yet implemented");
}

void findpts_el_setup_2(struct findpts_el_data_2 *fd,
                        const unsigned n[2],
                        const unsigned npt_max)
{
  const unsigned nr=n[0], ns=n[1];
  const unsigned tot = 8*ns + 4*nr;
  unsigned d,i, lag_size[2];

  fd->npt_max = npt_max;
  fd->p = tmalloc(struct findpts_el_pt_2, npt_max*2);

  fd->n[0]=nr, fd->n[1]=ns;
  for(d=0;d<2;++d) lag_size[d] = gll_lag_size(fd->n[d]);

  fd->z[0]        = tmalloc(double,lag_size[0]+lag_size[1]
                                   +7*(nr+ns) + tot +
                                   work_size(nr,ns,npt_max));
  fd->z[1]        = fd->z[0]+nr;
  fd->lag_data[0] = fd->z[1]+ns;
  fd->lag_data[1] = fd->lag_data[0]+lag_size[0];
  fd->wtend[0]    = fd->lag_data[1]+lag_size[1];
  fd->wtend[1]    = fd->wtend[0]+6*nr;
  fd->sides       = fd->wtend[1]+6*ns;
  fd->work        = fd->sides + tot;

  fd->side_init = 0;
  
  for(d=0;d<2;++d) {
    double *wt=fd->wtend[d]; unsigned n=fd->n[d];
    lobatto_nodes(fd->z[d],n);
    fd->lag[d] = gll_lag_setup(fd->lag_data[d],n);
    fd->lag[d](wt    , fd->lag_data[d],n,2,-1);
    fd->lag[d](wt+3*n, fd->lag_data[d],n,2, 1);
    
    wt[0]=1; for(i=1;i<n;++i) wt[i]=0;
    wt+=3*n; { for(i=0;i<n-1;++i) wt[i]=0; } wt[i]=1;
  }

  for(d=0;d<2;++d)
    fd->edge[0].x[d]    = fd->sides +    d *ns, \
    fd->edge[0].dxdn[d] = fd->sides + (2+d)*ns, \
    fd->edge[1].x[d]    = fd->sides + (4+d)*ns, \
    fd->edge[1].dxdn[d] = fd->sides + (6+d)*ns; \

  for(d=0;d<2;++d)
    fd->edge[2].x[d] = 0, /* will point to user data */
    fd->edge[2].dxdn[d] = fd->sides + 8*ns +    d *nr,
    fd->edge[3].x[d] = 0, /* will point to user data */
    fd->edge[3].dxdn[d] = fd->sides + 8*ns + (2+d)*nr;
}

void findpts_el_free_2(struct findpts_el_data_2 *fd)
{
  free(fd->p);
  free(fd->z[0]);
}

typedef void compute_edge_data_fun(struct findpts_el_data_2 *fd);

/* work[2*(nr+ns)] */
static void compute_edge_data_r(struct findpts_el_data_2 *fd)
{
  const unsigned nr = fd->n[0], ns=fd->n[1], nrsm1 = nr*(ns-1);
  unsigned d;
  double *work = fd->work, *out = fd->sides + 8*ns;
  memcpy(work   , fd->wtend[1]+  ns, ns*sizeof(double));
  memcpy(work+ns, fd->wtend[1]+4*ns, ns*sizeof(double));
  for(d=0;d<2;++d) {
    tensor_mxm(work+2*ns,nr, fd->x[d],ns, work,2);
    memcpy(out+   d *nr, work+2*ns      , nr*sizeof(double));
    memcpy(out+(2+d)*nr, work+2*ns+nr   , nr*sizeof(double));
    fd->edge[2].x[d] = fd->x[d];
    fd->edge[3].x[d] = fd->x[d] + nrsm1;
  }
}

/* work[4*(nr+ns)] */
static void compute_edge_data_s(struct findpts_el_data_2 *fd)
{
  const unsigned nr = fd->n[0], ns=fd->n[1];
  unsigned d;
  double *work = fd->work, *out = fd->sides;
  memcpy(work     , fd->wtend[0]     , 2*nr*sizeof(double));
  memcpy(work+2*nr, fd->wtend[0]+3*nr, 2*nr*sizeof(double));
  for(d=0;d<2;++d) {
    tensor_mtxm(work+4*nr,ns, fd->x[d],nr, work,4);
    memcpy(out+   d *ns, work+4*nr     , ns*sizeof(double));
    memcpy(out+(2+d)*ns, work+4*nr+  ns, ns*sizeof(double));
    memcpy(out+(4+d)*ns, work+4*nr+2*ns, ns*sizeof(double));
    memcpy(out+(6+d)*ns, work+4*nr+3*ns, ns*sizeof(double));
  }
}

static const struct findpts_el_gedge_2 *get_edge(
  struct findpts_el_data_2 *fd, unsigned ei)
{
  const unsigned mask = 1<<(ei/2);
  if((fd->side_init&mask)==0) {
    compute_edge_data_fun *const fun[2] = {
      compute_edge_data_s,
      compute_edge_data_r
    };
    fun[ei/2](fd);
    fd->side_init |= mask;
  }
  return &fd->edge[ei];
}

void findpts_el_2(struct findpts_el_data_2 *fd, unsigned npt, const double tol)
{
  fail(1,"findpts_el_2 not yet implemented");
}

void findpts_el_eval_2(
        double *const out_base, const unsigned out_stride,
  const double *const   r_base, const unsigned   r_stride, const unsigned pn,
  const double *const in, struct findpts_el_data_2 *fd)
{
  fail(1,"findpts_el_eval_2 not yet implemented");
}
