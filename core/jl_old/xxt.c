#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#ifdef MPI
#  include <string.h>
#  include <math.h>
#  include <mpi.h>
#endif

#include "errmem.h"
#include "types.h"
#include "sort.h"
#include "sparse_cholesky.h"
#ifdef MPI
#  include "minmax.h"
#  include "crystal.h"
#  include "gs.h"
#endif
#include "tuple_list.h"

#define xxt_free xxt_jl_free
#define xxt_solve xxt_jl_solve
#define xxt_stats xxt_jl_stats

/*
  portable log base 2
  
  does a binary search to find leading order bit
  
  UINT_BITS = number of bits in a uint
  BITS(0) = UINT_BITS
  BITS(i) = half of BITS(i-1), rounded up
  MASK(i) = bitmask with BITS(i) 1's followed by BITS(i) 0's
*/

static unsigned lg(uint v)
{
  unsigned r = 0;
#define UINT_BITS (sizeof(uint)*CHAR_BIT)
#define BITS(i) ((UINT_BITS+(1<<(i))-1)>>(i))
#define MASK(i) ((((uint)1<<BITS(i)) - 1) << BITS(i))
#define CHECK(i) if((BITS(i)!=1) && (v&MASK(i))) v>>=BITS(i), r+=BITS(i)  
  CHECK(1); CHECK(2); CHECK(3); CHECK(4); CHECK(5); CHECK(6); CHECK(7);
  CHECK(8); CHECK(9); /* this covers up to 1024-bit uints */
  if(v&2) ++r;
  return r;
#undef UINT_BITS
#undef BITS
#undef MASK
#undef CHECK
}

typedef struct {
  uint n, *Arp, *Aj; real *A;
} csr_mat;

#ifdef MPI

typedef struct {

  /* communication */

  MPI_Comm comm;
  uint pid, np;  /* proc id, number of procs */
  uint pcoord;   /* coordinate in communication tree */ 
  unsigned plevels; /* # of stages of communication */
  sint *pother;     /* let p = pother[i], then during stage i of fan-in,
                           if p>=0, receive from p
                           if p< 0, send to (-p-1)
                       fan-out is just the reverse ...
                       on proc 0, pother is never negative
                       on others, pother is negative for the last stage only */
  MPI_Request *mpireq;
  MPI_Status  *mpistatus;
  
  /* separators */

  unsigned nsep;  /* number of separators */
  uint *sep_size; /* # of dofs on each separator,
                     ordered from the bottom to the top of the tree:
                     separator 0 is the bottom-most one (dofs not shared)
                     separator nsep-1 is the root of the tree */

  unsigned null_space;
  real *share_weight;

  /* vector sizes */

  uint un;        /* user's vector size */
  
  /* xxt_solve works with "condensed" vectors;
     same dofs as in user's vectors, but no duplicates and no Dirichlet nodes,
     and also ordered topologically (children before parents) according to the
     separator tree */
  
  uint cn;        /* size of condensed vectors */
  sint *perm_u2c; /* permutation from user vector to condensed vector,
                     p=perm_u2c[i]; xu[i] = p=-1 ? 0 : xc[p];          */
  uint ln, sn;    /* xc[0 ... ln-1] are not shared   (ln=sep_size[0])
                     xc[ln ... ln+sn-1] are shared
                     ln+sn = cn                    */
  
  uint xn;        /* # of columns of x = sum_i(sep_size[i]) - sep_size[0] */

  /* data */
  sparse_cholesky_data fac_A_ll;
  csr_mat              A_sl;
  uint *Xp; real *X;   /* column i of X starts at X[Xp[i]] */
  
  /* execution buffers */
  real *vl, *vc, *vx, *combuf;
} xxt_data;

#else

typedef struct {
  unsigned null_space;
  uint un, cn;
  sint *perm_u2c;
  sparse_cholesky_data fac_A_ll;
  real *vc;
} xxt_data;

#endif

#ifdef MPI

/*
  for the binary communication tree, the procs are divided in half
  at each level, with the second half always the larger
  
  e.g., for np = 13:

       +------13-------+
       |               |
   +---6---+       +---7---+
   |       |       |       |
 +-3-+   +-3-+   +-3-+   +-4-+
 1   2   1   2   1   2   2   2
    1^1     1^1     1^1 1^1 1^1

  plevels is the number of levels in the tree
    = np==1 ? 1 : ( lg(np-1)+2 )

  labelling the nodes with proc id's gives this communication tree:

       +-------0-------+
       |               |
   +---0---+       +---6---+
   |       |       |       |
 +-0-+   +-3-+   +-6-+   +-9-+
 0   1   3   4   6   7   9   b
    1^2     4^5     7^8 9^a b^c

  consider proc 7 (pid = 7);
  pcoord gives the position of the leaf labelled 7:
    Root Right Left Right Left -> RRLRL -> 11010
    so pcoord = 11010 binary
  note the parent coordinate can be found by bit shifting right
    (i.e. dividing by 2)
*/

/* sets: pcoord, nsep, plevels, pother, mpireq, mpistatus */
static void locate_proc(xxt_data *data)
{
  const uint id = data->pid;
  uint n = data->np, c=1, odd=0, base=0;
  unsigned level=0;
  while(n>1) {
    ++level;
    odd=(odd<<1)|(n&1);
    c<<=1, n>>=1;
    if(id>=base+n) c|=1, base+=n, n+=(odd&1);
  }
  data->pcoord=c;
  data->nsep = level+1;
  data->plevels = data->nsep-1;
  data->pother = tmalloc(sint,data->plevels);
  data->mpireq    = tmalloc(MPI_Request,data->plevels);
  data->mpistatus = tmalloc(MPI_Status ,data->plevels);
  for(level=0;level<data->plevels;++level) {
    if((c&1)==1) {
      uint targ = id - (n-(odd&1));
      data->pother[level]=-(sint)(targ+1);
      data->plevels = level+1;
      break;
    } else {
      data->pother[level]=id+n;
      c>>=1, n=(n<<1)+(odd&1), odd>>=1;
    }
  }
}

/* the tuple list describing the condensed dofs:
   [(separator level, share count, global id)] */
static const unsigned dof_mi=2, dof_ml=1;
static const unsigned dof_level=0, dof_count=1, dof_id=2;

/* determine the size of each separator;
   sums the separator sizes following the fan-in, fan-out comm. pattern
   uses the share-counts to avoid counting dofs more than once */
/* sets: xn, sep_size, ln, sn */
static void discover_sep_sizes(xxt_data *data, tuple_list *dof, buffer *buf)
{
  const unsigned ns=data->nsep, nl=data->plevels;
  const uint n = dof->n;
  float *v, *recv;
  unsigned i,lvl; uint j;
  MPI_Status status;

  buffer_reserve(buf,2*ns*sizeof(float));
  v=buf->ptr, recv=v+ns;
    
  for(i=0;i<ns;++i) v[i]=0;
  for(j=0;j<n;++j) v[dof->vi[dof_mi*j+dof_level]]
                      +=1/(float)dof->vi[dof_mi*j+dof_count];

  /* fan-in */
  for(lvl=0;lvl<nl;++lvl) {
    sint other = data->pother[lvl];
    unsigned s = ns-(lvl+1);
    if(other<0) {
      MPI_Send(v   +lvl+1,s,MPI_FLOAT,-other-1,s,data->comm);
    } else {
        MPI_Recv(recv+lvl+1,s,MPI_FLOAT,other,s,data->comm,&status);
        for(i=lvl+1;i<ns;++i) v[i]+=recv[i];
    }
  }
  /* fan-out */
  for(;lvl;) {
    sint other = data->pother[--lvl];
    unsigned s = ns-(lvl+1);
    if(other<0)
      MPI_Recv(v+lvl+1,s,MPI_FLOAT,-other-1,s,data->comm,&status);
    else
      MPI_Send(v+lvl+1,s,MPI_FLOAT,other,s,data->comm);
  }

  data->xn=0;    
  data->sep_size = tmalloc(uint,ns);
  for(i=0;i<ns;++i) { 
    uint s=v[i]+.1f;
    data->sep_size[i]=s;
    data->xn+=s;
  }
  data->ln=data->sep_size[0];
  data->sn=data->cn-data->ln;
  data->xn-=data->ln;
}

/* assuming [A,Aend) is sorted,
   removes 0's and any duplicate entries,
   returns new end */
static ulong *unique_nonzero(ulong *A, ulong *Aend)
{
  if(Aend==A) return A;
  else {
    ulong *end = Aend-1, last=*end, *p=A,*q=A,v=0;
    *end = 1;
    while(*q==0) ++q;   /*  *q==0 => q!=end since *end==0 */
    *end = 0;
    while(q!=end) {
      v=*q++, *p++=v; 
      while(*q==v) ++q; /*  *q==v => q!=end since *end==0 */
    }
    if(last!=v) *p++=last;
    return p;
  }
}

static void merge_sep_ids(xxt_data *data, ulong *sep_id, ulong *other,
                          ulong *work, unsigned s0)
{
  const unsigned ns = data->nsep;
  unsigned s;
  ulong *p=sep_id, *q=other;
  for(s=s0;s<ns;++s) {
    ulong *end;
    uint size = data->sep_size[s];
    memcpy(work+2*size,p,size*sizeof(ulong));
    memcpy(work+3*size,q,size*sizeof(ulong));
    sort_long(work+2*size,2*size,1, work, work+4*size);
    end = unique_nonzero(work,work+2*size);
    memcpy(p,work,(end-work)*sizeof(ulong));
    p+=size, q+=size;
  }
}

static void init_sep_ids(xxt_data *data, tuple_list *dof, ulong *xid)
{
  const unsigned ns=data->nsep;
  const uint n=data->cn, *sep_size=data->sep_size;
  unsigned s=1;
  uint i, size;
  if(ns==1) return;
  size=sep_size[s];
  for(i=data->ln;i<n;++i) {
    unsigned si = dof->vi[dof_mi*i+dof_level];
    while(s!=si) {
      memset(xid,0,size*sizeof(ulong));
      xid+=size;
      size=data->sep_size[++s];
    }
    *xid++ = dof->vl[i], --size;
  }
  while(s!=ns) {
    memset(xid,0,size*sizeof(ulong));
    xid+=size;
    size=data->sep_size[++s];
  }
}

static void find_perm_x2c(uint ln, uint cn, const ulong *cid,
                          uint xn, const ulong *xid, sint *perm)
{
  const ulong *cid_end = cid+cn, *xid_end = xid+xn; uint i=ln;
  cid+=ln;
  while(cid!=cid_end) {
    ulong v=*cid;
    while(*xid!=v) ++xid, *perm++ = -1;
    *perm++ = i++, ++cid, ++xid;
  }
  while(xid!=xid_end) ++xid, *perm++ = -1;
}

/* sets: perm_x2c */
static sint* discover_sep_ids(xxt_data *data, tuple_list *dof, buffer *buf)
{
  const unsigned ns=data->nsep, nl=data->plevels;
  const uint xn=data->xn, *sep_size=data->sep_size;
  ulong *xid, *recv, *work, *p;
  unsigned lvl;
  uint size,ss;
  MPI_Status status;
  sint *perm_x2c;
  
  size=0; for(lvl=1;lvl<ns;++lvl) if(sep_size[lvl]>size) size=sep_size[lvl];
  buffer_reserve(buf,(2*xn+6*size)*sizeof(ulong));
  xid=buf->ptr, recv=xid+xn, work=recv+xn;
  
  init_sep_ids(data,dof,xid);

  if(nl) {
    /* fan-in */
    p=xid, size=xn;
    for(lvl=0;lvl<nl;++lvl) {
      sint other = data->pother[lvl];
      if(other<0) {
        MPI_Send(p   ,size*sizeof(ulong),MPI_UNSIGNED_CHAR,
                 -other-1,size,data->comm);
      } else {
        MPI_Recv(recv,size*sizeof(ulong),MPI_UNSIGNED_CHAR,
                 other,size,data->comm,&status);
        merge_sep_ids(data,p,recv,work,lvl+1);
      }
      ss=data->sep_size[lvl+1];
      if(ss>=size || lvl==nl-1) break;
      p+=ss, size-=ss;
    }
    /* fan-out */
    for(;;) {
      sint other = data->pother[lvl];
      if(other<0)
        MPI_Recv(p,size*sizeof(ulong),MPI_UNSIGNED_CHAR,
                 -other-1,size,data->comm,&status);
      else
        MPI_Send(p,size*sizeof(ulong),MPI_UNSIGNED_CHAR,
                 other,size,data->comm);
      if(lvl==0) break;
      ss=data->sep_size[lvl];
      p-=ss, size+=ss, --lvl;
    }
  }

#if 0
  printf("xid%d:",data->pid);
  { uint i; for(i=0;i<xn;++i) printf(" %d",(int)xid[i]); }
  printf("\n");
#endif
  
  perm_x2c=tmalloc(sint,xn);
  find_perm_x2c(data->ln,data->cn,(ulong*)dof->vl, xn,xid, perm_x2c);

#if 0
  printf("x2c %d:",data->pid);
  { uint i; for(i=0;i<xn;++i) printf(" %d",(int)perm_x2c[i]); }
  printf("\n");
#endif

  return perm_x2c;
}

static void apply_QQt(xxt_data *data, real *v, uint n, uint tag)
{
  const unsigned nl=data->plevels;
  real *p=v, *recv=data->combuf;
  unsigned lvl, nsend=0;
  uint size=n, ss;
  MPI_Status status;
  
  if(n==0 || nl==0) return;

  tag=tag*2+0;
  /* fan-in */
  for(lvl=0;lvl<nl;++lvl) {
    sint other = data->pother[lvl];
    if(other<0) {
      MPI_Send(p   ,size*sizeof(real),MPI_UNSIGNED_CHAR,
               -other-1,tag,data->comm);
    } else {
      uint i;
      MPI_Recv(recv,size*sizeof(real),MPI_UNSIGNED_CHAR,
               other   ,tag,data->comm,&status);
      for(i=0;i<size;++i) p[i]+=recv[i];
    }
    ss=data->sep_size[lvl+1];
    if(ss>=size || lvl==nl-1) break;
    p+=ss, size-=ss;
  }
  /* fan-out */
  for(;;) {
    sint other = data->pother[lvl];
    if(other<0) {
      MPI_Recv (p,size*sizeof(real),MPI_UNSIGNED_CHAR,
                -other-1,tag,data->comm,&status);
    } else {
      MPI_Isend(p,size*sizeof(real),MPI_UNSIGNED_CHAR,
                other   ,tag,data->comm,&data->mpireq[nsend++]);
    }
    if(lvl==0) break;
    ss=data->sep_size[lvl];
    p-=ss, size+=ss, --lvl;
  }
  if(nsend) MPI_Waitall(nsend,data->mpireq,data->mpistatus);
}

static real sum(xxt_data *data, real v, uint n, uint tag)
{
  const unsigned nl=data->plevels;
  real r;
  unsigned lvl,nsend=0;
  uint size=n, ss;
  MPI_Status status;

  tag=tag*2+1;
  if(n==0 || nl==0) return v;
  /* fan-in */
  for(lvl=0;lvl<nl;++lvl) {
    sint other = data->pother[lvl];
    if(other<0) {
      MPI_Send(&v,sizeof(real),MPI_UNSIGNED_CHAR,
               -other-1,tag,data->comm);
    } else {
      MPI_Recv(&r,sizeof(real),MPI_UNSIGNED_CHAR,
               other   ,tag,data->comm,&status);
      v+=r;
    }
    ss=data->sep_size[lvl+1];
    if(ss>=size || lvl==nl-1) break;
    size-=ss;
  }
  /* fan-out */
  for(;;) {
    sint other = data->pother[lvl];
    if(other<0) {
      MPI_Recv (&v,sizeof(real),MPI_UNSIGNED_CHAR,
                -other-1,tag,data->comm,&status);
    } else {
      MPI_Isend(&v,sizeof(real),MPI_UNSIGNED_CHAR,
                other   ,tag,data->comm,&data->mpireq[nsend++]);
    }
    if(lvl==0) break;
    ss=data->sep_size[lvl];
    size+=ss, --lvl;
  }
  if(nsend) MPI_Waitall(nsend,data->mpireq,data->mpistatus);
  return v;
}

#endif

/* sorts an array of ids, removes 0's and duplicates;
   just returns the permutation */
static uint unique_ids(uint n, const ulong *id, sint *perm, buffer *buf)
{
  uint *p, i, un=0; ulong last=0;
  buffer_reserve(buf,2*n*sizeof(sort_data_long));
  p = buf->ptr;
  index_sort_long(id,n,1, p, buf->ptr);
  for(i=0;i<n;++i) {
    uint j = p[i]; ulong v = id[j];
    if(v==0) perm[j]=-1;
    else {
      if(v!=last) last=v, ++un;
      perm[j]=un-1;
    }
  }
  return un;
}

#ifdef MPI

/* given user's list of dofs (as id's)
   uses gather-scatter to find share-count and separator # for each
   outputs as a list, sorted topologically (children before parents)
                      according to the sep. tree (and without duplicates),
           as well as the permutation to get there from the user's list */
/* sets: un, cn, perm_u2c */
static void discover_dofs(xxt_data *data, uint n, const ulong *id,
                          tuple_list *dof, crystal_data *crystal)
{
  const uint pcoord = data->pcoord, ns=data->nsep;
  sint *perm;
  uint i, cn, *p, *pi;
  gs_data *gs; real *v;
  buffer *buf=&crystal->all->buf;
  
  data->un = n;
  data->perm_u2c = perm = tmalloc(sint,n);
  data->cn = cn = unique_ids(n,id,perm,buf);
  tuple_list_init_max(dof,dof_mi,dof_ml,0,cn); dof->n=cn;
  for(i=0;i<n;++i) if(perm[i]>=0) dof->vl[perm[i]]=id[i];

  gs = gs_data_setup(cn,(ulong*)dof->vl,1,crystal);
  
  buf = &crystal->all->buf;
  buffer_reserve(buf,cn*sizeof(real));
  v = buf->ptr;

  for(i=0;i<cn;++i) v[i]=pcoord;
  gs_op(v,GS_OP_BPR,gs);
  for(i=0;i<cn;++i) dof->vi[dof_mi*i+dof_level]=ns-1-lg((uint)v[i]);
  
  for(i=0;i<cn;++i) v[i]=1;
  gs_op(v,GS_OP_ADD,gs);
  for(i=0;i<cn;++i) dof->vi[dof_mi*i+dof_count]=v[i]+.1;
  
  gs_data_free(gs);

  buf = &crystal->all->buf;
  buffer_reserve(buf,umax_3(2*cn*sizeof(sort_data),
                            cn*sizeof(uint)+cn*dof_mi*sizeof(sint),
                            cn*sizeof(uint)+cn*dof_ml*sizeof(slong)));
  p=buf->ptr;
  index_sort((uint*)&dof->vi[dof_level],cn,dof_mi, p, buf->ptr);
  tuple_list_permute(dof,p,p+cn);
  pi = p+cn; for(i=0;i<cn;++i) pi[p[i]]=i;
  for(i=0;i<n;++i) if(perm[i]>=0) perm[i]=pi[perm[i]];

#if 0
  printf("id    %d:",crystal->id);
  for(i=0;i<n;++i) printf(" %d", id[i]);
  printf("\n");
  printf("perm  %d:",crystal->id);
  for(i=0;i<n;++i) printf(" %d", perm[i]);
  printf("\n");
  printf("id    %d:",crystal->id);
  for(i=0;i<cn;++i) printf(" %d", dof->vl[i]);
  printf("\n");
  printf("level %d:",crystal->id);
  for(i=0;i<cn;++i) printf(" %d", dof->vi[dof_mi*i+dof_level]);
  printf("\n");
  printf("count %d:",crystal->id);
  for(i=0;i<cn;++i) printf(" %d", dof->vi[dof_mi*i+dof_count]);
  printf("\n");
#endif
}

static real inner(const real *u, const real *v, unsigned n)
{
  const real *u_end = u+n;
  real sum = 0;
  while(u!=u_end) { sum += *u++ * *v++; }
  return sum;
}

/* vl += A_ls * vs */
static void apply_p_Als(real *vl, xxt_data *data, const real *vs, uint ns)
{
  const uint *Arp = data->A_sl.Arp,
             *Aj  = data->A_sl.Aj;
  const real *A   = data->A_sl.A;
  uint i,p,pe;
  for(i=0;i<ns;++i)
    for(p=Arp[i],pe=Arp[i+1];p!=pe;++p)
      vl[Aj[p]]+=A[p]*vs[i];
}

/* vs -= A_sl * vl */
static void apply_m_Asl(real *vs, uint ns, xxt_data *data, const real *vl)
{
  const uint *Arp = data->A_sl.Arp,
             *Aj  = data->A_sl.Aj;
  const real *A   = data->A_sl.A;
  uint i,p,pe;
  for(i=0;i<ns;++i)
    for(p=Arp[i],pe=Arp[i+1];p!=pe;++p)
      vs[i]-=A[p]*vl[Aj[p]];
}

/* returns a column of S : vs = -S(0:ei-1,ei) */
static void apply_S_col(real *vs, xxt_data *data, csr_mat *A_ss, uint ei,
                        real *vl)
{
  const uint ln=data->ln;
  const uint *Asl_rp = data->A_sl.Arp, *Ass_rp = A_ss->Arp,
             *Asl_j  = data->A_sl.Aj,  *Ass_j  = A_ss->Aj;
  const real *Asl    = data->A_sl.A,   *Ass    = A_ss->A;
  uint i,p,pe;
  for(i=0;i<ei;++i) vs[i]=0;
  for(p=Ass_rp[ei],pe=Ass_rp[ei+1];p!=pe;++p) {
    uint j=Ass_j[p];
    if(j>=ei) break;
    vs[j]=-Ass[p];
  }
  for(i=0;i<ln;++i) vl[i]=0;
  for(p=Asl_rp[ei],pe=Asl_rp[ei+1];p!=pe;++p) vl[Asl_j[p]]=-Asl[p];
  sparse_cholesky_solve(vl,&data->fac_A_ll,vl);
  apply_m_Asl(vs,ei,data,vl);
}

static void apply_S(real *Svs, uint ns, xxt_data *data, csr_mat *A_ss,
                    const real *vs, real *vl)
{
  const uint ln=data->ln;
  const uint *Ass_rp = A_ss->Arp,
             *Ass_j  = A_ss->Aj;
  const real *Ass    = A_ss->A;
  uint i, p,pe;
  for(i=0;i<ns;++i) {
    real sum=0;
    for(p=Ass_rp[i],pe=Ass_rp[i+1];p!=pe;++p) {
      uint j=Ass_j[p];
      if(j>=ns) break;
      sum+=Ass[p]*vs[j];
    }
    Svs[i]=sum;
  }
  for(i=0;i<ln;++i) vl[i]=0;
  apply_p_Als(vl,data,vs,ns);
  sparse_cholesky_solve(vl,&data->fac_A_ll,vl);
  apply_m_Asl(Svs,ns,data,vl);
}

/* vx = X' * vs */
static void apply_Xt(real *vx, uint nx, const xxt_data *data, const real *vs)
{
  const real *X = data->X; const uint *Xp = data->Xp;
  uint i; for(i=0;i<nx;++i) vx[i]=inner(vs,X+Xp[i],Xp[i+1]-Xp[i]);
}

/* vs = X * vx */
static void apply_X(real *vs, uint ns, const xxt_data *data,
                    const real *vx, uint nx)
{
  const real *X = data->X; const uint *Xp = data->Xp;
  uint i,j;
  for(i=0;i<ns;++i) vs[i]=0;
  for(i=0;i<nx;++i) {
    const real v = vx[i];
    const real *x = X+Xp[i]; uint n=Xp[i+1]-Xp[i];
    for(j=0;j<n;++j) vs[j]+=x[j]*v;
  }
}

static void allocate_X(xxt_data *data, sint *perm_x2c)
{
  uint xn=data->xn;
  uint i,h=0;
  if(data->null_space && xn) --xn;
  data->Xp = tmalloc(uint,xn+1);
  data->Xp[0]=0;
  for(i=0;i<xn;++i) {
    if(perm_x2c[i]!=-1) ++h;
    data->Xp[i+1]=data->Xp[i]+h;
  }
  data->X = tmalloc(real,data->Xp[xn]);
}

static void orthogonalize(xxt_data *data, csr_mat *A_ss, sint *perm_x2c,
                          buffer *buf)
{
  uint ln=data->ln, sn=data->sn, xn=data->xn;
  real *vl, *vs, *vx, *Svs;
  uint i,j;

  allocate_X(data,perm_x2c);
  
  buffer_reserve(buf,(ln+2*sn+xn)*sizeof(real));
  vl=buf->ptr, vs=vl+ln, Svs=vs+sn, vx=Svs+sn;

  if(data->null_space && xn) --xn;
  for(i=0;i<xn;++i) {
    uint ns=data->Xp[i+1]-data->Xp[i];
    sint ui = perm_x2c[i];
    real ytsy, *x;
    /*if(data->pid==0) printf("xxt : column %d    %d\n",i,data->pid);*/
    if(ui == -1) {
      for(j=0;j<i;++j) vx[j]=0;
    } else {
      ui-=ln;
      apply_S_col(vs, data,A_ss, ui, vl);
      apply_Xt(vx,i, data, vs);
    }
    apply_QQt(data,vx,i,xn-i);
    apply_X(vs,ns, data, vx,i);
    if(ui!=-1) vs[ui]=1;
    apply_S(Svs,ns, data,A_ss, vs, vl);
    ytsy = inner(vs,Svs,ns);
    ytsy = sum(data,ytsy,i+1,xn-(i+1));
    if(ytsy<EPS) ytsy=0; else ytsy = 1/sqrtr(ytsy);
    x=&data->X[data->Xp[i]];
    for(j=0;j<ns;++j) x[j]=ytsy*vs[j];
  }
}

#endif

/* produces CSR matrix from Yale-like format, summing duplicates */
static void condense_matrix(tuple_list *mat, uint nr, csr_mat *out, buffer *buf)
{
  uint *p, k, nz=mat->n;
  tuple_list_sort(mat,1,buf);
  tuple_list_sort(mat,0,buf);
  
  p = (uint*)mat->vi;
  for(k=0;k+1<nz;++k,p+=2) if(p[0]==p[2] && p[1]==p[3]) break;
  if(++k<nz) {
    uint i=p[0],j=p[1];
    real *pr = &mat->vr[k-1];
    const real *qr = &mat->vr[k];
    const uint *q = p+2;
    for(;k<nz;++k,q+=2,++qr) {
      if(i==q[0]&&j==q[1]) *pr += *qr, --mat->n;
      else p+=2, p[0]=i=q[0],p[1]=j=q[1], *(++pr)=*qr;
    }
  }
  
  nz=mat->n;
  out->n=nr;
  out->Arp = tmalloc(uint,nr+1+mat->n);
  out->Aj = out->Arp+nr+1;
  out->A = tmalloc(real,mat->n);
  for(k=0;k<nr;++k) out->Arp[k]=0;
  for(k=0;k<nz;++k) 
    out->Arp[mat->vi[2*k]]++,
    out->Aj[k]=mat->vi[2*k+1],
    out->A[k]=mat->vr[k];
  nz=0; for(k=0;k<=nr;++k) { uint t=out->Arp[k]; out->Arp[k]=nz, nz+=t; }
}

static void separate_matrix(
  uint nz, const uint *Ai, const uint *Aj, const real *A,
  const sint *perm, uint ln, uint sn,
  csr_mat *out_ll, csr_mat *out_sl, csr_mat *out_ss,
  buffer *buf
)
{
  uint k,n;
  tuple_list mat_ll, mat_sl, mat_ss;
  tuple_list_init_max(&mat_ll,2,0,1,2*nz);
  tuple_list_init_max(&mat_sl,2,0,1,2*nz);
  tuple_list_init_max(&mat_ss,2,0,1,2*nz);
  for(k=0;k<nz;++k) {
    sint i=perm[Ai[k]], j=perm[Aj[k]];
    if(i<0 || j<0 || Aj[k]<Ai[k] || A[k]==0) continue;
    if(j<i) { sint k=j; j=i,i=k; }
    if((uint)i<ln) {
      if((uint)j<ln) {
        n=mat_ll.n++,mat_ll.vi[2*n]=i,mat_ll.vi[2*n+1]=j,mat_ll.vr[n]=A[k];
        if(i!=j)
        n=mat_ll.n++,mat_ll.vi[2*n]=j,mat_ll.vi[2*n+1]=i,mat_ll.vr[n]=A[k];
      } else
        n=mat_sl.n++,mat_sl.vi[2*n]=j-ln,mat_sl.vi[2*n+1]=i,mat_sl.vr[n]=A[k];
    } else {
      n=mat_ss.n++,mat_ss.vi[2*n]=i-ln,mat_ss.vi[2*n+1]=j-ln,mat_ss.vr[n]=A[k];
      if(i!=j)
      n=mat_ss.n++,mat_ss.vi[2*n]=j-ln,mat_ss.vi[2*n+1]=i-ln,mat_ss.vr[n]=A[k];
    }
  }
  condense_matrix(&mat_ll,ln,out_ll,buf);
  condense_matrix(&mat_sl,sn,out_sl,buf);
  condense_matrix(&mat_ss,sn,out_ss,buf);
  tuple_list_free(&mat_ll);
  tuple_list_free(&mat_sl);
  tuple_list_free(&mat_ss);
}

#ifdef MPI

xxt_data *xxt_setup(uint n, const ulong *id,
                    uint nz, const uint *Ai, const uint *Aj, const real *A,
                    uint null_space, crystal_data *crystal)
{
  xxt_data *data = tmalloc(xxt_data,1);
  sint *perm_x2c;
  tuple_list dof;
  csr_mat A_ll, A_ss;

  MPI_Comm_dup(crystal->comm,&data->comm);
  data->pid = crystal->id;
  data->np  = crystal->num;
  locate_proc(data);

  data->null_space=null_space;

  discover_dofs(data,n,id,&dof,crystal);
  discover_sep_sizes(data,&dof,&crystal->all->buf);
#if 0
  if(data->pid>=0)
  printf("xxt %d: un=%d cn=%d ln=%d sn=%d xn=%d proc=%d\n",
         data->pid, data->un, data->cn, data->ln, data->sn, data->xn, data->pid);
#endif
  perm_x2c = discover_sep_ids(data,&dof,&crystal->all->buf);
  if(data->null_space) {
    uint i; real count = 0;
    for(i=0;i<data->cn;++i) count+=1/(real)dof.vi[dof_mi*i+dof_count];
    count=1/sum(data,count,data->xn,0);
    data->share_weight=tmalloc(real,data->cn);
    for(i=0;i<data->cn;++i)
      data->share_weight[i]=count/dof.vi[dof_mi*i+dof_count];
  }
  tuple_list_free(&dof);

  if(!data->null_space || data->xn!=0) {
    separate_matrix(nz,Ai,Aj,A,data->perm_u2c,
                    data->ln,data->sn,
                    &A_ll,&data->A_sl,&A_ss,
                    &crystal->all->buf);
  } else {
    separate_matrix(nz,Ai,Aj,A,data->perm_u2c,
                    data->ln-1,1,
                    &A_ll,&data->A_sl,&A_ss,
                    &crystal->all->buf);
  }                

#if 0
  if(data->pid==1) {
    csr_mat *m[] = {&A_ll,&data->A_sl, &A_ss};
    char *name[] = {"A_ll",     "A_sl","A_ss"};
    uint mi,i,p,pe;
    for(mi=0;mi<3;++mi) {
      printf("%s:\n", name[mi]);
      for(i=0;i<m[mi]->n;++i) {
        for(p=m[mi]->Arp[i],pe=m[mi]->Arp[i+1];p!=pe;++p)
          printf(" (%d,%d) %g\n",i,m[mi]->Aj[p],m[mi]->A[p]);
      }
    }
  }
#endif

  sparse_cholesky_factor(A_ll.n,A_ll.Arp,A_ll.Aj,A_ll.A,
                         &data->fac_A_ll, &crystal->all->buf);
  free(A_ll.Arp); free(A_ll.A);

  data->vl = tmalloc(real,(data->ln+data->cn+2*data->xn)*sizeof(real));
  data->vc = data->vl+data->ln;
  data->vx = data->vc+data->cn;
  data->combuf = data->vx+data->xn;

  orthogonalize(data,&A_ss,perm_x2c,&crystal->all->buf);
  free(A_ss.Arp); free(A_ss.A);
  free(perm_x2c);

  return data;
}

void xxt_solve(real *x, xxt_data *data, const real *b)
{
  uint cn=data->cn, un=data->un, ln=data->ln, sn=data->sn, xn=data->xn;
  real *vl=data->vl, *vc=data->vc, *vx=data->vx;
  uint i;
  for(i=0;i<cn;++i) vc[i]=0;
  for(i=0;i<un;++i) {
    sint p=data->perm_u2c[i];
    if(p>=0) vc[p]+=b[i];
  }
  if(xn>0 && (!data->null_space || xn>1)) {
    if(data->null_space) --xn;
    sparse_cholesky_solve(vc,&data->fac_A_ll,vc);
    apply_m_Asl(vc+ln,sn, data, vc);
    apply_Xt(vx,xn, data, vc+ln);
    apply_QQt(data,vx,xn,0);
    apply_X(vc+ln,sn, data, vx,xn);
    for(i=0;i<ln;++i) vl[i]=0;
    apply_p_Als(vl, data, vc+ln,sn);
    sparse_cholesky_solve(vl,&data->fac_A_ll,vl);
    for(i=0;i<ln;++i) vc[i]-=vl[i];
  } else {
    sparse_cholesky_solve(vc,&data->fac_A_ll,vc);
    if(data->null_space) {
      if(xn==0) vc[ln-1]=0;
      else if(sn==1) vc[ln]=0;
    }
  }
  if(data->null_space) {
    real s=0;
    for(i=0;i<cn;++i) s+=data->share_weight[i]*vc[i];
    s = sum(data,s,data->xn,0);
    for(i=0;i<cn;++i) vc[i]-=s;
  }
  for(i=0;i<un;++i) {
    sint p=data->perm_u2c[i];
    x[i] = p>=0 ? vc[p] : 0;
  }
}

void xxt_stats(xxt_data *data)
{
  int a,b; uint xcol;
  if(data->pid==0) {
    unsigned s;
    printf("xxt: separator sizes on %d =",data->pid);
    for(s=0;s<data->nsep;++s) printf(" %d",(int)data->sep_size[s]);
    printf("\n");
    printf("xxt: shared dofs on %d = %d\n",data->pid,data->sn);
  }
  a=data->ln;
  MPI_Reduce(&a,&b,1,MPI_INT,MPI_MAX,0,data->comm);
  if(data->pid==0) printf("xxt: max non-shared dofs = %d\n",b);
  a=data->sn;
  MPI_Reduce(&a,&b,1,MPI_INT,MPI_MAX,0,data->comm);
  if(data->pid==0) printf("xxt: max shared dofs = %d\n",b);
  xcol=data->xn; if(xcol&&data->null_space) --xcol;
  a=xcol;
  MPI_Reduce(&a,&b,1,MPI_INT,MPI_MAX,0,data->comm);
  if(data->pid==0) printf("xxt: max X cols = %d\n",b);
  a=data->Xp[xcol]*sizeof(real);
  MPI_Reduce(&a,&b,1,MPI_INT,MPI_MAX,0,data->comm);
  if(data->pid==0) printf("xxt: max X size = %d bytes\n",b);
}

void xxt_free(xxt_data *data)
{
  MPI_Comm_free(&data->comm);
  free(data->pother);
  free(data->mpireq);
  free(data->mpistatus);
  free(data->sep_size);
  free(data->perm_u2c);
  if(data->null_space) free(data->share_weight);
  sparse_cholesky_free(&data->fac_A_ll);
  free(data->A_sl.Arp); free(data->A_sl.A);
  free(data->Xp); free(data->X);
  free(data->vl);
  free(data);
}

#else

xxt_data *xxt_setup(uint n, const ulong *id,
                    uint nz, const uint *Ai, const uint *Aj, const real *A,
                    uint null_space, void *crystal)
{
  xxt_data *data = tmalloc(xxt_data,1);
  buffer buf;
  csr_mat A_ll, A_sl, A_ss;

  buffer_init(&buf,1024);
  data->null_space=null_space;
  data->un = n;
  data->perm_u2c = tmalloc(sint,n);
  data->cn = unique_ids(n,id,data->perm_u2c,&buf);
  
  separate_matrix(nz,Ai,Aj,A,data->perm_u2c,
                  data->cn-(null_space?1:0),(null_space?1:0),
                  &A_ll,&A_sl,&A_ss,
                  &buf);
#if 0
  if(data->pid==0) {
    csr_mat *m[] = {&A_ll, &A_sl, &A_ss};
    char *name[] = {"A_ll","A_sl","A_ss"};
    uint mi,i,p,pe;
    for(mi=0;mi<3;++mi) {
      printf("%s:\n", name[mi]);
      for(i=0;i<m[mi]->n;++i) {
        for(p=m[mi]->Arp[i],pe=m[mi]->Arp[i+1];p!=pe;++p)
          printf(" (%d,%d) %g\n",i,m[mi]->Aj[p],m[mi]->A[p]);
      }
    }
  }
#endif

  sparse_cholesky_factor(A_ll.n,A_ll.Arp,A_ll.Aj,A_ll.A,
                         &data->fac_A_ll, &buf);
  free(A_ll.Arp); free(A_ll.A);
  free(A_sl.Arp); free(A_sl.A);
  free(A_ss.Arp); free(A_ss.A);

  data->vc = tmalloc(real,data->cn*sizeof(real));

  buffer_free(&buf);
  
  return data;
}

void xxt_solve(real *x, xxt_data *data, const real *b)
{
  uint cn=data->cn, un=data->un;
  real *vc=data->vc;
  uint i;
  for(i=0;i<cn;++i) vc[i]=0;
  for(i=0;i<un;++i) {
    sint p=data->perm_u2c[i];
    if(p>=0) vc[p]+=b[i];
  }
  sparse_cholesky_solve(vc,&data->fac_A_ll,vc);
  if(data->null_space) {
    real s=0;
    vc[cn-1]=0;
    for(i=0;i<cn;++i) s+=vc[i];
    s /= cn;
    for(i=0;i<cn;++i) vc[i]-=s;
  }
  for(i=0;i<un;++i) {
    sint p=data->perm_u2c[i];
    x[i] = p>=0 ? vc[p] : 0;
  }
}

void xxt_stats(xxt_data *data)
{
  printf("xxt: separator sizes on 0 = %d\n",data->cn);
}

void xxt_free(xxt_data *data)
{
  free(data->perm_u2c);
  sparse_cholesky_free(&data->fac_A_ll);
  free(data->vc);
  free(data);
}

#endif

/*--------------------------------------------------------------------------
   FORTRAN Interface
  --------------------------------------------------------------------------*/

#include "fname.h"

#define xxtsetup   FORTRAN_NAME(xxtsetup,XXTSETUP)
#define xxtsolve   FORTRAN_NAME(xxtsolve,XXTSOLVE)
#define xxtstats   FORTRAN_NAME(xxtstats,XXTSTATS)
#define xxtfree    FORTRAN_NAME(xxtfree ,XXTFREE)

static xxt_data **handle_array = 0;
static int handle_max = 0;
static int handle_n = 0;

#ifdef MPI
crystal_data *fcrystal_handle(sint h);
#endif

void xxtsetup(sint *handle, const sint *crystal_handle,
              const sint *n, const slong id[],
              const sint *nz, const sint Ai[], const sint Aj[], const real A[],
              const sint *null_space)
{
#ifdef MPI
  crystal_data *crystal = fcrystal_handle(*crystal_handle);
#else
  void *crystal = 0;
#endif
  if(handle_n==handle_max)
    handle_max+=handle_max/2+1,
    handle_array=trealloc(xxt_data*,handle_array,handle_max);
  handle_array[handle_n]=xxt_setup(*n,(const ulong*)id,
                                   *nz,(const uint*)Ai,(const uint*)Aj,A,
                                   *null_space,crystal);
  *handle = handle_n++;
}

void xxtsolve(const sint *handle, real x[], const real b[])
{
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle])
    failwith("invalid handle to xxtsolve");
  xxt_solve(x,handle_array[*handle],b);
}

void xxtstats(const sint *handle)
{
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle])
    failwith("invalid handle to xxtstats");
  xxt_stats(handle_array[*handle]);
}

void xxtfree(sint *handle)
{
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle])
    failwith("invalid handle to xxtfree");
  xxt_free(handle_array[*handle]);
  handle_array[*handle] = 0;
}


