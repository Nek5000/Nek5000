/*--------------------------------------------------------------------------
   Parallel Point Finding
          
   Initializing the data:
     unsigned ndim = 3;
     unsigned nel;        // number of elements
     const unsigned n[3]; // number of nodes in r, s, t directions
     const real *xm[3];   // n[0]*n[1]*n[2]*nel x,y,z coordinates
     real tol = 0.01;     // how far point is allowed to be outside element
                          //   relative to element size
     unsigned max_size = n[0]*n[1]*n[2]*nel; // maximum size of hash table
     crystal_data *crystal = ... ; // communication handle (see crystal.*)
                          
     pfindpt_data *data = pfindpt_setup(3,xm,n,nel,max_size,tol,crystal);
     
   Using the data:
     tuple_list list;
     
     for point j,
       list.vi[list.mi*j + 0] = processer number
       list.vi[list.mi*j + 1] = element number
       list.vi[list.mi*j + 2] = code
       list.vi[list.mi*j + 3] = open for the user
       ...
       list.vi[list.mi*j + list.mi-1] = open for the user
       
       list.vr[list.mr*j + 0] = distance
       list.vr[list.mr*j + 1] = x
       list.vr[list.mr*j + 2] = y
       list.vr[list.mr*j + 3] = z when ndim==3
       list.vr[list.mr*j + ndim+1] = r
       list.vr[list.mr*j + ndim+2] = s
       list.vr[list.mr*j + ndim+3] = t when ndim==3
       list.vr[list.mr*j + 2*ndim+1] = open for the user
       ...
       list.vr[list.mr*j + list.mr-1] = open for the user

    pfindpt(data, &list, guess);
  
    - On input, only the xyz fields are used; the rest are set as output.
      The exception is if guess is non-zero, then element number and parametric
      coords will be used as an initial guess.
    - The code is set as follows:
         0 : normal
        -1 : point not within mesh (to within given tolerance)
         1 : point either exactly on element boundary, or outside mesh
             (but within given tolerance); in this case the returned distance
             can be used to test if the point is really outside the mesh

  To transfer points:
  
    int dynamic = 1;
    pfindpt_transfer(data, &list, dynamic);
    
    - This is just a small wrapper over transfer()
   
   When done:
     pfindpt_free(&data);
     
  --------------------------------------------------------------------------*/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#ifdef MPI
#  include <string.h>
#  include <mpi.h>
#endif

#include "errmem.h"
#include "types.h"
#include "minmax.h"
#include "poly.h"
#include "tensor.h"
#include "findpt.h"
#include "tuple_list.h"

#ifdef MPI
#  include "crystal.h"
#  include "transfer.h"
#else
typedef void crystal_data;
#endif


#ifdef MPI

static int iceil(real x) { return ceilr(x); }
static int ifloor(real x) { return floorr(x); }

typedef struct  {
  real bnd[6], fac[3];
  unsigned n, ncell;
  uint *cell_offset, *cell;
} hash_data;

static void hash_free(hash_data *h)
{
  free(h->cell_offset), free(h->cell);
}

static void hash_range(const hash_data *h, const real *bnd, unsigned d,
                       unsigned *ia, unsigned *ib)
{
  const real a = bnd[d*2  ], b = bnd[d*2+1];
  const int      i0 = ifloor( (a - h->bnd[d*2]) * h->fac[d] );
  const unsigned i1 = iceil(  (b - h->bnd[d*2]) * h->fac[d] );
  *ia = imax_2(0,i0);
  *ib = umin_2(i1,h->n);
  if(*ib == *ia) ++(*ib);
}

static void hash_size_calc(hash_data *h, unsigned ndim, const real *bnd,
                           const crystal_data *crystal)
{
  unsigned d; uint nn; double out[3], in[3];
   
  for(d=0;d<ndim;++d) out[d]=bnd[d*2];
  MPI_Allreduce(out,in,3,MPI_DOUBLE,MPI_MIN,crystal->comm);
  for(d=0;d<ndim;++d) h->bnd[d*2]=in[d];
  
  for(d=0;d<ndim;++d) out[d]=bnd[d*2+1];
  MPI_Allreduce(out,in,3,MPI_DOUBLE,MPI_MAX,crystal->comm);
  for(d=0;d<ndim;++d) h->bnd[d*2+1]=in[d];

  h->n = iceil(pow(crystal->num,1./ndim));
  for(d=0;d<ndim;++d) h->fac[d] = h->n/(h->bnd[d*2+1] - h->bnd[d*2]);
/* 
  h->fac[1] = h->n/(h->bnd[3] - h->bnd[2]);
  h->fac[2] = h->n/(h->bnd[5] - h->bnd[4]);
*/  
  nn=(uint)h->n*h->n; if(ndim==3) nn*=h->n;
  h->ncell = (nn-1)/crystal->num+1;
  if(crystal->id+(h->ncell-1)*crystal->num>=nn) --h->ncell;
  h->cell_offset = tmalloc(uint,h->ncell+1);
}

static void hash_build(hash_data *h, unsigned ndim, const real *bnd,
                       crystal_data *crystal)
{
  tuple_list tl; sint *work;
  unsigned range[6]={0,1,0,1,0,1};
  uint i, j, k; uint count=1;

  hash_size_calc(h,ndim,bnd,crystal);

  for(i=0;i<ndim;++i)
    hash_range(h,bnd,i,&range[i*2],&range[i*2+1]),
    count*=range[i*2+1]-range[i*2];

  tuple_list_init_max(&tl,2,0,0,count);
  tl.n=count; work=tl.vi;
  for(k=range[4];k!=range[5];++k) {
    for(j=range[2];j!=range[3];++j) {
      for(i=range[0];i!=range[1];++i) {
        uint index = ((uint)k*h->n+j)*h->n+i;
        *work++ = index;
        *work++ = index%crystal->num;
      }
    }
  }

  transfer(1,&tl,1,crystal);

  /* sort on hash index (equiv. to cell) */
  tuple_list_sort(&tl,0,&crystal->all->buf);

  h->cell=tmalloc(uint,tl.n);
  j=0; h->cell_offset[0]=0; h->cell_offset[j+1]=0;
  work = tl.vi;
  for(i=tl.n;i;--i,work+=2) {
    const uint cell = work[0]/crystal->num,
               src  = work[1];
    for(k=j+1;k<=cell;++k) h->cell_offset[k+1]=h->cell_offset[k];
    j=cell;
    h->cell[h->cell_offset[j+1]++] = src;
  }  
  for(k=j+1;k<h->ncell;++k) h->cell_offset[k+1]=h->cell_offset[k];
  tuple_list_free(&tl);
}

#endif

typedef struct {
  unsigned ndim;
  void *fd;
  findpt_func findpt;
#ifdef MPI
  hash_data hash;
  crystal_data *crystal;
  tuple_list tl_hash, /* vi:{id,p,cell}    vd:{x,y,z}      */
             tl_in,   /* vi:{id,p,sp}      vd:{x,y,z}      */
             tl_out;  /* vi:{id,p,el,code} vd:{r,s,t,dist} */
#endif
} pfindpt_data;

#define IW_ID   0
#define IW_PROC 1
#define IW_CELL 2
#define IW_SRCP 2
#define IW_EL   2
#define IW_CODE 3

#define DW_X 0
#define DW_R 0
#define DW_DIST ndim

pfindpt_data *pfindpt_setup(unsigned ndim, const real *const*xw,
                            const unsigned *n, uint nel,
                            uint max_hash_size, real bbox_tol,
                            crystal_data *crystal)
{
  pfindpt_data *p = tmalloc(pfindpt_data,1);
  const real *bnd;
  p->ndim = ndim;
#ifdef MPI
  p->crystal = crystal;
#endif
  if(ndim==2) {
    p->fd = findpt_setup_2(xw,n,nel,max_hash_size,bbox_tol);
    p->findpt = (findpt_func)&findpt_2;
    bnd = findpt_allbnd_2((findpt_data_2*)p->fd);
  } else if(ndim==3) {
    p->fd = findpt_setup_3(xw,n,nel,max_hash_size,bbox_tol);
    p->findpt = (findpt_func)&findpt_3;
    bnd = findpt_allbnd_3((findpt_data_3*)p->fd);
  } else
    bnd=0,fail("%s: pfindpt_setup: parameter ndim=%u invalid",__FILE__,ndim);
#ifdef MPI
  hash_build(&p->hash,ndim,bnd,p->crystal);
  tuple_list_init(&p->tl_hash,3,0,ndim);
  tuple_list_init(&p->tl_in,3,0,ndim);
  tuple_list_init(&p->tl_out,4,0,ndim+1);
#endif
  return p;
}

void pfindpt_free(pfindpt_data *p)
{
  if(p->ndim==2) findpt_free_2(p->fd); else findpt_free_3(p->fd);
#ifdef MPI
  hash_free(&p->hash);
  tuple_list_free(&p->tl_hash);
  tuple_list_free(&p->tl_in);
  tuple_list_free(&p->tl_out);
#endif
  free(p);
}

#define I_PROC 0
#define I_EL   1
#define I_CODE 2

#define D_DIST 0
#define D_X    1
#define D_R    (D_X+ndim)

#define GROW(a,b) (a+a/2+1>b?a+a/2+1:b)

void pfindpt_transfer(pfindpt_data *p, tuple_list *list, int dynamic)
{
#ifdef MPI
  transfer(dynamic,list,I_PROC,p->crystal);
#endif  
}

void pfindpt(pfindpt_data *p, tuple_list *list, int guess)
{
  const int lmi=list->mi, lmr=list->mr;
  const int ndim=p->ndim; int d;
  uint i, id;
  sint *ri; real *rr;
#ifdef MPI
  sint *oi; real *or;
  const uint np=p->crystal->num;
  p->tl_hash.n = p->tl_in.n = p->tl_out.n = 0;
  if(p->tl_hash.max<list->n)
    tuple_list_resize(&p->tl_hash,GROW(p->tl_hash.max,list->n));
  oi = p->tl_hash.vi, or = p->tl_hash.vr;
  id = p->crystal->id;
#else
  id = 0;
#endif
  ri=list->vi,rr=list->vr;
  for(i=0;i<list->n;++i,ri+=lmi,rr+=lmr) {
    ri[I_PROC]=id;
    ri[I_CODE]=p->findpt(p->fd,&rr[D_X],guess,
                         (uint*)&ri[I_EL],&rr[D_R],&rr[D_DIST]);
#ifdef MPI
    if(ri[I_CODE]!=0) { /* point not found, or point on boundary */
      int hi; uint hash=0;
      for(d=ndim-1;d>=0;--d) {
        hi = ifloor((rr[D_X+d]-p->hash.bnd[d*2])*p->hash.fac[d]);
        if(hi<0 || (unsigned)hi>=p->hash.n) { ri[I_CODE]=-1; break; }
        hash = hash*p->hash.n+hi;
      }
      if(d>=0) continue;
      oi[IW_ID]=i;
      oi[IW_PROC]=hash%np;
      oi[IW_CELL]=hash/np;
      memcpy(or,&rr[D_X],ndim*sizeof(real));
      oi+=3, or+=ndim, ++p->tl_hash.n;
    }
#endif
  }

#ifdef MPI  

  transfer(1,&p->tl_hash,IW_PROC,p->crystal);
  
  ri=p->tl_hash.vi,rr=p->tl_hash.vr;
  oi=p->tl_in.vi,or=p->tl_in.vr;
  for(i=0;i<p->tl_hash.n;++i,ri+=3,rr+=ndim) {
    uint *proc = p->hash.cell + p->hash.cell_offset[ri[IW_CELL]],
         *pend = p->hash.cell + p->hash.cell_offset[ri[IW_CELL]+1];
    for(;proc!=pend;++proc) {
      if((sint)*proc==ri[IW_PROC]) continue;
      if(p->tl_in.n==p->tl_in.max) {
        tuple_list_grow(&p->tl_in);
        oi=p->tl_in.vi+3*p->tl_in.n;
        or=p->tl_in.vr+ndim*p->tl_in.n;
      }
      oi[IW_ID]=ri[IW_ID];
      oi[IW_SRCP]=ri[IW_PROC];
      oi[IW_PROC]=*proc;
      memcpy(or,rr,ndim*sizeof(real));
      oi+=3, or+=ndim, ++p->tl_in.n;
    }
  }
  
  transfer(1,&p->tl_in,IW_PROC,p->crystal);

  if(p->tl_out.max<p->tl_in.n)
    tuple_list_resize(&p->tl_out,GROW(p->tl_out.max,p->tl_in.n));
  ri=p->tl_in.vi,rr=p->tl_in.vr;
  oi=p->tl_out.vi,or=p->tl_out.vr;
  for(i=0;i<p->tl_in.n;++i,ri+=3,rr+=ndim) {
    oi[IW_CODE]=p->findpt(p->fd,&rr[DW_X],0,
                          (uint*)&oi[IW_EL],&or[DW_R],&or[DW_DIST]);
    if(oi[IW_CODE]==-1) continue;
    oi[IW_ID]=ri[IW_ID];
    oi[IW_PROC]=ri[IW_SRCP];
    oi+=4, or+=ndim+1, ++p->tl_out.n;
  }

  transfer(1,&p->tl_out,IW_PROC,p->crystal);

  ri=p->tl_out.vi,rr=p->tl_out.vr;
  for(i=0;i<p->tl_out.n;++i,ri+=4,rr+=ndim+1) {
    oi = list->vi+lmi*ri[IW_ID];
    or = list->vr+lmr*ri[IW_ID];
    if(oi[I_CODE]!=0) {
      if(ri[IW_CODE]==1 && oi[I_CODE]==1 && rr[DW_DIST]>=or[D_DIST]) continue;
      oi[I_PROC] = ri[IW_PROC];
      oi[I_CODE] = ri[IW_CODE];
      or[D_DIST] = rr[DW_DIST];
      oi[I_EL]   = ri[IW_EL];
      memcpy(&or[D_R],&rr[DW_R],ndim*sizeof(real));
    }
  }
  
#endif
}

void pfindpt_weights(pfindpt_data *p, const real *r)
{
  if(p->ndim==2) findpt_weights_2(p->fd,r);
            else findpt_weights_3(p->fd,r);
}

real pfindpt_eval(pfindpt_data *p, const real *u)
{
  if(p->ndim==2) return findpt_eval_2(p->fd,u);
            else return findpt_eval_3(p->fd,u);
}

