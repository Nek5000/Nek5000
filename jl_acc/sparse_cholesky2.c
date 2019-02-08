#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "types.h"
#include "minmax.h"
#include "errmem.h"
#include "sort.h"

/* factors: L is in CSR format
            D is a diagonal matrix stored as a vector
   actual factorization is:
   
                  -1      T
     A   = (I-L) D   (I-L)
     
      -1        -T        -1
     A   = (I-L)   D (I-L)
     
   (triangular factor is unit diagonal; the diagonal is not stored)
*/
typedef struct {
  uint n, *Lrp, *Lj;
  real *L, *D;
} sparse_cholesky_data;

/*
  symbolic factorization: finds the sparsity structure of L
  
  uses the concept of elimination tree:
    the parent of node j is node i when L(i,j) is the first
      non-zero in column j below the diagonal (i>j)
    L's structure is discovered row-by-row; the first time
      an entry in column j is set, it must be the parent

  the nonzeros in L are the nonzeros in A + paths up the elimination tree  
  
  linear in the number of nonzeros of L
*/
static void factor_symbolic(uint n, const uint *Arp, const uint *Aj,
                            sparse_cholesky_data *out, uint *work)
{
  uint *visit, *parent, *sorted;
  uint *Lrp, *Lj;
  uint i,nz=0;

  out->n=n;
  
  /* sorted needs 2*n; work needs 4*n */
  visit = work, parent = visit+n, sorted=parent+n;

  for(i=0;i<n;++i) {
    uint p=Arp[i], pe=Arp[i+1];
    visit[i]=i, parent[i]=n;
    for(;p!=pe;++p) {
      uint j=Aj[p]; if(j>=i) break;
      for(;visit[j]!=i;j=parent[j]) {
        ++nz, visit[j]=i;
        if(parent[j]==n) { parent[j]=i; break; }
      }
    }
  }
  
  Lrp=out->Lrp=tmalloc(uint,n+1+nz);
  Lj =out->Lj =Lrp+n+1;

  Lrp[0]=0;
  for(i=0;i<n;++i) {
    uint p=Arp[i], pe=Arp[i+1], count=0, *Ljr=&Lj[Lrp[i]];
    visit[i]=i;
    for(;p!=pe;++p) {
      uint j=Aj[p]; if(j>=i) break;
      for(;visit[j]!=i;j=parent[j]) Ljr[count++]=j, visit[j]=i;
    }
    sort(Ljr,count,1,sorted,sorted+count);
    memcpy(Ljr,sorted,count*sizeof(uint));
    Lrp[i+1]=Lrp[i]+count;
  }
}

/*
  numeric factorization:
  
  L is built row-by-row, using:    ( ' indicates transpose )

  
  [ A  r ]  = [ (I-L)   ] [ D^(-1)  ] [ (I-L)' -s ]
  [ r' a ]    [  -s'  1 ] [     1/d ] [         1 ]
            
            = [ A   (I-L) D^(-1) (-s)  ]
              [ r'  s' D^(-1) s + 1/d  ]
              
  so, if r' is the next row of A, up to but excluding the diagonal,
  then the next row of L, s', obeys
  
     r = - (I-L) D^(-1) s
 
  let y = (I-L)^(-1) (-r)
  then s = D y, and d = 1/(s' y)
  
*/
static void factor_numeric(uint n, const uint *Arp, const uint *Aj,
                           const real *A,
                           sparse_cholesky_data *out,
                           uint *visit, real *y)
{
  const uint *Lrp=out->Lrp, *Lj=out->Lj;
  real *D, *L;
  uint i;
  
  D=out->D=tmalloc(real,n+Lrp[n]);
  L=out->L=D+n;
  
  for(i=0;i<n;++i) {
    uint p,pe; real a=0;
    visit[i]=n;
    for(p=Lrp[i],pe=Lrp[i+1];p!=pe;++p) {
      uint j=Lj[p]; y[j]=0, visit[j]=i;
    }
    for(p=Arp[i],pe=Arp[i+1];p!=pe;++p) {
      uint j=Aj[p];
      if(j>=i) { if(j==i) a=A[p]; break; }
      y[j]=-A[p];
    }
    for(p=Lrp[i],pe=Lrp[i+1];p!=pe;++p) {
      uint q,qe,j=Lj[p]; real lij,yj=y[j];
      for(q=Lrp[j],qe=Lrp[j+1];q!=qe;++q) {
        uint k=Lj[q]; if(visit[k]==i) yj+=L[q]*y[k];
      }
      y[j]=yj;
      L[p]=lij=D[j]*yj;
      a-=yj*lij;
    }
    D[i]=1/a;
  }
}

/* x = A^(-1) b;  works when x and b alias */
void sparse_cholesky_solve(real *x, const sparse_cholesky_data *fac, real *b)
{
  const uint n=fac->n, *Lrp=fac->Lrp, *Lj=fac->Lj;
  const real *L=fac->L, *D=fac->D;
  uint i, p,pe;
  for(i=0;i<n;++i) {
    real xi = b[i];
    for(p=Lrp[i],pe=Lrp[i+1];p!=pe;++p) xi+=L[p]*x[Lj[p]];
    x[i]=xi;
  }
  for(i=0;i<n;++i) x[i]*=D[i];
  for(i=n;i;) {
    real xi = x[--i];
    for(p=Lrp[i],pe=Lrp[i+1];p!=pe;++p) x[Lj[p]]+=L[p]*xi;
  }
}

void sparse_cholesky_factor(uint n, const uint *Arp, const uint *Aj,
                            const real *A,
                            sparse_cholesky_data *out, buffer *buf)
{
  const uint n_uints_as_reals = (n*sizeof(uint)+sizeof(real)-1)/sizeof(real);
  buffer_reserve(buf,umax_2(4*n*sizeof(uint),
                            (n_uints_as_reals+n)*sizeof(real)));
  factor_symbolic(n,Arp,Aj,out,buf->ptr);
  factor_numeric(n,Arp,Aj,A,out,buf->ptr,n_uints_as_reals+(real*)buf->ptr);
}

void sparse_cholesky_free(sparse_cholesky_data *fac)
{
  free(fac->Lrp); fac->Lj=fac->Lrp=0;
  free(fac->D);   fac->L =fac->D  =0;
}

