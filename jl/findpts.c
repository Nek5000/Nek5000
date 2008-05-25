/*------------------------------------------------------------------------------
  
  FORTRAN interface for pfindpt
  
    integer ch
    call crystal_new(ch)
  
    integer h
    integer ndim ! = 2 or 3
    integer nr, ns, nt, nel
    real xm1(nr,ns,nt,nel), ym1(nr,ns,nt,nel), zm1(nr,ns,nt,nel)
    real tolerance

    call findpts_new(h,ch,ndim,xm1,ym1,zm1,nr,ns,nt,nel,tolerance)
  
  The returned handle will use the given crystal router handle for
  communication, and will be able to locate lists of points that are
  within the given mesh to the specified tolerance relative to element size
  (e.g. tolerance = 0.1 or 1.0e-10)
  
  The list of points must have the format:
  
    integer mi ! >= 3
    integer mr ! >= 1 + 2*ndim
    integer n, max
    integer vi(mi,max)
    real    vr(mr,max)
  
  For point j, j in [1 ... n], n <= max,
    vi(1,j) = processor number (0 to P-1)
    vi(2,j) = element number (1 to vi(1,j)'s nel)
    vi(3,j) = code (-1, 0, 1; explained below)
    vr(1,j) = distance (from located point to given point)
    vr(2,j) = x  (input)
    vr(3,j) = y  (input)
    vr(4,j) = z  (input; only when ndim=3)
    vr(ndim+2,j) = r   (output)
    vr(ndim+3,j) = s   (output)
    vr(ndim+4,j) = t   (output; only when ndim=3)

  To locate points:
  
    call findpts(h, n, vi, mi, vr, mr, guess)
  
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
  
    call findpts_transfer(h,n,max,vi,mi,vr,mr)
    
    - This is just a small wrapper over crystal_transfer
    - The error condition (more incoming points than max) is indicated by
        n .eq. max + 1    on return;
      so, to ignore errors one must use:
    
      call findpts_transfer(h,n,max,vi,mi,vr,mr)
      if(n.eq. max+1) n=max
  
  To evaluate scalar fields for point j (must be local):
  
    real scalar1(nr,ns,nt,nel), scalar2(nr,ns,nt,nel)
    real val1, val2
    real r, s, t
    integer el

    r = vr(ndim+2,j)
    s = vr(ndim+3,j)
    if(ndim.eq.3) t = vr(ndim+4,j)
    el = vi(1,j)
    
    call findpts_weights(h, r,s,t)
    call findpts_eval(h, val1, scalar1(1,1,1,el)) ! sets val1
    call findpts_eval(h, val2, scalar2(1,1,1,el)) ! sets val2
    ...
    
  ----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#ifdef MPI
#  include <mpi.h>
#endif

#include "fname.h"
#include "errmem.h"
#include "types.h"
#ifdef MPI
#  include "crystal.h"
#endif
#include "tuple_list.h"
#include "pfindpt.h"

#define findpts_new      FORTRAN_NAME(findpts_new,FINDPTS_NEW)
#define findpts_done     FORTRAN_NAME(findpts_done,FINDPTS_DONE)
#define findpts_transfer FORTRAN_NAME(findpts_transfer,FINDPTS_TRANSFER)
#define findpts          FORTRAN_NAME(findpts,FINDPTS)
#define findpts_weights  FORTRAN_NAME(findpts_weights,FINDPTS_WEIGHTS)
#define findpts_eval     FORTRAN_NAME(findpts_eval,FINDPTS_EVAL)

static pfindpt_data **handle=0;
static unsigned n=0, max=0;

#ifdef MPI
crystal_data *fcrystal_handle(sint h);
#endif

void findpts_new(sint *h, const sint *ch, const sint *ndim,
                 const real xm1[], const real ym1[], const real zm1[],
                 const sint *nr, const sint *ns, const sint *nt,
                 const sint *nel, const real *tol)
{
  const real *xw[3] = {xm1,ym1,zm1};
  unsigned nm[3] = {*nr,*ns,*nt};
#ifdef MPI
  crystal_data *crystal = fcrystal_handle(*ch);
#else
  void *crystal = 0;
#endif
  if(n==max) max+=max/2+1,handle=trealloc(pfindpt_data*,handle,max);
  if(*ndim==2) nm[2]=1;
  handle[n] = pfindpt_setup(*ndim, xw, nm, *nel, *nel*nm[0]*nm[1]*nm[2], *tol,
                            crystal);
  *h=n++;
}

static pfindpt_data *findpts_handle(sint h)
{
  if((unsigned)h>=n || handle[h]==0) failwith("invalid findpts handle");
  return handle[h];
}

void findpts_done(sint *h)
{
  pfindpt_data *p = findpts_handle(*h);
  handle[*h]=0;
  pfindpt_free(p);
}

void findpts_transfer(const sint *h, sint *n, const sint *max, sint vi[],
                      const sint *mi, real vr[], const sint *mr)
{
  tuple_list tl = {*mi,0,*mr, *n,*max, vi,0,vr};
  pfindpt_transfer(findpts_handle(*h),&tl,0);
  *n = tl.n;
}

#define I_EL   1
void findpts(const sint *h, const sint *n, sint vi[], const sint *mi,
             real vr[], const sint *mr, const sint *guess)
{
  uint i; sint *ri;
  tuple_list tl = {*mi,0,*mr, *n,*n, vi,0,vr};
  if(*guess) {
    ri = &vi[I_EL];
    for(i=0;i<tl.n;++i,ri+=tl.mi) --*ri;
  }
  pfindpt(findpts_handle(*h),&tl,*guess);
  ri = &vi[I_EL];
  for(i=0;i<tl.n;++i,ri+=tl.mi) ++*ri;
}

void findpts_weights(const sint *h, const real *r, const real *s, const real *t)
{
  const real ar[3] = {*r,*s,*t};
  pfindpt_weights(findpts_handle(*h),ar);
}

void findpts_eval(const sint *h, real *out, const real u[])
{
  *out = pfindpt_eval(findpts_handle(*h), u);
}

