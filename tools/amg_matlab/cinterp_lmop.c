#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "cinterp_common.h"

/*--------------------------------------------------------------------------
   sparse add
   
   y += alpha * x
   
   the sparse vector x is added to y
   it assumed that yi and xi are sorted,
   and that xi is a subset of yi
--------------------------------------------------------------------------*/
static void sp_add(mwSize yn, const mwIndex *yi, double *y, double alpha,
                   mwSize xn, const mwIndex *xi, const double *x)
{
  const mwIndex *xe = xi+xn;
  mwIndex iy;
  if(yn==0) return; iy = *yi;
  for(;xi!=xe;++xi,++x) {
    mwIndex ix = *xi;
    while(iy<ix) ++y, iy=*(++yi);
    *y++ += alpha * (*x), iy=*(++yi);
  }
}

static int interp_lmop(const sp_mat *T, const sp_mat_c *A, const double *u,
                        const sp_log_mat_c *X_skel)
{
  mwSize nf = X_skel->m, nc = X_skel->n;
  mwSize max_nz=0, max_Q;
  mwIndex j;
  
  double *sqv1, *sqv2;
  double *Q, *QQt;
  
  for(j=0;j<nc;++j) {
    mwSize nz=X_skel->jc[j+1]-X_skel->jc[j];
    if(nz>max_nz) max_nz=nz;
  }
  max_Q = (max_nz*(max_nz+1))/2;
  
  if(!(sqv1=mem_alloc((2*max_nz + max_Q + max_nz*max_nz)*sizeof(double))))
    return 0;
  sqv2 = sqv1+max_nz, Q = sqv2+max_nz, QQt = Q+max_Q;

  { mwIndex nz=T->jc[nf]; for(j=0;j<nz;++j) T->pr[j]=0; }
  for(j=0;j<nc;++j) {
    const mwIndex *Qi = &X_skel->ir[X_skel->jc[j]];
    mwSize nz = X_skel->jc[j+1]-X_skel->jc[j];
    mwIndex m,k;
    double *qk = Q;
    double uj2 = u[j]*u[j];
    for(k=0;k<nz*nz;++k) QQt[k]=0;
    for(k=0;k<nz;++k,qk+=k) {
      double alpha;
      mwIndex s = Qi[k];
      /* sqv1 := R_(k+1) A e_s */
      sp_restrict_sorted(sqv1, k+1,Qi, A->jc[s+1]-A->jc[s],
        &A->ir[A->jc[s]], &A->pr[A->jc[s]]);
      /* sqv2 := Q^t A e_s */
      mv_utt(sqv2, k,Q, sqv1);
      /* qk := Q Q^t A e_s */ 
      mv_ut(qk, k,Q, sqv2);
      /* alpha := ||(I-Q Q^t A)e_s||_A^2 = (A e_s)^t (I-Q Q^t A)e_s */
      alpha = sqv1[k];
      for(m=0;m<k;++m) alpha -= sqv1[m] * qk[m];
      /* qk := Q e_(k+1) = alpha^{-1/2} (I-Q Q^t A)e_s */
      alpha = -1.0 / sqrt(alpha);
      for(m=0;m<k;++m) qk[m] *= alpha;
      qk[k] = -alpha;
      /* QQt := QQt + qk qk^t */
      for(m=0;m<=k;++m) {
        mwIndex i, mnz = m*nz; double qkm = qk[m];
        for(i=0;i<=k;++i) QQt[mnz+i] += qkm * qk[i];
      }
    }
    /* T := T + u_j^2 QQt */
    qk=QQt;
    for(k=0;k<nz;++k,qk+=nz) {
      mwIndex i = Qi[k], ti = T->jc[i];
      sp_add(T->jc[i+1]-ti,&T->ir[ti],&T->pr[ti], uj2, nz,Qi,qk);
    }
  }
  mem_free(sqv1);
  return 1;
}

/* A, u, X_skel, T_skel */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  sp_mat T;
  sp_mat_c A;
  const double *u;
  sp_log_mat_c X_skel, T_skel;
  if(nrhs!=4) { mexWarnMsgTxt("Four inputs required."); return; }
  if(!mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
    mexWarnMsgTxt("First input not a sparse, double, real matrix."); return;
  }
  if(mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
    mexWarnMsgTxt("Second input not a full, double, real array."); return;
  }
  if(!mxIsSparse(prhs[2])) {
    mexWarnMsgTxt("Third input not a sparse matrix."); return;
  }
  if(!mxIsSparse(prhs[3])) {
    mexWarnMsgTxt("Four input not a sparse matrix."); return;
  }
  A.m = mxGetM(prhs[0]), A.n = mxGetN(prhs[0]);
  A.ir = mxGetIr(prhs[0]), A.jc = mxGetJc(prhs[0]), A.pr = mxGetPr(prhs[0]);
  u = mxGetPr(prhs[1]);
  X_skel.m = mxGetM(prhs[2]), X_skel.n = mxGetN(prhs[2]);
  X_skel.ir = mxGetIr(prhs[2]), X_skel.jc = mxGetJc(prhs[2]);
  T_skel.m = mxGetM(prhs[3]), T_skel.n = mxGetN(prhs[3]);
  T_skel.ir = mxGetIr(prhs[3]), T_skel.jc = mxGetJc(prhs[3]);
  if(A.m!=A.n) { mexWarnMsgTxt("A not square."); return; }
  if(A.m!=X_skel.m) { mexWarnMsgTxt("rows(A) != rows(X_skel)"); return; }
  if(A.m!=T_skel.m || A.n!=T_skel.n) {
    mexWarnMsgTxt("A, T_skel must be the same shape."); return;
  }
  if(mxGetM(prhs[1])!=X_skel.n || mxGetN(prhs[1])!=1) {
    mexWarnMsgTxt("u not a column vector, or rows(u) != cols(X_skel)"); return;
  }
  if(!(plhs[0]=mxCreateSparse(T_skel.m,T_skel.n,T_skel.jc[T_skel.n],mxREAL)))
    { mexWarnMsgTxt("Out of memory."); return; }
  T.m=T_skel.m, T.n=T_skel.n;
  T.ir = mxGetIr(plhs[0]), T.jc = mxGetJc(plhs[0]), T.pr = mxGetPr(plhs[0]);
  memcpy(T.ir,T_skel.ir,T_skel.jc[T.n]*sizeof(mwIndex));
  memcpy(T.jc,T_skel.jc,(T.n+1)*sizeof(mwIndex));
  interp_lmop(&T,&A,u,&X_skel);
}
