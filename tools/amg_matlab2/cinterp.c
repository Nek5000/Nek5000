#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "cinterp_common.h"

static void interp(const sp_mat *X, const sp_mat_c *A, const sp_mat_c *B,
                   const double *u, const double *lambda)
{
  mwSize nf = X->m, nc = X->n;
  mwSize max_nz=0, max_Q;
  mwIndex j;
  
  double *sqv1, *sqv2;
  double *Q, *QQt;
  
  for(j=0;j<nc;++j) {
    mwSize nz=X->jc[j+1]-X->jc[j];
    if(nz>max_nz) max_nz=nz;
  }
  max_Q = (max_nz*(max_nz+1))/2;
  
  sqv1 = mem_alloc((2*max_nz + max_Q)*sizeof(double));
  sqv2 = sqv1+max_nz, Q = sqv2+max_nz;

  for(j=0;j<nc;++j) {
    mwIndex xjc = X->jc[j];
    const mwIndex *Qi = &X->ir[xjc];
    mwSize nz = X->jc[j+1]-xjc;
    mwIndex m,k;
    double *qk = Q;
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
    }
    /* sqv1 := R B e_j */
    sp_restrict_sorted(sqv1, nz,Qi, B->jc[j+1]-B->jc[j],
      &B->ir[B->jc[j]], &B->pr[B->jc[j]]);
    /* sqv1 := R (B e_j + u_j lambda) */
    for(k=0;k<nz;++k) sqv1[k] += u[j]*lambda[Qi[k]];
    /* sqv2 := Q^t (B e_j + u_j lambda) */
    mv_utt(sqv2, nz,Q, sqv1);
    /* X e_j := Q Q^t (B e_j + u_j lambda) */
    mv_ut(&X->pr[xjc], nz,Q, sqv2);
  }
  mem_free(sqv1);
}

/* A, B, u, lambda, X_skel */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  sp_mat X;
  sp_mat_c A, B;
  sp_log_mat_c X_skel;
  const double *u, *lambda;
  if(nrhs!=5) { mexWarnMsgTxt("Five inputs required."); return; }
  if(!mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
    mexWarnMsgTxt("First input not a sparse, double, real matrix."); return;
  }
  if(!mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
    mexWarnMsgTxt("Second input not a sparse, double, real matrix."); return;
  }
  if(mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
    mexWarnMsgTxt("Third input not a full, double, real array."); return;
  }
  if(mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
    mexWarnMsgTxt("Fourth input not a full, double, real array."); return;
  }
  if(!mxIsSparse(prhs[4])) {
    mexWarnMsgTxt("Fifth input not a sparse matrix."); return;
  }
  A.m = mxGetM(prhs[0]), A.n = mxGetN(prhs[0]);
  A.ir = mxGetIr(prhs[0]), A.jc = mxGetJc(prhs[0]), A.pr = mxGetPr(prhs[0]);
  B.m = mxGetM(prhs[1]), B.n = mxGetN(prhs[1]);
  B.ir = mxGetIr(prhs[1]), B.jc = mxGetJc(prhs[1]), B.pr = mxGetPr(prhs[1]);
  u = mxGetPr(prhs[2]);
  lambda = mxGetPr(prhs[3]);
  X_skel.m = mxGetM(prhs[4]), X_skel.n = mxGetN(prhs[4]);
  X_skel.ir = mxGetIr(prhs[4]), X_skel.jc = mxGetJc(prhs[4]);
  if(A.m!=A.n) { mexWarnMsgTxt("A not square."); return; }
  if(A.m!=B.m) { mexWarnMsgTxt("rows(A) != rows(B)"); return; }
  if(A.m!=X_skel.m) { mexWarnMsgTxt("rows(A) != rows(X_skel)"); return; }
  if(B.n!=X_skel.n) { mexWarnMsgTxt("cols(B) != cols(X_skel)"); return; }
  if(mxGetM(prhs[2])!=B.n || mxGetN(prhs[2])!=1) {
    mexWarnMsgTxt("u not a column vector, or rows(u) != cols(B)"); return;
  }
  if(mxGetM(prhs[3])!=A.m || mxGetN(prhs[3])!=1) {
    mexWarnMsgTxt("lambda not a column vector, or rows(lambda) != rows(A)");
    return;
  }
  plhs[0] = mxCreateSparse(X_skel.m,X_skel.n,X_skel.jc[X_skel.n],mxREAL);
  X.m=X_skel.m, X.n=X_skel.n;
  X.ir = mxGetIr(plhs[0]), X.jc = mxGetJc(plhs[0]), X.pr = mxGetPr(plhs[0]);
  memcpy(X.ir,X_skel.ir,X_skel.jc[X.n]*sizeof(mwIndex));
  memcpy(X.jc,X_skel.jc,(X.n+1)*sizeof(mwIndex));
  interp(&X,&A,&B,u,lambda);
}
