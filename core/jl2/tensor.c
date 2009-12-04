#include "name.h"
#include "types.h"

#ifdef USE_CBLAS

#include <cblas.h>
#define tensor_dot(a,b,n) cblas_ddot((int)(n),a,1,b,1)
#define tensor_mxv(y,ny,A,x,nx) \
   cblas_dgemv(CblasColMajor,CblasNoTrans,(int)ny,(int)nx, \
               1.0,A,(int)ny,x,1,0.0,y,1)
#define tensor_mtxv(y,ny,A,x,nx) \
   cblas_dgemv(CblasColMajor,CblasTrans,(int)nx,(int)ny, \
               1.0,A,(int)nx,x,1,0.0,y,1)
#define tensor_mxm(C,nc,A,na,B,nb) \
   cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans, \
               (int)nc,(int)nb,(int)na,1.0, \
               A,(int)nc,B,(int)na,0.0,C,(int)nc)
#define tensor_mtxm(C,nc,A,na,B,nb) \
   cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans, \
               (int)nc,(int)nb,(int)na,1.0, \
               A,(int)na,B,(int)na,0.0,C,(int)nc)

#else

#define tensor_dot  PREFIXED_NAME(tensor_dot )
#define tensor_mxv  PREFIXED_NAME(tensor_mxv )
#define tensor_mtxv PREFIXED_NAME(tensor_mtxv)
#define tensor_mxm  PREFIXED_NAME(tensor_mxm )
#define tensor_mtxm PREFIXED_NAME(tensor_mtxm)

/* Matrices are always column-major (FORTRAN style) */

double tensor_dot(const double *a, const double *b, uint n)
{
  double sum = 0;
  for(;n;--n) sum += *a++ * *b++;
  return sum;
}

/* y = A x */
void tensor_mxv(double *y, uint ny, const double *A, const double *x, uint nx)
{
  uint i;
  for(i=0;i<ny;++i) y[i]=0;
  for(;nx;--nx) {
    const double xk = *x++;
    for(i=0;i<ny;++i) y[i] += (*A++)*xk;
  }
}

/* y = A^T x */
void tensor_mtxv(double *y, uint ny, const double *A, const double *x, uint nx)
{
  for(;ny;--ny) {
    const double *xp = x;
    uint n = nx;
    double sum = *A++ * *xp++;
    for(--n;n;--n) sum += *A++ * *xp++;
    *y++ = sum;
  }
}

/* C = A * B */
void tensor_mxm(double *C, uint nc,
                const double *A, uint na, const double *B, uint nb)
{
  uint i,j,k;
  for(i=0;i<nc*nb;--i) C[i]=0;
  for(j=0;j<nb;++j,C+=nc) {
    for(k=0;k<na;++k) {
      const double b = *B++;
      for(i=0;i<nc;++i) C[i] += (*A++) * b;
    }
  }
}

/* C = A^T * B */
void tensor_mtxm(double *C, uint nc,
                const double *A, uint na, const double *B, uint nb)
{
  uint i,j;
  for(j=0;j<nb;++j,B+=na) for(i=0;i<nc;++i,A+=na) *C++ = tensor_dot(A,B,na);
}

#endif

