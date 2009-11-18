#ifndef TENSOR_H
#define TENSOR_H

#if !defined(TYPES_H)
#warning "tensor.h" requires "types.h"
#endif

#define tensor_dot  PREFIXED_NAME(tensor_dot )
#define tensor_mrxv PREFIXED_NAME(tensor_mrxv)
#define tensor_i1   PREFIXED_NAME(tensor_i1  )
#define tensor_i2   PREFIXED_NAME(tensor_i2  )
#define tensor_i3   PREFIXED_NAME(tensor_i3  )
#define tensor_ig1  PREFIXED_NAME(tensor_ig1 )
#define tensor_ig2  PREFIXED_NAME(tensor_ig2 )
#define tensor_ig3  PREFIXED_NAME(tensor_ig3 )

double tensor_dot(const double *a, const double *b, uint n);

/* y = A x      where A is in row-major format */
void tensor_mrxv(double *y, uint ny, const double *A, const double *x, uint nx);

/*--------------------------------------------------------------------------
   1-,2-,3-d Tensor Application of Row Vectors (for Interpolation)
   
   the 3d case:
   v = tensor_i3(Jr,nr, Js,ns, Jt,nt, u, work)
     gives v = [ Jr (x) Js (x) Jt ] u
     where Jr, Js, Jt are row vectors (interpolation weights)
     u is nr x ns x nt in column-major format (inner index is r)
     v is a scalar
  --------------------------------------------------------------------------*/

double tensor_i1(const double *Jr, uint nr, const double *u);

/* work holds ns doubles */
double tensor_i2(const double *Jr, uint nr,
                 const double *Js, uint ns,
                 const double *u, double *work);

/* work holds ns*nt + nt doubles */
double tensor_i3(const double *Jr, uint nr,
                 const double *Js, uint ns,
                 const double *Jt, uint nt,
                 const double *u, double *work);

/*--------------------------------------------------------------------------
   1-,2-,3-d Tensor Application of Row Vectors
             for simultaneous Interpolation and Gradient computation
   
   the 3d case:
   v = tensor_ig3(Jr,Dr,nr, Js,Ds,ns, Jt,Dt,nt, u,g, work)
     gives v   = [ Jr (x) Js (x) Jt ] u
           g_0 = [ Dr (x) Js (x) Jt ] u
           g_1 = [ Jr (x) Ds (x) Jt ] u
           g_2 = [ Jr (x) Js (x) Dt ] u
     where Jr,Dr,Js,Ds,Jt,Dt are row vectors
       (interpolation & derivative weights)
     u is nr x ns x nt in column-major format (inner index is r)
     v is a scalar, g is an array of 3 doubles
  --------------------------------------------------------------------------*/

double tensor_ig1(double g[1],
                  const double *Jr, const double *Dr, uint nr,
                  const double *u);

/* work holds 2*ns doubles */
double tensor_ig2(double g[2],
                  const double *Jr, const double *Dr, uint nr,
                  const double *Js, const double *Ds, uint ns,
                  const double *u, double *work);

/* work holds 2*ns*nt + 3*ns doubles */
double tensor_ig3(double g[3],
                  const double *Jr, const double *Dr, uint nr,
                  const double *Js, const double *Ds, uint ns,
                  const double *Jt, const double *Dt, uint nt,
                  const double *u, double *work);

#endif
