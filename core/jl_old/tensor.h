#ifndef TENSOR_H
#define TENSOR_H

/* requires "types.h" */

#ifndef TYPES_H
#warning "tensor.h" requires "types.h"
#endif

/*--------------------------------------------------------------------------
   1-,2-,3-d Tensor Application
   
   the 3d case:
   tensor_f3(R,mr,nr, S,ms,ns, T,mt,nt, u,v, work1,work2)
     gives v = [ R (x) S (x) T ] u
     where R is mr x nr, S is ms x ns, T is mt x nt,
       each in row- or column-major format according to f := r | c
     u is nr x ns x nt in column-major format (inner index is r)
     v is mr x ms x mt in column-major format (inner index is r)
  --------------------------------------------------------------------------*/

void tensor_c1(const real *R, unsigned mr, unsigned nr, 
               const real *u, real *v);
void tensor_r1(const real *R, unsigned mr, unsigned nr, 
               const real *u, real *v);

/* work holds mr*ns reals */
void tensor_c2(const real *R, unsigned mr, unsigned nr,
               const real *S, unsigned ms, unsigned ns,
               const real *u, real *v, real *work);
void tensor_r2(const real *R, unsigned mr, unsigned nr,
               const real *S, unsigned ms, unsigned ns,
               const real *u, real *v, real *work);

/* work1 holds mr*ns*nt reals,
   work2 holds mr*ms*nt reals */
void tensor_c3(const real *R, unsigned mr, unsigned nr,
               const real *S, unsigned ms, unsigned ns,
               const real *T, unsigned mt, unsigned nt,
               const real *u, real *v, real *work1, real *work2);
void tensor_r3(const real *R, unsigned mr, unsigned nr,
               const real *S, unsigned ms, unsigned ns,
               const real *T, unsigned mt, unsigned nt,
               const real *u, real *v, real *work1, real *work2);

/*--------------------------------------------------------------------------
   1-,2-,3-d Tensor Application of Row Vectors (for Interpolation)
   
   the 3d case:
   v = tensor_i3(Jr,nr, Js,ns, Jt,nt, u, work)
   same effect as tensor_r3(Jr,1,nr, Js,1,ns, Jt,1,nt, u,&v, work1,work2):
     gives v = [ Jr (x) Js (x) Jt ] u
     where Jr, Js, Jt are row vectors (interpolation weights)
     u is nr x ns x nt in column-major format (inner index is r)
     v is a scalar
  --------------------------------------------------------------------------*/

real tensor_i1(const real *Jr, unsigned nr, const real *u);

/* work holds ns reals */
real tensor_i2(const real *Jr, unsigned nr,
               const real *Js, unsigned ns,
               const real *u, real *work);

/* work holds ns*nt + nt reals */
real tensor_i3(const real *Jr, unsigned nr,
               const real *Js, unsigned ns,
               const real *Jt, unsigned nt,
               const real *u, real *work);

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
     v is a scalar, g is an array of 3 reals
  --------------------------------------------------------------------------*/

real tensor_ig1(const real *Jr, const real *Dr, unsigned nr,
                const real *u, real *g);

/* work holds 2*ns reals */
real tensor_ig2(const real *Jr, const real *Dr, unsigned nr,
                const real *Js, const real *Ds, unsigned ns,
                const real *u, real *g, real *work);

/* work holds 2*ns*nt + 3*ns reals */
real tensor_ig3(const real *Jr, const real *Dr, unsigned nr,
                const real *Js, const real *Ds, unsigned ns,
                const real *Jt, const real *Dt, unsigned nt,
                const real *u, real *g, real *work);

#endif

