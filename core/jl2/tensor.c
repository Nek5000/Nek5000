#include "name.h"
#include "types.h"

#define tensor_dot  PREFIXED_NAME(tensor_dot )
#define tensor_mrxv PREFIXED_NAME(tensor_mrxv)
#define tensor_i1   PREFIXED_NAME(tensor_i1  )
#define tensor_i2   PREFIXED_NAME(tensor_i2  )
#define tensor_i3   PREFIXED_NAME(tensor_i3  )
#define tensor_ig1  PREFIXED_NAME(tensor_ig1 )
#define tensor_ig2  PREFIXED_NAME(tensor_ig2 )
#define tensor_ig3  PREFIXED_NAME(tensor_ig3 )

double tensor_dot(const double *a, const double *b, uint n)
{
  double sum = 0;
  for(;n;--n) sum += *a++ * *b++;
  return sum;
}

/* y = A x      where A is in row-major format */
void tensor_mrxv(double *y, uint ny, const double *A,
              const double *x, uint nx)
{
  for(;ny;--ny) {
    const double *xp = x;
    uint n = nx;
    double sum = *A++ * *xp++;
    for(--n;n;--n) sum += *A++ * *xp++;
    *y++ = sum;
  }
}


/*--------------------------------------------------------------------------
   1-,2-,3-d Tensor Application of Row Vectors (for Interpolation)
   
   the 3d case:
   v = tensor_i3(Jr,nr, Js,ns, Jt,nt, u, work)
     gives v = [ Jr (x) Js (x) Jt ] u
     where Jr, Js, Jt are row vectors (interpolation weights)
     u is nr x ns x nt in column-major format (inner index is r)
     v is a scalar
  --------------------------------------------------------------------------*/

double tensor_i1(const double *Jr, uint nr, const double *u)
{
  return tensor_dot(Jr,u,nr);
}

/* work holds ns doubles */
double tensor_i2(const double *Jr, uint nr,
                 const double *Js, uint ns,
                 const double *u, double *work)
{
  tensor_mrxv(work,ns, u, Jr,nr);
  return tensor_dot(Js,work,ns);
}

/* work holds ns*nt + nt doubles */
double tensor_i3(const double *Jr, uint nr,
                 const double *Js, uint ns,
                 const double *Jt, uint nt,
                 const double *u, double *work)
{
  double *work2 = work+nt;
  tensor_mrxv(work2,ns*nt,   u,     Jr,nr);
  tensor_mrxv(work ,nt   ,   work2, Js,ns);
  return tensor_dot(Jt,work,nt);
}

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
                  const double *u)
{
  *g =   tensor_dot(Dr,u,nr);
  return tensor_dot(Jr,u,nr);
}

/* work holds 2*ns doubles */
double tensor_ig2(double g[2],
                  const double *Jr, const double *Dr, uint nr,
                  const double *Js, const double *Ds, uint ns,
                  const double *u, double *work)
{
  double *a = work, *ar = a+ns;
  tensor_mrxv(a ,ns, u, Jr,nr);
  tensor_mrxv(ar,ns, u, Dr,nr);
  g[0] = tensor_dot(Js,ar,ns);
  g[1] = tensor_dot(Ds,a ,ns);
  return tensor_dot(Js,a ,ns);
}

/* work holds 2*ns*nt + 3*ns doubles */
double tensor_ig3(double g[3],
                  const double *Jr, const double *Dr, uint nr,
                  const double *Js, const double *Ds, uint ns,
                  const double *Jt, const double *Dt, uint nt,
                  const double *u, double *work)
{
  const uint nsnt = ns*nt;
  double *a = work, *ar = a+nsnt, *b = ar+nsnt, *br = b+ns, *bs = br+ns;
  tensor_mrxv(a ,nsnt, u , Jr,nr);
  tensor_mrxv(ar,nsnt, u , Dr,nr);
  tensor_mrxv(b ,nt  , a , Js,ns);
  tensor_mrxv(br,nt  , ar, Js,ns);
  tensor_mrxv(bs,nt  , a , Ds,ns);
  g[0] = tensor_dot(Jt,br,nt);
  g[1] = tensor_dot(Jt,bs,nt);
  g[2] = tensor_dot(Dt,b ,nt);
  return tensor_dot(Jt,b ,nt);
}
