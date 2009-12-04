#ifndef OBBOX_H
#define OBBOX_H

#if !defined(TYPES_H) || !defined(NAME_H)
#warning "obbox.h" requires "types.h" and "name.h"
#endif

#define obbox_2_calc  PREFIXED_NAME(obbox_2_calc)
#define obbox_3_calc  PREFIXED_NAME(obbox_3_calc)

/*--------------------------------------------------------------------------
   Oriented and axis-aligned bounding box computation for spectral elements
   
   Usage:
   
     double x[n][nt][ns][nr], y[n][nt][ns][nr], z[n][nt][ns][nr];
     obbox_3 ob[n];

     unsigned mr=4*nr, ms=4*ns, mt=4*nt;
     double tol = 1e-6;
     obbox_3_calc(ob, x,y,z, nr,ns,nt,n, mr,ms,mt, tol);
     
   The parameters mr,ms,mt specify number of points to use in computing
   bounds (see lob_bnd.h). It is expected that mr>nr, etc. For reasonable
   quality, a factor of at least 2 is recommended.
  
   tol is a relative amount by which to expand the bounding box.
   This would accommodate, e.g., rounding errors.
  
   The axis aligned bounds for a given element are
     ob[i].x.min <= x <= ob[i].x.max
     ob[i].y.min <= y <= ob[i].y.max
     ob[i].z.min <= z <= ob[i].z.max

   The oriented bounding box is given by
     (-1,-1,-1)^T <= ob[i].A * (x - ob[i].c0) <= (1,1,1)
   
   where the matrix is row-major format,
     dx = x - c0[0], dy = y - c0[1], dz = z - c0[2]
     -1 <= r[0] = A[0]*dx + A[1]*dy + A[2]*dz <= 1
     -1 <= r[1] = A[3]*dx + A[4]*dy + A[5]*dz <= 1
     -1 <= r[2] = A[6]*dx + A[7]*dy + A[8]*dz <= 1

   Also, ob[i].A * (x - ob[i].c0) should be a reasonable seed for Newton's.
    
  --------------------------------------------------------------------------*/

#ifndef LOB_BND_H
struct dbl_range { double min, max; };
#endif

typedef struct { double c0[2], A[4];
                 struct dbl_range x,y; } obbox_2;

typedef struct { double c0[3], A[9];
                 struct dbl_range x,y,z; } obbox_3;

void obbox_2_calc(obbox_2 *out,
  const double *x, const double *y,
  unsigned nr, unsigned ns, uint n,
  unsigned mr, unsigned ms, double tol);

void obbox_3_calc(obbox_3 *out,
  const double *x, const double *y, const double *z,
  unsigned nr, unsigned ns, unsigned nt, uint n,
  unsigned mr, unsigned ms, unsigned mt, double tol);

#endif

