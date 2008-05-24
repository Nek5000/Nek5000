#include <math.h>
#define NRANSI
/*#include "nrutil.h"*/

#define SIGN(a,b)  ((b) >= 0.0 ? fabs(a) : -fabs(a))

void tqli(double d[], double e[], int n, double **z)
{
  double pythag(double a, double b);
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  
  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) 
	  {error_msg_fatal("Too many iterations in tqli");}
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=pythag(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	  /* calc eigenvectors iff z=I passed in */
	  if (z)
	    {
	      for (k=1;k<=n;k++) {
		f=z[k][i+1];
		z[k][i+1]=s*z[k][i]+c*f;
		z[k][i]=c*z[k][i]-s*f;
	      }
	    }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}

#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software (-4. */
