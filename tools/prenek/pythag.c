#include <math.h>
#define NRANSI
/* #include "nrutil.h" */

static double sqrarg;

#define SQR(a)   ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software (-4. */

