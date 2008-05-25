#include <stdio.h>
#include <float.h>
#include "types.h"
#include "poly.h"

int main()
{
  int i, n=4;
  real z[50], w[50];
  lobatto_nodes(z,n);
  lobatto_weights(z,w,n);
  /*
  for(i=0;i<n;++i)
    printf("%+20.*Lg\t%+20.*Lg\n",LDBL_DIG+1,(long double)z[i],
                                  LDBL_DIG+1,(long double)w[i]);
  */
  for(i=0;i<n;++i)
    printf("%+20.*g\t%+20.*g\n",DBL_DIG+1,(double)z[i],
                                DBL_DIG+1,(double)w[i]);
  return 0;
}

