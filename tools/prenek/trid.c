
double *
tri_calc_evec(double *d, double *sd, int n, double lambda)
{
  int i, off;
  double *vec, *tri;
  double beta;


  /* allocate space for eigenvecto and tridiag matrix */
  vec = (double *) bss_malloc((n+1)*FLOAT_LEN);
  rvec_zero(vec,(n+1));
  tri = (double *) bss_malloc(3*n*FLOAT_LEN);
  rvec_zero(tri,3*n);

  /* create tri from diag and sub-diag */
  off = 0;
  *(tri+off++) = 0.0;
  *(tri+off++) = *(d+1) - lambda;
  *(tri+off)   = *(sd+2);
  for(i=1; i<(n-1); i++)
    {
      off = 3*i;
      *(tri+off++) = *(sd+i+1);
      *(tri+off++) = *(d+i+1) - lambda;
      *(tri+off)   = *(sd+i+2);
    }
  off = 3*(n-1);
  *(tri+off++) = *(sd+n);
  *(tri+off++) = *(d+n) - lambda;
  *(tri+off)   = 0.0;

  /* Down */
  for(i=0;i<(n-1);i++)
    {
      off = 3*i;
      beta = *(tri+off+1);
      if (fabs(beta)<1.0e-10)
	{printf("%d  %.16e\n",i,beta);}
      *(tri+off+4) -= *(tri+off+2) * *(tri+off+3)/beta;

    }

  for(i=0;i<(n-1);i++)
    {
      off = 3*i;
      *(tri+off+2) /= *(tri+off+1);
    }

  vec[n-2]=*(tri+3*(n-2)+2);
  vec[n-1]=-1.0;

  /* up */
  for(i=(n-2);i>0;i--)
    {*(vec+i-1) = -*(tri+3*i-1) * *(vec+i);}

  normalize(vec,n);
  for(i=n;i>0;i--)
    {vec[i]=-1.0*vec[i-1];}
  vec[0]=0.0;

  bss_free(tri);
  return(vec);
}


void
trid_slv(double *d, double *sd, int n, double lambda, double *b)
{
  int i, off;
  double *vec, *tri;
  double beta;


  for(i=0;i<n;i++)
    {b[i]=b[i+1];}
  b[n]=0.0;

  /* allocate space for eigenvector and tridiag matrix */
  tri = (double *) bss_malloc(3*n*FLOAT_LEN);
  rvec_zero(tri,3*n);

  /* create tri from diag and sub-diag */
  off = 0;
  *(tri+off++) = 0.0;
  *(tri+off++) = *(d+1) - lambda;
  *(tri+off)   = *(sd+2);
  for(i=1; i<(n-1); i++)
    {
      off = 3*i;
      *(tri+off++) = *(sd+i+1);
      *(tri+off++) = *(d+i+1) - lambda;
      *(tri+off)   = *(sd+i+2);
    }
  off = 3*(n-1);
  *(tri+off++) = *(sd+n);
  *(tri+off++) = *(d+n) - lambda;
  *(tri+off)   = 0.0;


/* Down */
for(i=0;i<(n-1);i++)
  {
     *(tri+3*i+2) /= *(tri+3*i+1);
     *(b+i) /= *(tri+3*i+1);
     *(b+i+1) += -*(tri+3*i+3)* *(b+i);
     *(tri+3*i+4) += -*(tri+3*i+3)* *(tri+3*i+2);
   }

/* last entry */     
*(b+n-1) /= *(tri+3*(n-1)+1);

/* up */
for(i=(n-1);i>0;i--)
  {*(b+i-1) += -*(tri+3*i-1) * *(b+i);}


normalize(b,n);
for(i=n;i>0;i--)
  {b[i]=b[i-1];}
b[0]=0.0;

free(tri);
}
