  z = (double **) bss_malloc((k+1)*sizeof(double *));

  for (i=1; i<=k; i++)
    {
      z[i] = (double *) bss_malloc((k+1)*REAL_LEN);
      rvec_zero(z[i],k+1);
      z[i][i] = 1.0;
    }


double *
tri_calc_evec(double *d, double *sd, int n, double lambda)
{
  int i, off;
  double *vec, *tri;
  double beta;


  /* allocate space for eigenvecto and tridiag matrix */
  vec = (double *) bss_malloc((n+1)*REAL_LEN);
  rvec_zero(vec,(n+1));
  tri = (double *) bss_malloc(3*n*REAL_LEN);
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
      printf("%d  %.16e\n",i,*(tri+off+1));
      printf("%d  %.16e\n",i,*(tri+off+4));
      beta = *(tri+off+1);
      if (fabs(beta)<1.0e-10)
	{printf("%d  %.16e\n",i,beta);}
      *(tri+off+2) /= beta;
      *(tri+off+4) -= *(tri+off+2) * *(tri+off+3);
      if (i==(n-2))
	{printf("%d  %.16e\n",n,*(tri+off+4));}
    }

  vec[n-2]=*(tri+3*(n-2)+2);
  vec[n-1]=-1.0;

  /* up */
  for(i=(n-2);i>0;i--)
    {*(vec+i-1) = -*(tri+3*i-1) * *(vec+i);}

  normalize(vec,n);
  for(i=n;i>0;i--)
    {vec[i]=-1.0*vec[i-1];}


  bss_free(tri);
  return(vec);
}






  int *sm_in, *sm_out;

  /* use smbwr.c to re-order initial matrix */
  sm_in  = extract_od(m->adj_list,lda);

  sm_bandwidth_reduction(sm_in,iptr,-1);
  for (i=0; i<lda; i++)
    {printf("%d ", iptr[i]); iptr[i] -=1;}
  printf("\n");



  /*

  rvec_set(rhs,1.0,n);

  for (i=0; i<n; i++)
    {
      denom = sin(((2.0*i-n)/(double) n)*PI/2.0);
      *tmp += denom;
      tmp++;
    }



  Random_Seed(1011);
  for (i=0; i<n; i++)
    {
      denom = Random_Integer(0,n-1);
      denom /= n;
      *tmp-= denom;
      tmp++;
    }

  root_n = (int) sqrt(n);

  if (root_n*root_n > n)
    {root_n--;}

  if (root_n*root_n != n)
    {error_msg_warning("Init_rhs() :: not square!");}

  if (root_n*root_n > n)
    {error_msg_fatal("Init_rhs() :: overshoot!");}

  for (i=0; i<root_n; i++)
    {
      for (j=0; j<root_n; j++)
	{
	  offset = abs(j-root_n/2);
	  denom = 1.0/(1.0*(offset+1));
	  denom /= n;
	  *tmp-= denom;
	  tmp++;
	}
    }



  return(rhs);



  ortho(rhs,n);


  double *rhs;


  rhs = (double *) bss_malloc(n*REAL_LEN);

  rvec_set(rhs,1.0/sqrt(n),n);

  return(rhs);
  */
double *v_add(double *x1, double *x2, int n);
double *v_const_mult(double x1, double *x2, int n);





/********************************sparse_matrix.c*******************************
Function: v_add()

Input : 
Output: 
Return: 
Description:   adds two vectors: x1+x2
*********************************sparse_matrix.c******************************/
double *v_add(double *x1, double *x2, int n)
{
  double *sum, *tmp;


  sum = tmp = (double *) bss_malloc(n*REAL_LEN);

  while (n--) {*tmp++ = *x1++ + *x2++;}

  return(sum);             
}



/********************************sparse_matrix.c*******************************
Function: v_const_mult()

Input : 
Output: 
Return: 
Description:   multiplies a vector by a real constant: alph*ax2_i.
*********************************sparse_matrix.c******************************/
double *v_const_mult(double alpha, double *x2, int n)
{
  double *sca, *tmp;


  sca = tmp = (double *) bss_malloc(n*REAL_LEN);

  while (n--) {*tmp++ = alpha * *x2++;}

  return(sca);             
}





        if (k==1)
      {r_k1 = r_k2 = p_k = p_k1 = r_k;}
    else
      {
	tp = r_k2;
	r_k2 = r_k1;
	r_k1 = r_k;
	if (k>3) {bss_free(tp);}
	tp = p_k1;
	p_k1=p_k;
	if (k>3) {bss_free(tp);}
	tp=x_k1;
	x_k1=x_k;
	if (k>2) {bss_free(tp);}
	beta = v_mult(r_k1,r_k1,n)/v_mult(r_k2,r_k2,n);
	tpcm = v_const_mult(beta,p_k1,n);
	p_k  = v_add(r_k1,tpcm,n);
	bss_free(tpcm);
	ortho(p_k,n);
      }

    alpha = v_mult(r_k1,r_k1,n)/v_mult(p_k,tpmm,n);
    tpcm=v_const_mult(alpha,p_k,n);
    x_k = v_add(x_k1,tpcm,n);
    bss_free(tpcm);
    tpcm=v_const_mult(-1.0*alpha,tpmm,n);
    r_k=v_add(r_k1,tpcm,n);
    bss_free(tpcm);
    bss_free(tpmm);

    res=v_mult(r_k,r_k,n);




  printf("\n{");
  for (i=0; i<k; i++)
    {
      printf("\n{");
      for (j=0; j<k; j++)
	{
	  if (i==j)
	    {printf("%.5lf",dv[i]); *(z+i*k+j) = dv[i];}
	  else if (i==(j-1))
	    {printf("%.5lf",ev[i+1]); *(z+i*k+j) = ev[i+1];}
	  else if (i==(j+1))
	    {printf("%.5lf",ev[i]); *(z+i*k+j) = ev[i];}
	  else 
	    {printf("%.5lf",0.0); *(z+i*k+j) = 0.0;}

	  if (j!=k-1)
	    {printf(",");}
	}
      printf("}");
      if (i!=k-1)
	{printf(",");}
    }
  printf("\n}\n");



  i=j=k;
  gaujord_ (z,&i,&j,dv_w);



  for (i=0; i<k; i++)
    {
      for (j=0; j<k; j++)
	{
	  printf("%.5lf", *(z+i*k+j));
	  if (j!=k-1)
	    {printf(", ");}
	}
      printf("\n");
    }

