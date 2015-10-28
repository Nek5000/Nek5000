/********************************sparse_matrix.c*******************************
Function: lanczos()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int *
8_lanczos(int *list, int n, int *ia, int *ja, double *vals)
{
  int i,j,k=0;
  double *r_k=NULL, *p_k=NULL, *x_k=NULL, *Ap_k=NULL;
  double alpha=0.0, beta=0.0, rr=0.0, rr_limit, new_rr;
  double alpha_old=0.0, beta_old=0.0;
  int *color, *ptr;
  int iter_limit;
  queue_ADT ln_q;
  struct lanczos_node *ln_ptr;
  double  *dv, *ev, *dv_w, *ev_w, dmax, dmin, **z, **Q;
  double *w_k, *tmp, *v_k, *v_k_old;
  double la, lb;
  int *companion;
  

  /* queue to hold lanczos vectors */
  ln_q = new_queue();

  /* initial rhs w/2-norm of ~1.0 and orthogonal to (1,1,...,1)^T */
  r_k = Init_rhs(n);
  p_k = (double *) bss_malloc(n*FLOAT_LEN);
  rvec_copy(p_k,r_k,n);
  rr  = hmt_ddot(r_k,r_k,n);

  /* x_k to initial guess (the 0.0 vector) */
  x_k  = Init_guess(n);

  /* space for Ap_k */
  Ap_k = (double *) bss_malloc(n*FLOAT_LEN);
  rvec_zero(Ap_k,n);
  
  /* set lanczos limits */
  rr_limit = rr*RES_LIMIT;
  iter_limit = MIN(n,CG_LIMIT);

  /* lanczos */
  v_k = (double *) bss_malloc(n*FLOAT_LEN);
  rvec_copy(v_k,r_k,n);
  v_k_old = (double *) bss_malloc(n*FLOAT_LEN);
  tmp     = (double *) bss_malloc(n*FLOAT_LEN);
  w_k     = (double *) bss_malloc(n*FLOAT_LEN);
  rvec_zero(v_k_old,n);
  rvec_zero(w_k,n);
  rvec_zero(tmp,n);

#ifdef DEBUG    
    printf("START LANCZOS :: res=%.16e, res_lim=%.16e,iter_lim=%d\n",rr,
	   rr_limit,iter_limit);
#endif

  while ((rr>rr_limit) && (k<iter_limit))
  {
    k++;


    ln_ptr = (struct lanczos_node *) bss_malloc(sizeof(struct lanczos_node));
    ln_ptr->v = (double *) bss_malloc(n*FLOAT_LEN);
    rvec_copy(ln_ptr->v,v_k,n);
    /*    normalize(ln_ptr->v,n); */


    if (j==1)
      {matrix_mult(w_k,v_k,n,ia,ja,vals);}
    else
      {
	matrix_mult(w_k,v_k,n,ia,ja,vals);
	hmt_daxpy(w_k,-1.0*lb,v_k_old,n);
      }
    la = hmt_ddot(w_k,v_k,n);
    hmt_daxpy(w_k,-1.0*la,v_k,n);
    lb = sqrt(hmt_ddot(w_k,w_k,n));


    rvec_copy(v_k_old,v_k,n);
    rvec_zero(v_k,n);
    hmt_daxpy(v_k,1.0/lb,w_k,n);


    matrix_mult(Ap_k,p_k,n,ia,ja,vals);
    alpha = rr/hmt_ddot(Ap_k,p_k,n);
    hmt_daxpy(x_k,alpha,p_k,n);
    hmt_daxpy(r_k,-1.0*alpha,Ap_k,n);
    new_rr = hmt_ddot(r_k,r_k,n);
    beta = new_rr/rr;
    /* should scale and then daxpy */
    rev_daxpy(r_k,beta,p_k,n);
    ortho(p_k,n);
    printf("B:norm p_k=%.16e\n",hmt_ddot(p_k,p_k,n));
    /* should we reorthogonalize the p_k's */
    /*
    normalize(p_k,n);
    printf("A:norm p_k=%.16e/n",hmt_ddot(p_k,p_k,n));
    */
    rr = new_rr;

    /* tridiagonal matrix calculation */

    if (k==1)
      {
	ln_ptr->d = 1.0/alpha;
	ln_ptr->e = 0.0;
	printf("d(%d)=%.16e\n",k,ln_ptr->d);
	printf("l(%d)=%.16e\n",k-1,ln_ptr->e);
      }
    else
      {
	ln_ptr->d = 1.0/alpha + beta_old/alpha_old;
	ln_ptr->e = sqrt(beta_old)/alpha_old;
	printf("d(%d)=%.16e\n",k,ln_ptr->d);
	printf("l(%d)=%.16e\n",k-1,ln_ptr->e);
      }

    enqueue(ln_q,ln_ptr);

    alpha_old = alpha;
    beta_old = beta;

    if (lb<=0.00000000001)
      {error_msg_warning("lanczos beta is zero!\n"); break;}

#ifdef DEBUG
    printf("iter#%d :: res=%.16e\n",k,rr);
#endif
  }

  /*
  bss_free(r_k);
  bss_free(p_k);
  */

#ifdef DEBUG
  /* how good is the solution? Check infinity norm. */
  matrix_mult(Ap_k,x_k,n,ia,ja,vals);
  r_k = Init_rhs(n);
  alpha=0.0;
  for (i=0; i<n; i++)
    {alpha = MAX(fabs(Ap_k[i]-r_k[i]),alpha);}

  printf("CHK CG :: inf-norm = %.16e\n",alpha);
  /*  bss_free(r_k); */
#endif

  dv   = (double *) bss_malloc((k+1)*FLOAT_LEN);
  dv_w = (double *) bss_malloc((k+1)*FLOAT_LEN);
  ev   = (double *) bss_malloc((k+1)*FLOAT_LEN);
  ev_w = (double *) bss_malloc((k+1)*FLOAT_LEN);
  Q    = (double **) bss_malloc((k+1)*sizeof(double *));
  Q[0] = (double *) bss_malloc((k+1)*FLOAT_LEN);
  rvec_zero(Q[0],n);
  for (i=1; i<=k; i++)
    {
      ln_ptr = dequeue(ln_q);
      Q[i] = ln_ptr->v;
      dv[i]  = ln_ptr->d;
      if (i>1)
	{
	  ev[i]  = ln_ptr->e;
	  /*	  printf("ev[%d]=%.16e\n",i,ev[i]); */
	}
    }

  if (len_queue(ln_q))
    {error_msg_warning("lanczos queue not exhasuted!");}

  dv[0]=ev[0]=ev[1]=0.0;

  /*  
  calc_ (dv,ev,dv_w,ev_w,&k,&dmax,&dmin,z); 
  printf("eval_min=%.16e\neval_max=%.16e\n",dmin,dmax);
  */
  z = (double **) bss_malloc((k+1)*sizeof(double *));

  for (i=0; i<=k; i++)
    {
      z[i] = (double *) bss_malloc((k+1)*FLOAT_LEN);
      rvec_zero(z[i],k+1);
      (z[i])[i] = 1.0;
    }

  tqli(dv,ev,k,z);

  beta = 99999999.0;
  j=-1;
  for (i=1; i<=k; i++)
    {
      beta_old = MIN(beta,dv[i]);
      if (beta_old != beta)
	{j=i; beta=beta_old;}
      printf("eval[%d]=%.16e\n",i,dv[i]);
    }
  printf("\n");

  printf("min eval[%d]=%.16e\n",j,dv[j]);
  beta = dv[j];

  for (i=1; i<=k; i++)
    {ev_w[i] = (z[i])[j]; printf("%.16e\n",ev_w[i]);}

  rvec_zero(dv_w,n);
  for (j=1; j<=k; j++)
    for (i=0; i<n; i++)
      {dv_w[i] += ev_w[j]*(Q[j])[i];}

  /*check A evec_min = eval_min evec_min */
  matrix_mult(Ap_k,dv_w,n,ia,ja,vals);
  printf("BETA=%.16e\n",beta);
  for (i=0; i<n; i++)
    {
      /*  printf("%.16e %.16e\n",beta*dv_w[i],Ap_k[i]); */
      if ((alpha=fabs(beta*dv_w[i]-Ap_k[i])) > 0.0)
	{printf("%d :: %.16e\n",i,alpha);}
    }

  printf("approximate eigenvector::\n");
  for (i=0; i<n; i++)
    {printf("%d :: %.16e\n",list[i]+1,dv_w[i]);}

  /*
  for (i=1; i<=k; i++)
    {
      for (j=1; j<=k; j++)
	{
	  printf("%.5lf", (z[i])[j]);
	  if (j!=k)
	    {printf(", ");}
	}
      printf("\n");
    }
  */

  /*  solution held in x_k !!!*/
  /*
  bss_free(x_k);
  bss_free(Ap_k);
  */

  /* end CG/Lanczos */
  printf("%d :: %d :: %d\n",n,sizeof(int *), sizeof(double *));

  color = ptr = (int *) bss_malloc(n*INT_LEN);
  companion = (int *) bss_malloc(n*INT_LEN);
  ivec_c_index(companion,n);
  rvec_sort_companion(dv_w,companion,n);

  printf("sorted approximate eigenvector::\n");
  for (i=0; i<n; i++)
    {printf("%d :: %.16e\n",companion[i],dv_w[i]);}

  printf("(");
  for (i=0; i<n/2; i++)
    {
      /*      printf("%d ",list[companion[i]]+1); */
      color[companion[i]] = -1;
    }
  /*  printf(") %d (",list[companion[n/2]]+1); */

  color[companion[n/2]] = 0;

  for (i=n/2; i<n; i++)
    {
      /*      printf("%d ",list[companion[i]]+1); */
      color[companion[i]] = 1;
    }
  printf(")\n");

#ifdef DEBUG_0
  printf("(");
  for (i=0; i<n/2; i++)
    {
      printf("%d ",list[i]+1);
      *ptr++ = -1;
    }
  printf(") %d (",list[n/2]+1);

  *ptr++ = 0;

  for (i=n/2+1; i<n; i++)
    {
      printf("%d ",list[i]+1);
      *ptr++ = 1;
    }
  printf(")\n");
#endif

  return(color);
}



/********************************sparse_matrix.c*******************************
Function: lanczos()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int *
old_lanczos(int *list, int n, int *ia, int *ja, double *vals)
{
  int *color;
  int *ptr, i;

  color = ptr = (int *) bss_malloc(INT_LEN * n);

  /*
  Random_Seed(7);
  for (i=0; i<n; i++)
    {
      printf("%d ",list[i]+1);
      *ptr++ = Random_Integer(-1,1);
    }
  printf("\n");
  */

  printf("(");
  for (i=0; i<n/2; i++)
    {
      printf("%d ",list[i]+1);
      *ptr++ = -1;
    }
  printf(") %d (",list[n/2]+1);

  *ptr++ = 0;

  for (i=n/2+1; i<n; i++)
    {
      printf("%d ",list[i]+1);
      *ptr++ = 1;
    }
  printf(")\n");

  return(color);

	  /* call lanczos on submatrix */
	  /* separate ev about median :: 1/2->lc :: 1/2->rc */
	  /* determine pot. separators */
	  /* extract crs submatrix corresponding to pot. separators */
	  /* use perp. ev to clean up separators */
	  /* special cleaning :: KL? */
	  /* label separators */
}


/********************************sparse_matrix.c*******************************
Function: lanczos()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int *
?_lanczos(int *list, int n, int *ia, int *ja, double *vals)
{
  int i,j,k=0;
  double *r_k=NULL, *p_k=NULL, *x_k=NULL, *Ap_k=NULL;
  double alpha=0.0, beta=0.0, rr=0.0, rr_limit, new_rr;
  double alpha_old=0.0, beta_old=0.0;
  int *color, *ptr;
  int iter_limit;
  queue_ADT ln_q;
  struct lanczos_node *ln_ptr;
  double  *dv, *ev, *dv_w, *ev_w, dmax, dmin, **z, **Q;
  double *w_k, *tmp, *v_k, *v_k_old;
  double la, lb;

  /* queue to hold lanczos vectors */
  ln_q = new_queue();

  /* initial rhs w/2-norm of ~1.0 and orthogonal to (1,1,...,1)^T */
  r_k = Init_rhs(n);
  p_k = (double *) bss_malloc(n*FLOAT_LEN);
  rvec_copy(p_k,r_k,n);
  rr  = hmt_ddot(r_k,r_k,n);

  /* x_k to initial guess (the 0.0 vector) */
  x_k  = Init_guess(n);

  /* space for Ap_k */
  Ap_k = (double *) bss_malloc(n*FLOAT_LEN);
  rvec_zero(Ap_k,n);
  
  /* set lanczos limits */
  rr_limit = rr*RES_LIMIT;
  iter_limit = MIN(n,CG_LIMIT);

  /* lanczos */
  v_k = (double *) bss_malloc(n*FLOAT_LEN);
  rvec_copy(v_k,r_k,n);
  v_k_old = (double *) bss_malloc(n*FLOAT_LEN);
  tmp     = (double *) bss_malloc(n*FLOAT_LEN);
  w_k     = (double *) bss_malloc(n*FLOAT_LEN);
  rvec_zero(v_k_old,n);
  rvec_zero(w_k,n);
  rvec_zero(tmp,n);

#ifdef DEBUG    
    printf("START LANCZOS :: res=%.16e, res_lim=%.16e,iter_lim=%d\n",rr,
	   rr_limit,iter_limit);
#endif

  while ((rr>rr_limit) && (k<iter_limit))
  {
    k++;


    ln_ptr = (struct lanczos_node *) bss_malloc(sizeof(struct lanczos_node));
    ln_ptr->v = (double *) bss_malloc(n*FLOAT_LEN);
    rvec_copy(ln_ptr->v,v_k,n);
    /*    normalize(ln_ptr->v,n); */


    if (j==1)
      {matrix_mult(w_k,v_k,n,ia,ja,vals);}
    else
      {
	matrix_mult(w_k,v_k,n,ia,ja,vals);
	hmt_daxpy(w_k,-1.0*lb,v_k_old,n);
      }
    la = hmt_ddot(w_k,v_k,n);
    hmt_daxpy(w_k,-1.0*la,v_k,n);
    lb = sqrt(hmt_ddot(w_k,w_k,n));


    rvec_copy(v_k_old,v_k,n);
    rvec_zero(v_k,n);
    hmt_daxpy(v_k,1.0/lb,w_k,n);


    matrix_mult(Ap_k,p_k,n,ia,ja,vals);
    alpha = rr/hmt_ddot(Ap_k,p_k,n);
    hmt_daxpy(x_k,alpha,p_k,n);
    hmt_daxpy(r_k,-1.0*alpha,Ap_k,n);
    new_rr = hmt_ddot(r_k,r_k,n);
    beta = new_rr/rr;
    /* should scale and then daxpy */
    rev_daxpy(r_k,beta,p_k,n);
    ortho(p_k,n);
    printf("B:norm p_k=%.16e\n",hmt_ddot(p_k,p_k,n));
    /* should we reorthogonalize the p_k's */
    /*
    normalize(p_k,n);
    printf("A:norm p_k=%.16e/n",hmt_ddot(p_k,p_k,n));
    */
    rr = new_rr;

    /* tridiagonal matrix calculation */

    if (k==1)
      {
	ln_ptr->d = 1.0/alpha;
	ln_ptr->e = 0.0;
	printf("d(%d)=%.16e\n",k,ln_ptr->d);
	printf("l(%d)=%.16e\n",k-1,ln_ptr->e);
      }
    else
      {
	ln_ptr->d = 1.0/alpha + beta_old/alpha_old;
	ln_ptr->e = sqrt(beta_old)/alpha_old;
	printf("d(%d)=%.16e\n",k,ln_ptr->d);
	printf("l(%d)=%.16e\n",k-1,ln_ptr->e);
      }

    enqueue(ln_q,ln_ptr);

    alpha_old = alpha;
    beta_old = beta;

    if (lb<=0.00000000001)
      {break;}

#ifdef DEBUG
    printf("iter#%d :: res=%.16e\n",k,rr);
#endif
  }

  bss_free(r_k);
  bss_free(p_k);

#ifdef DEBUG
  /* how good is the solution? Check infinity norm. */
  matrix_mult(Ap_k,x_k,n,ia,ja,vals);
  r_k = Init_rhs(n);
  alpha=0.0;
  for (i=0; i<n; i++)
    {alpha = MAX(fabs(Ap_k[i]-r_k[i]),alpha);}

  printf("CHK CG :: inf-norm = %.16e\n",alpha);
  bss_free(r_k);
#endif

  dv   = (double *) bss_malloc((k+1)*FLOAT_LEN);
  dv_w = (double *) bss_malloc((k+1)*FLOAT_LEN);
  ev   = (double *) bss_malloc((k+1)*FLOAT_LEN);
  ev_w = (double *) bss_malloc((k+1)*FLOAT_LEN);
  Q    = (double **) bss_malloc((k+1)*sizeof(double *));
  for (i=1; i<=k; i++)
    {
      ln_ptr = dequeue(ln_q);
      Q[i] = ln_ptr->v;
      dv[i]  = ln_ptr->d;
      if (i>1)
	{
	  ev[i]  = ln_ptr->e;
	  /*	  printf("ev[%d]=%.16e\n",i,ev[i]); */
	}
    }

  dv[0]=ev[0]=ev[1]=0.0;


  /*  
  calc_ (dv,ev,dv_w,ev_w,&k,&dmax,&dmin,z); 
  printf("eval_min=%.16e\neval_max=%.16e\n",dmin,dmax);
  */

  z = (double **) bss_malloc((k+1)*sizeof(double *));

  for (i=0; i<=k; i++)
    {
      z[i] = (double *) bss_malloc((k+1)*FLOAT_LEN);
      rvec_zero(z[i],k+1);
      (z[i])[i] = 1.0;
    }


  tqli(dv,ev,k,z);
  beta = 10000.0;
  j=-1;
  for (i=1; i<=k; i++)
    {
      beta_old = MIN(beta,dv[i]);
      if (beta_old != beta)
	{j=i; beta=beta_old;}
      printf("eval[%d]=%.16e\n",i,dv[i]);
    }
  printf("\n");

  printf("min eval[%d]=%.16e\n",j,dv[j]);
  beta = dv[j];

  for (i=1; i<=k; i++)
    {ev_w[i] = (z[i])[j]; printf("%.16e\n",ev_w[i]);}


  rvec_zero(dv_w,n);
  for (j=1; j<=k; j++)
    for (i=0; i<n; i++)
      {dv_w[i] += ev_w[j]*(Q[j])[i];}

  
  matrix_mult(Ap_k,dv_w,n,ia,ja,vals);

  printf("BETA=%.16e\n",beta);
  for (i=0; i<n; i++)
    {
      printf("%.16e %.16e\n",beta*dv_w[i],Ap_k[i]);
      /*
      if ((alpha=fabs(beta*dv_w[i]-Ap_k[i])) > 0.0)
	printf("%d :: %.16e\n",i,alpha);
     */
    }

  printf("approximate eigenvector::\n");
  for (i=0; i<n; i++)
    {printf("%d :: %.16e\n",list[i]+1,dv_w[i]);}


  /*later check A evec_min = eval_min evec_min */

  for (i=1; i<=k; i++)
    {
      for (j=1; j<=k; j++)
	{
	  printf("%.5lf", (z[i])[j]);
	  if (j!=k)
	    {printf(", ");}
	}
      printf("\n");
    }

  /*  solution held in x_k !!!*/
  bss_free(x_k);
  bss_free(Ap_k);


  /* end CG/Lanczos */
  color = ptr = (int *) bss_malloc(INT_LEN * n);

  printf("(");
  for (i=0; i<n/2; i++)
    {
      printf("%d ",list[i]+1);
      *ptr++ = -1;
    }
  printf(") %d (",list[n/2]+1);

  *ptr++ = 1;

  for (i=n/2+1; i<n; i++)
    {
      printf("%d ",list[i]+1);
      *ptr++ = 1;
    }
  printf(")\n");

  return(color);
}
