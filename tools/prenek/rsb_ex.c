/********************************sparse_matrix.c*******************************
Module Name: sparse_matrix.c
Module Info:

Author:  Henry M. Tufo III

e-mail:  hmt@@cs.brown.edu (lgrind ... so excuse the double <at>)

sn-mail: Division of Applied Mathematics, 
	 Brown University,
         Box F
	 Providence, RI 02912

Tel:	 (401) 863-7666


Last Modification: 3.27.96
*********************************sparse_matrix.c******************************/


/********************************sparse_matrix.c*******************************
NOTES ON USAGE: 

*********************************sparse_matrix.c******************************/


/********************************sparse_matrix.c*******************************
FILE FORMAT: 
------------------------------ Begin File -------------------------------------

------------------------------ End   File -------------------------------------

Note: 
*********************************sparse_matrix.c******************************/

/* C modules for I/O etc. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
/*#include <nx.h>*/


/* mine : const before types! */
#include "const.h"
#include "types.h"
#include "error.h"
#include "ivec.h"
#include "bss_malloc.h"
#include "adj_list.h"
#include "queue.h"
#include "stack.h"
#include "rsb_driver.h"



struct lanczos_node {
  double d;
  double e;
  double *v;
};


typedef struct id{
  int lda, sda, nnz;
  char ty[80];
  char desc[80];
  int  *permutation;
  int  gs_handle;
  adj_node_ptr adj_list;
  int     *row_lens;
  int    **row_indices;
  double **row_vals;
  int     *col_lens;
  int    **col_indices;
  double **col_vals;
} matrix_id;


/* sparse matrix template */
typedef struct id  *matrix_id_ptr;


struct tree_node{
  int depth;
  int rank;
  int num_act;
  int *active;
  int num_sep;
  int *sep;
  REAL *rhs;
  struct tree_node *lc, *rc;
  int lb, ub;
};


/* defines */
#define   MAX_LEN       101             /* max input line length        */
#define   MAX_LDA       100000          /* max size of matrix           */
#define   FLD_WTH       6               /* printing info                */

#define   MAX_MAT       100
#define   HEADER_LEN    4

#define   SPARSE        0.25

#define   TRI           1
#define   BI            2

#define CG_LIMIT       10000
#define CG_MR_LIMIT    250
#define INIT_CG_LIMIT  10

#define INIT_RES_LIMIT 1.0e-02 /*pow(2.0,-20.0)*/
#define RES_LIMIT      1.0e-14 /*pow(2.0,-20.0)*/
#define RES_ORTHOG     1.0e-05 /*pow(2.0,-20.0)*/

#define PI             acos(-1.0)


/* c defined function proto. Definition in respective .h file */
/* some conflicts amongst include files ... so comment out   */
/* Note: return values for fprintf and printf are ignored ... */
double sqrt(double x);
double fabs(double x);
/*
FILE *fopen(char *name,char *mode);
int fclose(FILE *fp);

int fprintf(FILE * stream,const char *format, ...);
int printf(const char *format, ...);
char *strtok(char *s,const char *ct);
char *fgets(char * s,int n,FILE *stream);
void rewind(FILE * stream);
void *malloc(size_t size);
int atoi(const char *s);
void free(void *p);
void  exit(int status);
*/

/* my protos */
int log_2_floor(int n);
void compress_adj_list(void);
matrix_id_ptr extract_sub_matrix(int *m,int n);
struct tree_node *new_tree_node();
int SMI_csc(int *ptrs, int *cols, double *vals, int lda, int type);
int SMI_csc_map(int *map, int *ptrs, int *cols, double *vals, int lda, 
		int type);


double *lanczos(double *rhs, int n, int *ia, int *ja, double *vals);
double lap_sub_matrix(int *list, int n, int **ia, int **ja, double **vals);
void ortho(double *r, int n);
double *Init_guess(int n);
double *Init_rhs(int n);
double hmt_ddot(double *x1, double *x2, int n);
void rev_daxpy(double *x1, double scale, double *x2, int n);
void hmt_daxpy(double *x1, double scale, double *x2, int n);
void matrix_mult(double *res, double *v, int n, int *ia, int *ja, double *aij);
void skew_matrix_mult(double *Ap_k, double *p_k, int n, int *ia, int *ja, 
		      double *vals, double lambda);
double *cg_min_res(double *rhs, int n, int *ia, int *ja, double *vals, 
		   double lambda);



void tri_calc_evec(double *d, double *sd, int n, double lambda, double *ret);
int trid_slv(double *d, double *sd, int n, double lambda, double *b);
double *trid_mult(double *d, double *sd, int n, double *b);



double median(double *fiedler, int n);

double *perturb(double *rhs, int n);
void check_row_sum(int n, int *ia, int *ja, double *vals);
int *separate_graph(int *list, int n, double *rhs);


/* Global Variables (visible)     */

/* Global Variables (not visible) */
extern matrix_id_ptr matrix_list[MAX_MAT];
extern matrix_id_ptr active_matrix;
extern int num_mat;

/* hack */
extern double res_limit;


extern int *pairs;
extern double *vals;
extern int sda;
extern double num_cuts;

extern double start_nnz;
extern double num_cuts_level;
extern double num_edges_level;
extern int sep_sz;

double wt_subg_level;
double wt_cuts_level;
double sep_sz_level;
int nnz_level;


int *color_proc, *color_sep;
queue_ADT done_q;
struct tree_node *root;


/********************************sparse_matrix.c*******************************
Function: Init_rhs_list()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
double *Init_rhs_list(double *rhs, int *list, int n)
{
  int i, j, root_n, numin, offset;
  double *tmp;
  double denom;


  tmp = rhs;
  for (i=j=0; i<n; i++)
    {j+= list[i];}
  j/=n;

  for (i=0; i<n; i++)
    {
      if (list[i]< j)
	{*tmp++ = -(list[i]+1.0);}
      else
	{*tmp++ =  (list[i]+1.0);}
    }

  /* srand(101001729); */
  /* *tmp++ = ((double) rand())/((double)RAND_MAX); */
  /*
  ortho(rhs,n);
  normalize(rhs,n);
  */
  return(rhs);
}



/********************************sparse_matrix.c*******************************
Function: lap_sub_matrix()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
double
lap_sub_matrix(int *list, register int n, int **ia, int **ja, double **vals)
{
  register int i, j, k, elm;
  register int *iptr_ia, *iptr_ja;
  register int *iptr_row, *iptr_list;
  double *dptr_vals, *v_hold, row_sum;
  int **rows;
  double **v;
  double *dptr, wt;
  int ct;


  error_msg_warning("lap_sub_matrix() :: begin\n");

  /* need len and indice information */
  rows = active_matrix->row_indices;
  v    = active_matrix->row_vals;

  /* get space for ia offsets and initialize to all zeros */
  *ia = iptr_ia = (int *) bss_malloc(INT_LEN*(n+2));
  ivec_zero(iptr_ia,n+2);

  /* get space for ja and set ia */
  *iptr_ia = 0;
  for (i=j=0; i<n; i++)
    {
      iptr_list = list;
      iptr_row = *(rows + *(iptr_list+i));

      k=0;
      ct=0;
      while ((elm=*iptr_row++)>=0)
	{
	  while ((*iptr_list>=0)&&(elm > *iptr_list))
	    {iptr_list++;}
	  if (*iptr_list<0) break;
	  if (elm == *iptr_list)
	    {k++; iptr_list++;}
	}

      if ((k<2) && (n>1))
	{error_msg_warning("lap_sub_matrix() :: disconnected graph?");}

      j+=k;
      *(iptr_ia + 1) = *(iptr_ia) + k;
      iptr_ia++;
    }
  *(iptr_ia + 1) = -1;

  /* ok ... get space for ja and vals */
  *ja = iptr_ja = (int *) bss_malloc(INT_LEN*j);
  ivec_set(iptr_ja,0,j);
  *vals = dptr = (double *) bss_malloc(REAL_LEN*j);

  /* finish laplacian map ==> (ia,ja,vals) */
  iptr_ia = *ia;
  /* rvec_neg_one(dptr,j); */
  rvec_zero(dptr,j);
  
  for (wt=0.0,i=0; i<n; i++)
    {
      iptr_list = list;
      k = *(iptr_list+i);
      iptr_row = *(rows + k);
      dptr_vals = *(v + k);

      j=0;
      row_sum = 0.0;
      while ((elm=*iptr_row++)>=0)
	{
	  while ((*iptr_list>=0)&&(elm>*iptr_list))
	    {j++; iptr_list++;}
	  if (*iptr_list<0) break;
	  if (elm == *iptr_list)
	    {
	      *iptr_ja++=j++;
	      if (elm==k)
		{v_hold = dptr;}
	      else
		{
		  *dptr = *dptr_vals;
		  row_sum += *dptr;
		}

	      /* old */
	      /* = (double) (*(iptr_ia+1) - *(iptr_ia) - 1); iptr_ia++;}*/
	      /* for spd */
	      /*
	      {*dptr = (double) (*(iptr_ia+1) - *(iptr_ia)); iptr_ia++;}
	      */

	      dptr++;
	      iptr_list++;
	    }
	  dptr_vals++;
	}
      *v_hold = -row_sum;
      wt+= -row_sum;
    }

  check_row_sum(n,*ia,*ja,*vals);

  error_msg_warning("lap_sub_matrix() :: end\n");

  return(wt);

  /* later :: chk to make sure graph is connected */

#ifdef DEBUG_0
  iptr_ia = ia+1;
  iptr_ja = ja;
  dptr   = vals;
  for (i=0; i<n; iptr_ia++, i++)
    {
      k = *(iptr_ia) - *(iptr_ia-1);
      for (j=0; j<k; j++)
	{printf("%3d ",*iptr_ja++);}
      printf("\n");

      for (j=0; j<k; j++)
	{printf("%3.0f ",*dptr++);}
      printf("\n");
    }
#endif
}



/********************************sparse_matrix.c*******************************
Function: 

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
double *
perturb(double *rhs, int n)
{
  int i;
  double med, *dptr1;

  return(rhs);

  /* slight perturbation? */  
  med=0.0;
  dptr1 = rhs;
  for (i=0; i<n; i++)
    {med += fabs(*dptr1++);}
  
  if (med<EPS)
    {
      printf("sum=%.16e, norm=%.16e\n",med,hmt_ddot(rhs,rhs,n));
      med=1.0;
      dptr1 = rhs;
      for (i=0; i<n; i++)
	{*dptr1++ += (med*rand())/((double)RAND_MAX);}
    }
  else
    {
      med/= 5.0*n;  
      dptr1 = rhs;
      srand(102031);
      for (i=0; i<n; i++)
	{*dptr1++ += (med*rand())/((double)RAND_MAX);}
    }
  return(rhs);
}



/********************************sparse_matrix.c*******************************
Function: 

Input : 
Output: 
Return: 
Description:  

Make sure that row sum of laplacian matrix is in fact zero as it should be.

*********************************sparse_matrix.c******************************/
void
check_row_sum(int n, int *ia, int *ja, double *vals)
{
  int i;
  double *b, *x;


  error_msg_warning("check_row_sum() :: begin\n");

  b = (double *) bss_malloc(n*REAL_LEN);
  x = (double *) bss_malloc(n*REAL_LEN);

  rvec_set(x,1.0,n);
  matrix_mult(b,x,n,ia,ja,vals);

  for (i=0; i<n; i++)
    {
      if (fabs(b[i])>1.0e-14) 
	{error_msg_fatal("check_row_sum() :: row sum %d = %g\n",i);}
    }


  bss_free(b);
  bss_free(x);

  error_msg_warning("check_row_sum() :: end\n");
}



/********************************sparse_matrix.c*******************************
Function: restricted_kl()

Input : 
Output: 
Return: 
Description:  

   well you know ... I'll get to it!
*********************************sparse_matrix.c******************************/
void
restricted_kl(int *list, int n, int *ia, int *ja, double *vals, int *color)
{
  ;
}



/********************************sparse_matrix.c*******************************
Function: brute_separate_graph()

Input : 
Output: 
Return: 
Description:  

bisect graph and highlight those which have elements with cut edges
   o <0 lhs
   o -2 lhs element has edge cut
   o >0 rhs
   o +2 lhs element has edge cut

instead of using lanczos why not use bandwidth reduction!

*********************************sparse_matrix.c******************************/
int *
brute_sep_graph(int *list, int n, int *ia, int *ja, double *vals, double *rhs, 
		int *ns)
{
  int i,j,k;
  int *color;
  double *x, *b, rs, nc=0.0;
  double *fiedler;
  int *map;
  int nnz, uh_nnz, bw;
  int *iptr, *sm_in, *sm_out;
  int v1,v2,v3;
  double d12,d13;

  error_msg_warning("brute_sep_graph() :: begin\n");
  color = (int *) bss_malloc(n*INT_LEN);

  if (n<1)
    {error_msg_fatal("brute_sep_graph() :: n=%d\n",n);}
  else if (n==1)
    {color[0] = -2; return(color);}
  else if (n==2)
    {color[0] = -2; color[1] =  2; return(color);}
  else
    {
      nnz = ia[n];
      uh_nnz = nnz - n;
      if (uh_nnz&1)
	{error_msg_fatal("brute_sep_graph_() :: matrix not symmetric!!\n");}
      uh_nnz>>=1;
      sm_in  = iptr = (int *) bss_malloc((2*(uh_nnz+1))*INT_LEN);
      sm_out = (int *) bss_malloc(n*INT_LEN);
      *iptr++ = n;
      *iptr++ = k = uh_nnz;
      for (i=j=0; i<n; i++)
	{
	  while (j<ia[i+1])
	    {
	      if (ja[j]>i)
		{*iptr++ = i+1; *iptr++ = ja[j]+1; k--;}
	      j++;
	    }
	}
      if (k)
	{error_msg_fatal("brute_sep_graph () :: pair list failure!\n");}

      bw=sm_bandwidth_reduction(sm_in,sm_out,-2);
      if (bw<=0)
	{error_msg_fatal("brute_sep_graph_ () :: sm_...() failure!\n");}

      v1 = sm_out[0]-1;
      v2 = sm_out[1]-1;
      v3 = sm_out[2]-1;

      for (i=0; i<ia[1]; i++)
	{ 
	  if (ja[i]==v3)
	    {break;}
	}

      /* found */
      if (i!=ia[1])
	{
	  if (rsb_dim==_2D)
	    {
	      d12  = pow(xc[list[v1]]-xc[list[v2]],2.0);
	      d12 += pow(yc[list[v1]]-yc[list[v2]],2.0);
	      d12 = sqrt(d12);
	      d13  = pow(xc[list[v1]]-xc[list[v3]],2.0);
	      d13 += pow(yc[list[v1]]-yc[list[v3]],2.0);
	      d13 = sqrt(d13);
	      if (d13<d12)
		{j=sm_out[1]; sm_out[1]=sm_out[2]; sm_out[2] = j;}
	    }
	  else
	    {
	      d12  = pow(xc[list[v1]]-xc[list[v2]],2.0);
	      d12 += pow(yc[list[v1]]-yc[list[v2]],2.0);
	      d12 += pow(zc[list[v1]]-zc[list[v2]],2.0);
	      d12 = sqrt(d12);
	      d13  = pow(xc[list[v1]]-xc[list[v3]],2.0);
	      d13 += pow(yc[list[v1]]-yc[list[v3]],2.0);
	      d13 += pow(zc[list[v1]]-zc[list[v3]],2.0);
	      d13 = sqrt(d13);
	      if (d13<d12)
		{j=sm_out[1]; sm_out[1]=sm_out[2]; sm_out[2] = j;}
	    }
	}

      for (i=0; i<n/2; i++)
	{color[sm_out[i]-1] = -2;}
      for (i=n/2; i<n; i++)
	{color[sm_out[i]-1] =  2;}
      bss_free(sm_in);
      bss_free(sm_out);
    }
      

  error_msg_warning("brute_sep_graph() :: end\n");
  return(color);
}


/********************************sparse_matrix.c*******************************
Function: separate_graph()

Input : 
Output: 
Return: 
Description:  

bisect graph and highlight those which have elements with cut edges
   o <0 lhs
   o -2 lhs element has edge cut
   o >0 rhs
   o +2 lhs element has edge cut


*********************************sparse_matrix.c******************************/
int *
sep_graph(int *list, int n, int *ia, int *ja, double *vals, double *rhs, int *ns)
{
  int i,j,k;
  int *color;
  double *x, *b, rs, nc=0.0;
  double *fiedler;
  int *map;


  error_msg_warning("sep_graph() :: begin\n");


  /* brute separate */
  if (n<=4)
    {return(brute_sep_graph(list,n,ia,ja,vals,rhs,ns));}

  color = (int *) bss_malloc(n*INT_LEN);
  map   = (int *) bss_malloc(n*INT_LEN);
  ivec_c_index(map,n);
  x     = (double *) bss_malloc(n*REAL_LEN);
  rvec_zero(x,n);
  b     = (double *) bss_malloc(n*REAL_LEN);
  rvec_zero(b,n);


  /* get fiedler vector associated w/subgraph */
  fiedler=lanczos(perturb(rhs,n), n, ia, ja, vals);

  /* determine lhs,rhs */
  rvec_sort_companion(fiedler,map,n);
  for (i=0; i<n/2; i++)
    {color[map[i]] = -1;}
  for (i=n/2; i<n; i++)
    {color[map[i]] =  1; x[map[i]]=1.0;}

  /* mark elms w/edges cut */
  matrix_mult(b,x,n,ia,ja,vals);
  for (i=0; i<n; i++)
    {
      rs = fabs(b[i]);
      if (rs>0.0) 
	{
	  color[i]*=2;
	  nc += rs;
	}
    }
  nc/=2;
  nc+=0.5;
  *ns = (int) nc;
  
  /* local clean up to minimize number of cut edges */
  restricted_kl(list,n,ia,ja,vals,color);

  bss_free(map);
  bss_free(x);
  bss_free(b);

  error_msg_warning("sep_graph() :: end\n");
  return(color);
}



/********************************sparse_matrix.c*******************************
Function: label_outflow()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int *
label_outflow(int *vertex, int nv, int *num_unique, int *num_outflow)
{
  int i, j, off, nvo, nvu;
  int *iptr;



  error_msg_warning("label_outflow() :: begin\n");


  *num_unique  = nvu = vertex[nv];
  *num_outflow = nvo = vertex[nv+1];

  error_msg_warning("nv=%d, nu=%d, no=%d\n",nv,nvu,nvo);

  iptr = (int *)bss_malloc(INT_LEN*(nvu+1));
  ivec_zero(iptr,nvu+1);

  for (i=j=0; i<nv; i++)
    {
      if (vertex[i] < 0)
	{
	  off = vertex[i] = -vertex[i];
	  if (!iptr[off])
	    {
	      iptr[off] = nvu-j;
	      j++;
	    }
	}
    }

  if (j!=nvo)
    {error_msg_fatal("label_outflow() :: found=%d, actual=%d\n",j,nvo);}

  error_msg_warning("label_outflow() :: end\n");

  return(iptr);
}



/********************************sparse_matrix.c*******************************
Function: det_sep()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int *
det_sep(int *list, int *color, int n, int *vertex, int nv, int *vertex_order, 
	int nvu, int ub, int *num_sep)
{
  int i, j, nc, off, ns, nsm, v;
  int *iptr;


  error_msg_warning("det_sep() :: begin\n");
  
  /* nsm should be an upper bound on the number of vertices to be labelled */
  nsm = n;
  nc = (rsb_dim==_2D) ? 4 : 8;
  nsm *= nc;
  nsm = MAX(nsm,nc);

  iptr = (int *) bss_malloc(INT_LEN*nsm);

  if (color==NULL)
    {
      for (ns=i=0; i<n; i++)
	{
	  off = nc*list[i];
	  for (j=0; j<nc; j++)
	    {
	      v = vertex[off+j];
	      if (vertex_order[v] == 0)
		{
#ifdef DEBUG		  
		  error_msg_warning("elm=%d,v=%d\n",off,v);
#endif
		  if (ns>=nsm)
		    {error_msg_fatal("oops ns about to exceed nsm=%d!\n",nsm);}
		  iptr[ns] = vertex_order[v] = ub-ns;
		  ns++;
		}
	    }
	}
    }
  else
    {
      for (i=0; i<n; i++)
	{
	  if (color[i] == -2)
	    {
	      off = nc*list[i];
	      for (j=0; j<nc; j++)
		{
		  v = vertex[off+j];
		  if (vertex_order[v] <= 0)
		    {
		      /* error_msg_warning("elm=%d,v=%d\n",off,v); */
		      vertex_order[v]--;
		    }
		}
	    }
	}


      for (ns=i=0; i<n; i++)
	{
	  if (color[i] == 2)
	    {
	      off = nc*list[i];
	      for (j=0; j<nc; j++)
		{
		  v = vertex[off+j];
		  if (vertex_order[v] < 0)
		    {
		      /* error_msg_warning("elm=%d,v=%d\n",off,v); */
		      if (ns>=nsm)
			{error_msg_fatal("oops ns about to exceed nsm=%d!\n",nsm);}

		      iptr[ns] = vertex_order[v] = ub-ns;
		      ns++;
		    }
		}
	    }
	}
    }

  /* reset the unlabelled vertices to 0 state */
  for (i=1; i<=nvu; i++)
    {if (vertex_order[i]<0) {vertex_order[i]=0;}}

  *num_sep = ns;
  ivec_sort(iptr,ns);

  error_msg_warning("det_sep() :: end\n");
  return(iptr);
}



/********************************sparse_matrix.c*******************************
Function: det_lc()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
struct tree_node *
det_lc(int *list, int *color, double *rhs, int n, int *vertex, int nv, 
       int *vertex_order, int nvu, int lb, int *num_verts)
{
  int i, j, k;
  int nl,nc,v,off;
  int *iptr;
  double *dptr;
  struct tree_node *t_ptr;
  

  error_msg_warning("det_lc() :: begin\n");

  t_ptr = new_tree_node();
  nc = (rsb_dim==_2D) ? 4 : 8;
  for (nl=i=0; i<n; i++)
    {
      if (color[i] < 0)
	{ 
	  nl++;
	  off = nc*list[i];
	  for (j=0; j<nc; j++)
	    {
	      v = vertex[off+j];
	      if (vertex_order[v] <= 0)
		{vertex_order[v]--;}
	    }
	}
    }


  t_ptr->num_act = nl;
  t_ptr->active = iptr = (int *) bss_malloc((nl+1)*sizeof(int));
  t_ptr->rhs    = dptr = (double *) bss_malloc(nl*sizeof(double));
  for (nl=i=0; i<n; i++)
    {
      if (color[i] < 0)
	{
	  nl++; 
	  *iptr++ = list[i];
	  *dptr++ = rhs[i];
	}
    }
  *iptr = -1;

  if (nl!=t_ptr->num_act)
    {error_msg_fatal("lc ct off!");}

  /* reset the unlabelled vertices to 0 state */
  for (nc=0,i=1; i<=nvu; i++)
    {if (vertex_order[i]<0) {nc++; vertex_order[i]=0;}}


  t_ptr->lb = lb;
  *num_verts = nc;
  t_ptr->ub  = nc+lb-1;


  error_msg_warning("det_lc() :: end\n");

  return(t_ptr);
}



/********************************sparse_matrix.c*******************************
Function: det_rc()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
struct tree_node *
det_rc(int *list, int *color, double *rhs, int n, int *vertex, int nv, 
       int *vertex_order, int nvu, int lb, int *num_verts)
{
  int i, j, k;
  int nl,nc,v,off;
  int *iptr;
  double *dptr;
  struct tree_node *t_ptr;
  

  error_msg_warning("det_rc() :: begin\n");

  t_ptr = new_tree_node();
  nc = (rsb_dim==_2D) ? 4 : 8;
  for (nl=i=0; i<n; i++)
    {
      if (color[i] > 0)
	{
	  nl++;
	  off = nc*list[i];
	  for (j=0; j<nc; j++)
	    {
	      v = vertex[off+j];
	      if (vertex_order[v] <= 0)
		{vertex_order[v]--;}
	    }
	}
    }

  t_ptr->num_act = nl;
  t_ptr->active = iptr = (int *) bss_malloc((nl+1)*sizeof(int));
  t_ptr->rhs    = dptr = (double *) bss_malloc(nl*sizeof(double));
  for (nl=i=0; i<n; i++)
    {
      if (color[i] > 0)
	{
	  nl++; 
	  *iptr++ = list[i];
	  *dptr++ = rhs[i];
	}
    }
  *iptr = -1;

  if (nl!=t_ptr->num_act)
    {error_msg_fatal("lc ct off!");}

  /* reset the unlabelled vertices to 0 state */
  for (nc=0,i=1; i<=nvu; i++)
    {if (vertex_order[i]<0) {nc++; vertex_order[i]=0;}}


  t_ptr->lb = lb;
  *num_verts = nc;
  t_ptr->ub  = nc+lb-1;

  error_msg_warning("det_rc() :: end\n");

  return(t_ptr);
}



/********************************sparse_matrix.c*******************************
Function: SMI_rsb()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
void
SMI_rsb(int mat, int *color, int level, int *in_vertex, int nv)
{
  matrix_id_ptr m;
  int i, depth, rank, lda, *iptr, *order;
  stack_ADT s;
  queue_ADT in_q, out_q, tmp_q;
  struct tree_node *t_ptr, *t_lc, *t_rc;
  int *tmp_color;
  int nl, nr, ns, j, k;
  double *dptr, *rhs;
  double *fill_level, dtmp;
  double wt_subg_old=0.0;
  double wt_cuts_old=0.0;
  int *ia, *ja, n;
  double *vals;
  int *list, *vertex_order, nvu, nvo;
  int *vertex;

  /* check valid matrix handle */
  if ((mat<1)||(mat>MAX_MAT)||((m=active_matrix=matrix_list[--mat]) == NULL))
    {error_msg_fatal("SMI_rsb() :: matrix handle is NULL!!!");}

  /* problem size */
  lda = m->lda;
  start_nnz = nnz_adj_list(m->adj_list,lda);
#ifdef DEBUG  
  error_msg_warning("Start RSB on %d elements\n",lda);
#endif

  /* check that it's symmetric */
  if (!is_adj_list_symm(m->adj_list,lda))
    {error_msg_fatal("SMI_rsb() :: not symmetric!!!");}

  /* compress ll into vector representation */
  compress_adj_list();

  /* queues/stack for building/processing tree */
  done_q = new_queue();  
  in_q   = new_queue();  
  out_q  = new_queue();  
  s = new_stack();  

  /* set tree root to original problem */
  /* later :: use -1 flag to indicate all lda rows */
  root = t_ptr = new_tree_node();
  t_ptr->num_act = lda;
  t_ptr->active = iptr = (int *) bss_malloc(INT_LEN*(lda+1));
  ivec_c_index(iptr,lda);
  *(iptr+lda) = -1;
  t_ptr->rhs = dptr = (double *) bss_malloc(lda*REAL_LEN);
  Init_rhs_list(dptr,iptr,lda);
  t_ptr->depth = 0;

  /* zero out vector and pre-label outflow bc vertices */
  vertex = (int *) bss_malloc(INT_LEN*(nv+2));
  ivec_copy(vertex,in_vertex,nv+2);
  vertex_order = label_outflow(vertex,nv,&nvu,&nvo);
  t_ptr->lb = 1;
  t_ptr->ub = nvu-nvo;

#ifdef DEBUG
  for (k=0; k<=nvu; k++)
    {error_msg_warning("%d %d\n",k,vertex_order[k]);}
#endif

  /* queue root */
  enqueue(out_q,t_ptr);

  /* init color ==> elm to proc map) */
  ivec_set(color,-1,lda); 

  /* init order (row->xxt ordering) */
  /*
  order = (int *) bss_malloc(lda*INT_LEN);
  ivec_zero(order,lda);
  */

  /* create binary tree and process all non-leaf nodes */
  num_cuts=0.0;
  depth = MIN(log_2_floor(lda),level);
  /*  fill_level = (double *)bss_malloc(depth*REAL_LEN); */
  for(i=0; i<depth; i++)
    {
      for (rank=len_queue(out_q),j=0; j<rank; j++)
	{
	  /* subproblem parameters */
	  t_ptr = dequeue(out_q);
	  t_ptr->depth = i;
	  t_ptr->rank  = j;
	  n     = t_ptr->num_act;
	  list  = t_ptr->active;
	  rhs   = t_ptr->rhs;

	  error_msg_warning("i=%d, j=%d, lb=%d, ub=%d\n",i,j,t_ptr->lb,t_ptr->ub);

#ifdef DEBUG
	  error_msg_warning("RSB :: (%d,%d) on %d elements\n",i,rank,n);
	  if (!rhs||!list) {error_msg_fatal("SMI_rsb() :: null rhs or list?\n");}
#endif

	  /* extract corresponding submatrix */
	  wt_subg_level += lap_sub_matrix(list, n, &ia, &ja, &vals);

	  /* rsb on subgraph */
	  tmp_color = sep_graph(list, n, ia, ja, vals, rhs,&ns);

#ifdef DEBUG
	  for (k=0; k<=nvu; k++)
	    {error_msg_warning("%d %d\n",k,vertex_order[k]);}
#endif
	  /* determine separator vertices */
	  k=0;
	  nr = t_ptr->ub;
	  t_ptr->sep = det_sep(list,tmp_color,n,vertex,nv,vertex_order,nvu,nr,&ns);
	  t_ptr->num_sep   = ns;

#ifdef DEBUG
	  k+=ns;
	  error_msg_warning("num_sep=%d\n",ns);

	  for (k=0; k<=nvu; k++)
	    {error_msg_warning("%d %d\n",k,vertex_order[k]);}
#endif	      
	  /* set up for rsb on left subgraph */
	  nr = t_ptr->lb;
	  t_ptr->lc = det_lc(list,tmp_color,rhs,n,vertex,nv,vertex_order,nvu,
			     nr,&nl);
	  enqueue(in_q,t_ptr->lc);

#ifdef DEBUG
	  k+=nl;
	  error_msg_warning("lb=%d, num_left=%d\n",t_ptr->lb,nl);
#endif

	  /* set up for rsb on right subgraph */
	  t_ptr->rc = det_rc(list,tmp_color,rhs,n,vertex,nv,vertex_order,nvu,
			     t_ptr->lb+nl,&nr);
	  enqueue(in_q,t_ptr->rc);

#ifdef DEBUG
	  k+=nr;
	  error_msg_warning("num_right=%d\n",nr);
#endif

	  /* hold the entire tree as we build it */
	  enqueue(done_q,t_ptr);

	  /* might want to keep these for later use */
	  bss_free(tmp_color);
	  bss_free(ia);
	  bss_free(ja);
	  bss_free(vals);
	}

      tmp_q = out_q;
      out_q = in_q;
      in_q  = tmp_q;
    }

  /* i=depth :: label leaves (stored in in_q) and queue on done_q */
  for (rank=len_queue(out_q),j=0; j<rank; j++)
    {
      t_ptr = dequeue(out_q);
      t_ptr->depth = i;
      t_ptr->rank  = j;
      n     = t_ptr->num_act;
      list  = t_ptr->active;

      error_msg_warning("i=%d, j=%d, lb=%d, ub=%d, n=%d, who=%d\n",i,j,
	     t_ptr->lb,t_ptr->ub,n,list[0]);
      /* process */
      nl = t_ptr->lb;
      nr = t_ptr->ub;
      if (nr<nl)
	{
	  t_ptr->sep = NULL;
	  t_ptr->num_sep   = 0;
	}
      else
	{
	  ns = 8;
	  t_ptr->sep = det_sep(list,NULL,n,vertex,nv,vertex_order,nvu,nr,&ns);
	  t_ptr->num_sep   = ns;
	}

      enqueue(done_q,t_ptr);
    }

  for (k=0; k<=nvu; k++)
    {error_msg_warning("%d %d\n",k,vertex_order[k]);}
      
  /* dump out needed xxt info */

  /* free tree */

#ifdef INFO
  printf("done_q :: %3d\n",len_queue(done_q));
  printf("in_q   :: %3d\n",len_queue(in_q));
  printf("out_q  :: %3d\n ",len_queue(out_q));
#endif

  if (ns!=nvu)
    {error_msg_warning("ns=%d, nvu=%d :: missing some vertices!\n",ns,nvu);}

  /* bss_free(fill_level); */
  /* bss_free(order); */

  if (len_queue(in_q))
    {error_msg_fatal("SMI_rsb() :: in_q not empty!\n");}
  free_queue(in_q);
 
 if (len_queue(out_q))
    {error_msg_fatal("SMI_rsb() :: out_q not empty!\n");}
  free_queue(out_q);

  free_stack(s);
}



/********************************sparse_matrix.c*******************************
Function: SMI_()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
SMI_extract_proc_map(int *color, int level)
{
  int i,j,k;
  int nl, nq,na,off;
  int *list;
  queue_ADT tmp_q, hold_q;
  struct tree_node *t_ptr;


  hold_q = new_queue();

  nl = 1;
  i = level;
  while (i--)
    {nl<<=1;}


  if (level==0)
    {off=0;}
  else
    {off = (1<<level)-1;}

  nq = len_queue(done_q);

  error_msg_warning("level=%d, nl=%d, off=%d, nq=%d\n",level,nl,off,nq);
  if ((off+nl) > nq)
    {error_msg_fatal("SMI_extract_proc() :: not that many in queue!\n");}
  
  for (i=0; i<off; i++)
    {
      t_ptr = (struct tree_node *) dequeue(done_q);
      enqueue(hold_q,t_ptr);
    }

  nl += off;
  for (k=2,i=off; i<nl; i++, k++)
    {
      t_ptr = (struct tree_node *) dequeue(done_q);
      enqueue(hold_q,t_ptr);
      
      na = t_ptr->num_act;
      list = t_ptr->active;
      printf("proc%2d has %2d elms :: ",k-2,na);
      for (j=0; j<na; j++)
	{
	  printf("%d ",list[j]);
	  color[list[j]] = k;
	}
      printf("\n");
    }
  printf("\n");


  for (i=nl; i<nq; i++)
    {
      t_ptr = (struct tree_node *) dequeue(done_q);
      enqueue(hold_q,t_ptr);
    }

  free_queue(done_q);
  done_q = hold_q;
}



/********************************sparse_matrix.c*******************************
Function: SMI_()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
SMI_rsb_clean()
{
  int i,j,k;
  int ns, depth, lda;
  struct tree_node *t_ptr;


  k=0;
  j=0;
  ns=0;

  root=NULL;
  while(len_queue(done_q))
    {
      t_ptr = dequeue(done_q);
      i = t_ptr->depth;
      /*
      if (i!=j)
	{printf("\n"); j++;}
      if (i!=depth)
	{printf("%d(%d)\t",t_ptr->num_sep,lda>>i); k+=(t_ptr->num_sep*(lda>>i));}
      */
      if (t_ptr->num_act)  {bss_free(t_ptr->active);}
      if (t_ptr->sep)      {bss_free(t_ptr->sep);}
      if (t_ptr->rhs)      {bss_free(t_ptr->rhs);}
      ns+=t_ptr->num_sep;
      bss_free(t_ptr);
    }

  /* free stacks, queues, ... */
  free_queue(done_q);

  done_q = NULL;
}


