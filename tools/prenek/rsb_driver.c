#include <stdio.h>
#include <math.h>


#include "const.h"
#include "types.h"
#include "queue.h"
#include "stack.h"
#include "bss_malloc.h"
#include "adj_list.h"
#include "ivec.h"
#include "rsb_driver.h"


int rsb_dim=-1, nnz=-1, uh_nnz=-1, matrix=-1, max_level = _MHC;
int max_repeats;
int bw=-1;
int seed=0;
REAL res_limit;
int  n, nv;
int  *ia,*ja,*pre_nek_color,*vertex;
static REAL *vals;
REAL *xc, *yc, *zc;
int max_proc;
int p,l;


static int smi_init_flag = FALSE;
    
static int ps, bs;


/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
Dump xxt elm/proc map as well as xxt firing order info
**********************************rsb_driver.c********************************/
#if UPCASE
void 
DUMP_XXT  (char *file_prefix, int *length)
#else
void 
dump_xxt_ (char *file_prefix, int *length)
#endif
{
  SMI_dump_xxt(file_prefix, *length);
}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
void 
FLUSH_IO  (void)
#else
void 
flush_io_ (void)
#endif
{
  fflush(stdout);
}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
void 
END_RSB  (void)
#else
void 
end_rsb_ (void)
#endif
{
  ia = ja = pre_nek_color = vertex = NULL;
  xc = yc = zc = vals = NULL;
  n = nv = nnz = uh_nnz = rsb_dim = bw = -1;
  seed = 0;
  max_level = _MHC;
  max_repeats = DEFAULT_MAX_REPEATS;
  max_proc  = 1<<_MHC;
  SMI_rsb_clean();

  printf("perm missing=%d, start w/%d missing\n",(perm_calls()-perm_frees())-ps,ps);
  printf("bss  missing=%d, start w/%d missing\n",(bss_calls()-bss_frees())-bs,bs);

  bss_init();
  perm_init();

  fflush(stdout);
}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
void 
CLEAR_RSB  (void)
#else
void 
clear_rsb_ (void)
#endif
{
  SMI_rsb_clean();

  printf("perm missing=%d, start w/%d missing\n",(perm_calls()-perm_frees())-ps,ps);
  printf("bss  missing=%d, start w/%d missing\n",(bss_calls()-bss_frees())-bs,bs);

  fflush(stdout);
}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
C_EXTRACT_PROC_MAP  (int *color)
#else
void
c_extract_proc_map_ (int *color)
#endif
{
  int n, level;
  float prs_read;


#if UPCASE
  PRS   ("Num Procs? :  $");
  KEYPAD (&prs_read);
#else
  prs_  ("Num Procs? : $");
  keypad_ (&prs_read);
#endif    


  n = (int) (prs_read+EPS);
  n = MIN(n,max_proc);
  n = MAX(n,1);
  level = log_2_floor(n);
  

  SMI_extract_proc_map(color,level);

#ifdef DEBUG
  error_msg_warning("extract_proc_map_() :: n=%d, Input=%g\n",n,prs_read);
#endif
}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
C_SET_LEVEL (void)
#else
void
c_set_level_(void)
#endif
{
  float prs_read;


#if UPCASE
  PRS  ("level? :  $");
  KEYPAD (&prs_read);
#else
  prs_  ("level? : $");
  keypad_ (&prs_read);
#endif    

  max_level = MIN(_MHC,(int) (prs_read+EPS));
  max_level = MAX(0,max_level);
  max_proc  = 1<<max_level;

#ifdef DEBUG
  error_msg_warning("c_set_level_() :: max_level=%d, Input=%g\n",max_level,
		    prs_read);
#endif
}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
void
C_SET_SEED (void)
#else
void
c_set_seed_(void)
#endif
{
  float prs_read;


#if UPCASE
  PRS  ("BWR Order? 2/1/0: $");
  KEYPAD (&prs_read);
#else
  prs_  ("BWR Order? 2/1/0: $");
  keypad_ (&prs_read);
#endif    

  seed = (int) (prs_read+EPS);

#ifdef DEBUG
  error_msg_warning("c_set_seed_() :: seed=%d, Input=%g\n",seed,prs_read);
#endif
}


/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
void
C_SET_MAX_LANCZOS (void)
#else
void
c_set_max_lanczos_(void)
#endif
{
  float prs_read;


#if UPCASE
  PRS  ("max ? $");
  KEYPAD (&prs_read);
#else
  prs_ ("max repeats? $");
  keypad_ (&prs_read);
#endif    

  max_repeats = (int) (prs_read+EPS);

#ifdef DEBUG
  error_msg_warning("c_set_max_repeats_() :: rep=%d, Input=%g\n",max_repeats,
		    prs_read);
#endif
}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
void
C_SET_MAX_REPEATS (void)
#else
void
c_set_max_repeats_(void)
#endif
{
  float prs_read;


#if UPCASE
  PRS  ("max repeats? $");
  KEYPAD (&prs_read);
#else
  prs_ ("max repeats? $");
  keypad_ (&prs_read);
#endif    

  max_repeats = (int) (prs_read+EPS);

#ifdef DEBUG
  error_msg_warning("c_set_max_repeats_() :: rep=%d, Input=%g\n",max_repeats,
		    prs_read);
#endif
}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
void
C_SET_TOLERANCE (void)
#else
void
c_set_tolerance_(void)
#endif
{
  float prs_read;


#if UPCASE
  PRS  ("Enter Residual: $");
  KEYPAD (&prs_read);
#else
  prs_ ("Enter Residual: $");
  keypad_ (&prs_read);
#endif    

  if ((prs_read<0.0)||(prs_read>=1.0))
    {
      if (rsb_dim==_2D)
	{res_limit=_2D_DEFAULT_RES_LIMIT;}
      else
	{res_limit=_3D_DEFAULT_RES_LIMIT;}
    }
  else if (prs_read==0.0)
    {res_limit=1.0e-15;}
  else
    {res_limit=prs_read;}

#ifdef DEBUG
  error_msg_warning("c_set_tolerance_() :: res_limit=%g, Input=%g\n",
		    res_limit,prs_read);
#endif
}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
void
C_PRS_DEFAULTS (void)
#else     
void
c_prs_defaults_(void)
#endif     
{
  float prs_read;


  prs_read = (float) res_limit;
#if UPCASE
  PRSI ("NEL $ ", &n);
  PRSI ("NV  $ ", &nv);
  PRSI ("Max Level $ ", &max_level);
  PRSI ("Max Repeats $ ", &max_repeats);
  PRSR ("Lanczos tolerance = $ ", &prs_read);
#else
  prsi_ ("NEL $ ", &n);
  prsi_ ("NV  $ ", &nv);
  prsi_ ("Max Level $ ", &max_level);
  prsi_ ("Max Repeats $ ", &max_repeats);
  prsr_ ("Lanczos tolerance = $ ", &prs_read);
#endif    

/*
  prsi_ ("Seed function by smbwr $ ", &seed);
  PRSI ("Seed function by smbwr $ ", &seed);
*/

}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
void 
BEGIN_RSB  (int *in_ia, int *in_ja, REAL *in_vals, REAL *in_xc, REAL *in_yc,
	    REAL *in_zc, int *in_color, int *in_n, int *in_vertex, 
	    int *in_nv, int *dim)
#else
void 
begin_rsb_ (int *in_ia, int *in_ja, REAL *in_vals, REAL *in_xc, REAL *in_yc,
	    REAL *in_zc, int *in_color, int *in_n, int *in_vertex, 
	    int *in_nv, int *dim)
#endif
{
  /* initialize malloc and sparse matrix packages */
  bss_init();
  perm_init();
  if (!smi_init_flag++)
    {SMI_init();}
  
  ps = perm_calls() - perm_frees();
  bs = bss_calls()  - bss_frees();

  if (ps||bs)
    {error_msg_fatal("begin_rsb() :: malloc stats bad %d,%d\n",ps,bs);}

  /* check csr input */
  if (!in_ia||!in_ja||!in_vals||!in_n)
    {error_msg_fatal("begin_rsb_() :: bad csr input!\n");}

  if (!in_xc||!in_yc||!in_zc)
    {error_msg_fatal("begin_rsb_() :: bad centroid input!\n");}

  if (!in_color||!in_vertex||!in_nv||!dim)
    {error_msg_fatal("begin_rsb_() :: bad misc input!\n");}
  
  /* hold elemental adjacency matrix and coordinates */
  n      = *in_n;
  ia     = in_ia;
  ja     = in_ja;
  vals   = in_vals;
  pre_nek_color  = in_color;
  xc     = in_xc;
  yc     = in_yc;
  zc     = in_zc;

  /* hold vertex numbering info */
  nv     = *in_nv;
  vertex = in_vertex;
  
  /* quick pseudo-symm check */
  nnz = ia[n]-1;
  uh_nnz = nnz - n;
  if (uh_nnz&1)
    {error_msg_fatal("begin_rsb_() :: matrix not symmetric!!\n");}
  uh_nnz>>=1;
  
  /* 0==>2D, 1==>3D */
  rsb_dim = *dim;
  if ((rsb_dim!=_2D)&&(rsb_dim!=_3D))
    {error_msg_fatal("begin_rsb_() :: dim must be _2D or _3D :: dim=%d\n",*dim);}

  if (rsb_dim==_2D)
    {res_limit=_2D_DEFAULT_RES_LIMIT;}
  else
    {res_limit=_3D_DEFAULT_RES_LIMIT;}  

  max_repeats = DEFAULT_MAX_REPEATS;

  bw = -1;
  seed = 0;
  max_level = _MHC;
  max_proc  = 1<<_MHC;
}



/*********************************rsb_driver.c*********************************
Function: 

Input : 
Output: 
Return: 
Description:  
**********************************rsb_driver.c********************************/
#if UPCASE
int
C_RSB  (void)
#else
int
c_rsb_ (void)
#endif
{
  int i, j, k, level;
  int *iptr, *sm_in, *sm_out;
  int mp;
  int success;


  /* determine max number of processors */
  level = log_2_floor(n);
     
  /* hack to walk through rsb */
  level = MIN(level,max_level);
  error_msg_warning("c_rsb_() :: max_proc=%d, level=%d\n",mp,level);

    
  /* use smbrc() to reorder. note fortran indexing 1...n -> 0...(n-1) */
  /* 0==>no reordering, 1==>ppv, 2==>try 'em all, default==>no reordering */

  /* create upper 1/2 pair list w/header for smbrc() */
  if ((seed==1)||(seed==2))
    {
      error_msg_fatal("c_rsb_() :: no sm_...() allowed!\n");
      sm_in  = iptr = (int *) bss_malloc((2*(uh_nnz+1))*INT_LEN);
      sm_out = (int *) bss_malloc(n*INT_LEN);
      *iptr++ = n;
      *iptr++ = k = uh_nnz;
      for (i=j=1; i<=n; i++)
	{
	  while (j<ia[i])
	    {
	      if (ja[j-1]>i)
		{*iptr++ = i; *iptr++ = ja[j-1]; k--;}
	      j++;
	    }
	}
      if (k)
	{error_msg_fatal("rsb_ () :: error creating smbwr() pair list\n");}
      bw=sm_bandwidth_reduction(sm_in,sm_out,-seed);
      matrix = SMI_csc_map(sm_out,ia,ja,vals,n,FULL);
#ifdef DEBUG      
      error_msg_warning("c_rsb_() :: bw=%d :: ",bw);
      for (i=0; i<n; i++)
	{printf("%d ",sm_out[i]);}
      printf("\n");
#endif
    }
  else
    {
#ifdef DEBUG      
      error_msg_warning("c_rsb_() :: using given ordering!\n");
#endif
      bw = (int) sqrt((REAL)n); 
      matrix = SMI_csc(ia,ja,vals,n,FULL);
    }

#ifdef DEBUG
  i = perm_calls() - perm_frees();
  j = bss_calls()  - bss_frees();
#endif
  success = SMI_rsb(matrix,pre_nek_color,level,vertex,nv);
#ifdef DEBUG
  printf("perm missing=%d, start w/%d missing\n",(perm_calls()-perm_frees())-i,i);
  printf("bss  missing=%d, start w/%d missing\n",(bss_calls()-bss_frees())-j,j);
#endif

  /* remap color if necessary */
  if ((seed==1)||(seed==2))
    {
      ivec_zero(sm_in,n);
      ivec_copy(sm_in,pre_nek_color,n);
      for (i=0; i<n; i++)
	{pre_nek_color[sm_out[i]-1] = sm_in[i];}
      bss_free(sm_in);
      bss_free(sm_out);
    }

#ifdef DEBUG
  SMI_fprint_id_hb_no_vals(matrix,FULL,"adj_list.hb");
#endif

  SMI_free_matrix(matrix);
  bss_stats();  
  perm_stats();

  if (success==SUCCESS)
    {return(0);}

  return(1);
}
