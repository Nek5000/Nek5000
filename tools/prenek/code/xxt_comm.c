#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>


#ifdef MPISRC
#include <mpi.h>
#endif


/* mission critical */
#ifdef MYMALLOC
#define PERM_MALLOC_BUF  16384
#define BSS_MALLOC_BUF   65536
#endif


#define MAX_MSG_BUF      4096
/*8192  32768  4096  65536  131072  262144  520192*/
#define THRESH    	0.2
#define N_HALF		4700
#define PRIV_BUF_SZ     45

#define REAL		1
#define INTEGER	        0

#define BYTE		8
#define BIT_0		0x1
#define BIT_1		0x2
#define BIT_2		0x4
#define BIT_3		0x8
#define BIT_4		0x10
#define BIT_5		0x20
#define BIT_6		0x40
#define BIT_7		0x80
#define BIT_8		0x100
#define BIT_9		0x200
#define BIT_10		0x400
#define BIT_11		0x800
#define BIT_12          0x1000
#define BIT_13          0x2000
#define BIT_14          0x4000
#define BIT_15          0x8000
#define BIT_16          0x10000 
#define BIT_17          0x20000 
#define BIT_18          0x40000 
#define BIT_19          0x80000
#define BIT_20          0x100000 
#define BIT_21          0x200000 
#define BIT_22          0x400000
#define BIT_23          0x800000 
#define BIT_24          0x1000000
#define BIT_25          0x2000000 
#define BIT_26          0x4000000 
#define BIT_27          0x8000000 
#define BIT_28          0x10000000  
#define BIT_29          0x20000000
#define BIT_30          0x40000000 
#define BIT_31          0x80000000
#define MAX_INT		0x7fffffff
#define MIN_INT		0x80000000
#define MAX_U_INT       0xffffffff
#define MIN_U_INT	0x0


#define FORTRAN		1
#define TO_FORT		1
#define C		0
#define TO_C		0
           
#define NON_UNIFORM     0
#define GL_MAX          1
#define GL_MIN          2
#define GL_MULT         3
#define GL_ADD          4
#define GL_B_XOR        5
#define GL_B_OR         6
#define GL_B_AND        7
#define GL_L_XOR        8
#define GL_L_OR         9
#define GL_L_AND        10


#define INT_LIST        2
#define BIT_MASK        1
#define NONE            0
#define ERROR		-1
#define NO_GS_TEMP      -2
#define NOT_ENOUGH_MEM  -3
#define BAD_ARGS        -4
#define PASS		1
#define FAIL		0

#define TRUE		1
#define FALSE		0


#define MAX_VEC		1674
#define FORMAT		30
#define MAX_COL_LEN    	100
#define MAX_LINE	FORMAT*MAX_COL_LEN
#define DELIM		" "
#define LINE		12
#define C_LINE		80


/* sorting args sorts.c */
#define   R_OPT		15
#define   NR_OPT	6
#define   SORT_OPT	6     
#define   MAX_STACK	50

#define LOG2(x)		log(x)/log(2)
#define SWAP(a,b)       temp=(a); (a)=(b); (b)=temp;
#define MAX(x,y)        ((x)>(y)) ? (x) : (y)
#define MIN(x,y)        ((x)<(y)) ? (x) : (y)

#define PTR_LEN		sizeof(char *)

#define INT_PTR_LEN	sizeof(int *)
#define INT_LEN		sizeof(int)
#define FLOAT 		float
#define FLOAT_LEN	sizeof(float)
#define CHAR_LEN        sizeof(char)
#define EPS             1.0e-6


/* -Dr8 as command line arg to icc */
#ifdef  r8
#undef  FLOAT
#undef  FLOAT_LEN 
#undef  EPS
#define FLOAT double
#define FLOAT_LEN	sizeof(double)
#define EPS             1.0e-12
#endif


/* local gs strength */
#define FULL         2
#define PARTIAL      1
#define NONE         0



typedef void (*vfp)();


typedef struct gather_scatter_id {
  int id;
 
  int nel_min;
  int nel_max;
  int nel_sum;
  int negl;
  int gl_max;
  int gl_min;
  int repeats;
  int ordered;
  int positive;
  double *vals;

  /* bit mask info */
  int mask_sz;
  int ngh_buf_sz;
  int *ngh_buf;
  int *neighbors;
  int num_nghs;

  int num_loads;

  /* repeats == true -> local info */
  int nel;         /* number of unique elememts */
  int *elms;       /* of size nel */
  int nel_total;
  int *local_elms; /* of size nel_total */
  int *companion;  /* of size nel_total */

  /* local info */
  int num_local_total;
  int local_strength;
  int num_local;
  int *num_local_reduce;
  int **local_reduce;
  int num_local_gop;
  int *num_gop_local_reduce;
  int **gop_local_reduce;

  /* pairwise info */
  int level;
  int num_pairs;
  int max_pairs;
  int *pair_list;
  int *msg_sizes;
  int **node_list;
  int len_pw_list;
  int *pw_elm_list;
  double *pw_vals;

#ifdef NXSRC
  int *msg_ids_in;
  int *msg_ids_out;
#endif
#ifdef MPISRC
  MPI_Request *msg_ids_in;
  MPI_Request *msg_ids_out;
#endif
  double *out;
  double *in;
  int msg_total;

  /* tree - crystal accumulator info */
  int max_left_over;
  int *pre;
  int *in_num;
  int *out_num;
  int **in_list;
  int **out_list;

  /* current memory status */
  int gl_bss_min;
  int gl_perm_min;

} gs_id;


/* global program control variables */
/* parallel.c */
static int p_init = 0;
static int my_id;
static int num_nodes;
static int boot_msg_buf[PRIV_BUF_SZ*INT_LEN];
static int floor_num_nodes, modfl_num_nodes, i_log2_num_nodes;
static double f_log2_num_nodes;
static unsigned int edge_bits[INT_LEN*BYTE];
static unsigned int edge_node[INT_LEN*BYTE];
static unsigned int h2l_edge_bits[INT_LEN*BYTE];
static unsigned int h2l_edge_node[INT_LEN*BYTE];
static int edge_not_pow_2;

/* ivec.c */
static int *offset_stack[2*MAX_STACK];
static int size_stack[MAX_STACK];

/* bss_malloc.c ...  stats */
static int    perm_req       = 0;
static int    num_perm_req   = 0;
static int    num_perm_frees = 0;
#ifdef MYMALLOC
static double perm_buf[PERM_MALLOC_BUF/8];
static double *perm_top;
#endif

static int    bss_req        = 0;
static int    num_bss_req    = 0;
static int    num_bss_frees  = 0;
#ifdef MYMALLOC
static double bss_buf[BSS_MALLOC_BUF/8];
static double *bss_top;
#endif

/* gs.c */
static int num_gs_ids = 0;

/* hack */
static gs_id *gs_handle = NULL;

static long msg_id[10];
  

/* the fortran interface */
void c_init_ (void);
void staged_iadd_ (int    *vals,int *level);
void staged_radd_ (double *vals,int *level);
void staged_gs_   (double *vals, double *work, int *level, int *segs);
void flush_io_ (void);
int gs_init_ (int *elms, int *nel, int *level);
void gs_gop_ (int *gs, double *vals);
void gs_free_ (int *gs);
void hmt_concat_ (double *vals, double *work, int *size);


/* my c prototypes */
void sgl_iadd(int    *vals, int level);
void sgl_radd(double *vals, int level);
void ssgl_radd(double *vals, double *work, int level, int *segs);
void fssgl_radd(double *vals, double *work, int level, int *segs);
void c_init (void);
void giop(int *vals, int *work, int n, int *oprs);
gs_id *gs_init(int *elms, int nel, int level);
void gs_gop(gs_id *gs, double *vals, char *operation);
void gs_free(gs_id *gs);
void gs_print_template(gs_id* gs, int who);
void hmt_concat(double *vals, double *work, int size);


void ivec_copy(int *arg1, int *arg2, int n);
void ivec_zero(int *arg1, int n);
void ivec_pos_one(int *arg1, int n);
void ivec_neg_one(int *arg1, int n);
void ivec_set(int *arg1, int arg2, int n);
int ivec_cmp(int *arg1, int *arg2, int n);
int ivec_lb(int *work, int n);
int ivec_ub(int *work, int n);
int ivec_sum(int *arg1, int n);
int ivec_u_sum(unsigned *arg1, int n);
int ivec_prod(int *arg1, int n);
vfp ivec_fct_addr(int type);
void ivec_non_uniform(int *arg1, int *arg2, int n, int *arg3);
void ivec_max(int *arg1, int *arg2, int n);
void ivec_min(int *arg1, int *arg2, int n);
void ivec_mult(int *arg1, int *arg2, int n);
void ivec_add(int *arg1, int *arg2, int n);
void ivec_xor(int *arg1, int *arg2, int n);
void ivec_or(int *arg1, int *arg2, int len);
void ivec_and(int *arg1, int *arg2, int len);
void ivec_lxor(int *arg1, int *arg2, int n);
void ivec_lor(int *arg1, int *arg2, int len);
void ivec_land(int *arg1, int *arg2, int len);
void ivec_or3 (int *arg1, int *arg2, int *arg3, int len);
void ivec_and3(int *arg1, int *arg2, int *arg3, int n);
int ivec_split_buf(int *buf1, int **buf2, int size);
void ivec_sort_companion(int *ar, int *ar2, int size);
void ivec_sort(int *ar, int size);
void SMI_sort(void *ar1, void *ar2, int size, int type);
int ivec_binary_search(int item, int *list, int n);
int ivec_linear_search(int item, int *list, int n);


void rvec_add(double *arg1, double *arg2, int n);
void rvec_copy(double *arg1, double *arg2, int n);
void rvec_zero(double *arg1, int n);

void *malloc(size_t size);
void  perm_init(void);
void *perm_malloc(size_t size);
int   perm_rem(void);
void  perm_free(void *ptr);
void  perm_stats(void);

void  bss_init(void);
void *bss_malloc(size_t size);
int   bss_rem(void);
void  bss_free(void *ptr);
void  bss_stats(void);



void error_msg_fatal(char *msg);
void error_msg_warning(char *msg);
void error_msg_warning_n0(char *msg);
void ierror_msg_warning(char *msg, int val);
void ierror_msg_warning_n0(char *msg, int val);

int div_ceil(int numin, int denom);
void set_bit_mask(int *bm, int len, int val);
int len_bit_mask(int num_items);
int ct_bits(char *ptr, int n);
void bm_to_proc(char *ptr, int p_mask, int *msg_list);
int len_buf(int item_size, int num_items);


/* PUBLIC - but not explicitly exported */
int gs_dump_ngh(gs_id *id, int loc_num, int *num, int *ngh_list);

/* PRIVATE */
static gs_id *gsi_check_args(int *elms, int nel, int level);
static void gsi_via_bit_mask(gs_id *gs);
static void gsi_via_int_list(gs_id *gs);
static void get_ngh_buf(gs_id *gs);
static void set_pairwise(gs_id *gs);
static void gs_gop_pairwise(gs_id *gs, double *in_vals, char *operation);
static void set_tree(gs_id *gs);
static void gs_gop_tree(gs_id *gs, double *vals, char *operation);
static void gs_gop_local(gs_id *gs, double *vals, char *operation);
static void gs_gop_local_in(gs_id *gs, double *vals, char *operation);
static void gs_gop_local_out(gs_id *gs, double *vals, char *operation);



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
c_init_ (void)
{
  c_init();
}



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
c_init (void)
{
  unsigned int edge;

/*  if (p_init++) return; */
  p_init++;

  if (sizeof(int) != sizeof(int *))
    {error_msg_fatal("Int/Int * Incompatible ... sorter hack!!!");}

  if ((sizeof(int) != sizeof(long))||(sizeof(int) != 4))
    {error_msg_fatal("Int/Long Incompatible or Int != 4 Bytes!");}

  if (sizeof(FLOAT)!=8)
    {error_msg_fatal("FLOAT/DOUBLE Incompatible?");}

  if (PRIV_BUF_SZ < 30)
    {error_msg_fatal("Check PRIV_BUF_SZ in const.h!");}
  
  my_id = mynode();
  num_nodes = numnodes();

  if (!num_nodes)
    {error_msg_fatal("Can't have no nodes!?!");}

  if (num_nodes> (MAX_INT >> 1))
  {error_msg_fatal("Can't have more then MAX_INT/2 nodes!!!");}

  bss_init();
  perm_init();

#ifdef INFO
  if (!my_id)
    {
      printf("BIT_31 = %u\n", BIT_31);
      printf("BIT_31 = %d\n", BIT_31);
      printf("MAX_INT  = %d\n", MAX_INT);
      printf("MIN_INT  = %d\n", MIN_INT);
    } 
#endif

  ivec_zero((int *)edge_bits,INT_LEN*BYTE);
  ivec_zero((int *)edge_node,INT_LEN*BYTE);
  ivec_zero((int *)h2l_edge_node,INT_LEN*BYTE);

  floor_num_nodes = 1;
  i_log2_num_nodes = modfl_num_nodes = 0;
  while (floor_num_nodes <= num_nodes)
    {
      edge_node[i_log2_num_nodes] = my_id ^ floor_num_nodes;
      edge_bits[i_log2_num_nodes] = floor_num_nodes;
      floor_num_nodes <<= 1; 
      i_log2_num_nodes++;
    }

  i_log2_num_nodes--;  
  floor_num_nodes >>= 1;
  modfl_num_nodes = (num_nodes - floor_num_nodes);

  /* low to high bits for contention free fan in/out */
  for (edge=0;edge<i_log2_num_nodes;edge++)
    {
      h2l_edge_node[edge] = edge_node[i_log2_num_nodes-edge-1];
    }
  h2l_edge_node[i_log2_num_nodes] = edge_node[i_log2_num_nodes];

  if ((my_id > 0) && (my_id <= modfl_num_nodes))
    {edge_node[i_log2_num_nodes]=edge_not_pow_2=((my_id|floor_num_nodes)-1);}
  else if (my_id >= floor_num_nodes)
    {
      ivec_zero((int *)edge_bits,INT_LEN*BYTE);
      ivec_zero((int *)edge_node,INT_LEN*BYTE);
      edge_node[i_log2_num_nodes]=edge_not_pow_2=((my_id^floor_num_nodes)+1);
    }
  else
    {edge_node[i_log2_num_nodes] = edge_not_pow_2 = 0;}

  f_log2_num_nodes = LOG2(num_nodes);
  if (((double)i_log2_num_nodes-f_log2_num_nodes) > EPS)
    {
      if (!my_id)
	{
      printf("i=%25.16f, f=%25.16f\n",1.0*i_log2_num_nodes,f_log2_num_nodes);
      printf("np1=%d, np2=%d\n",num_nodes,numnodes());
    }
      error_msg_warning_n0("Number of Processors is not a power of two!");
    }
  
  error_msg_warning_n0("did c_init()!!!");
}







/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
staged_iadd_ (int *gl_num, int *level)
{
  sgl_iadd(gl_num,*level);
}



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
sgl_iadd(int *vals, int level)
{
  int edge, type, dest, source, len, mask, ceil;
  long msg_id;
  int tmp, *work;


  /* all msgs will be of the same length */
  work = &tmp;
  len = INT_LEN;

  if (level > i_log2_num_nodes)
    {error_msg_fatal("sgl_add() :: level too big?");}

  if (level<=0)
    {return;}

  /* implement the mesh fan in/out exchange algorithm */
  if (my_id<floor_num_nodes)
    {
      mask = 0;
      for (edge = 0; edge < level; edge++)
	{
	  if (!(my_id & mask))
	    {
	      source = dest = edge_node[edge];
	      type = 10001 + my_id + (num_nodes*edge);
	      if (my_id > dest)
		{
		  msg_id = isend(type,vals,len,dest,0);
		  msgwait(msg_id);
		}
	      else
		{
		  type =  type - my_id + source;
		  msg_id = irecv(type,work,len);
		  msgwait(msg_id);
		  vals[0] += work[0];
		}
	    }
	  mask <<= 1;
	  mask += 1;
	}
    }

  if (my_id<floor_num_nodes)
    {
      mask >>= 1;
      for (edge = 0; edge < level; edge++)
	{
	  if (!(my_id & mask))
	    {
	      source = dest = edge_node[level-edge-1];
	      type = 10001 + my_id + (num_nodes*edge);
	      if (my_id < dest)
		{
		  msg_id = isend(type,vals,len,dest,0);
		  msgwait(msg_id);
		}
	      else
		{
		  type =  type - my_id + source;
		  msg_id = irecv(type,work,len);
		  msgwait(msg_id);
		  vals[0] = work[0];
		}
	    }
	  mask >>= 1;
	}
    }
}



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void 
flush_io_ (void)
{
  fflush(stdout);
}  



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
staged_radd_ (double *gl_num, int *level)
{
  sgl_radd(gl_num,*level);
}



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
sgl_radd(double *vals, int level)
{
  int edge, type, dest, source, len, mask, ceil;
  long msg_id;
  double tmp, *work;


  /* all msgs will be of the same length */
  work = &tmp;
  len = FLOAT_LEN;

  if (level > i_log2_num_nodes)
    {error_msg_fatal("sgl_add() :: level too big?");}

  if (level<=0)
    {return;}

  /* implement the mesh fan in/out exchange algorithm */
  if (my_id<floor_num_nodes)
    {
      mask = 0;
      for (edge = 0; edge < level; edge++)
	{
	  if (!(my_id & mask))
	    {
	      source = dest = edge_node[edge];
	      type = 1000001 + my_id + (num_nodes*edge);
	      if (my_id > dest)
		{
		  msg_id = isend(type,vals,len,dest,0);
		  msgwait(msg_id);
		}
	      else
		{
		  type =  type - my_id + source;
		  msg_id = irecv(type,work,len);
		  msgwait(msg_id);
		  vals[0] += work[0];
		}
	    }
	  mask <<= 1;
	  mask += 1;
	}
    }

  if (my_id<floor_num_nodes)
    {
      mask >>= 1;
      for (edge = 0; edge < level; edge++)
	{
	  if (!(my_id & mask))
	    {
	      source = dest = edge_node[level-edge-1];
	      type = 1000001 + my_id + (num_nodes*edge);
	      if (my_id < dest)
		{
		  msg_id = isend(type,vals,len,dest,0);
		  msgwait(msg_id);
		}
	      else
		{
		  type =  type - my_id + source;
		  msg_id = irecv(type,work,len);
		  msgwait(msg_id);
		  vals[0] = work[0];
		}
	    }
	  mask >>= 1;
	}
    }
}  




/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
 
ii+1 entries in seg :: 0 .. ii

******************************************************************************/
void 
staged_gs_ (register double *vals, register double *work, register int *level,
	    register int *segs)
{
  ssgl_radd(vals, work, *level, segs);
}



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
 
ii+1 entries in seg :: 0 .. ii

******************************************************************************/
void 
ssgl_radd(register double *vals, register double *work, register int level, 
	  register int *segs)
{
  register int edge, type, dest, source, mask;
  register int stage_n;

#ifdef DEBUG
  if (level > i_log2_num_nodes)
    {error_msg_fatal("sgl_add() :: level > log_2(P)!!!");}
#endif

  /* all msgs are *NOT* the same length */
  /* implement the mesh fan in/out exchange algorithm */
  for (mask=0, edge=0; edge<level; edge++, mask++)
    {
      stage_n = (segs[level] - segs[edge]);
      if (stage_n && !(my_id & mask))
	{
	  dest = edge_node[edge];
	  type = 100001 + my_id + (num_nodes*edge);
	  if (my_id>dest)
	    {csend(type, vals+segs[edge],stage_n<<3,dest,0);}
	  else
	    {
	      type =  type - my_id + dest;
	      crecv(type,work,stage_n<<3);
	      rvec_add(vals+segs[edge], work, stage_n); 
/*            daxpy(vals+segs[edge], work, stage_n); */
	    }
	}
      mask <<= 1;
    }
  mask>>=1;
  for (edge=0; edge<level; edge++)
    {
      stage_n = (segs[level] - segs[level-1-edge]);
      if (stage_n && !(my_id & mask))
	{
	  dest = edge_node[level-edge-1];
	  type = 10000001 + my_id + (num_nodes*edge);
	  if (my_id<dest)
	    {csend(type,vals+segs[level-1-edge],stage_n<<3,dest,0);}
	  else
	    {
	      type =  type - my_id + dest;
	      crecv(type,vals+segs[level-1-edge],stage_n<<3);
	    }
	}
      mask >>= 1;
    }
}  

int 
hmt_xor_ (register int *i1, register int *i2)
{
  return(*i1^*i2);
}




/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
 
ii+1 entries in seg :: 0 .. ii

******************************************************************************/
void 
old_ssgl_radd(double *vals, double *work, int level, int *segs)
{
  int edge, type, dest, source, len, mask, ceil;
  int offset, stage_n;

#ifdef DEBUG
  if (level > i_log2_num_nodes)
    {error_msg_fatal("sgl_add() :: level > log_2(P)!!!");}
#endif

/*
  ivec_neg_one(msg_id,10);
*/

  /* all msgs are *NOT* the same length */
  /* implement the mesh fan in/out exchange algorithm */
  mask = 0;
  for (edge = 0; edge < level; edge++)
    {
      stage_n = (segs[level] - segs[edge]);
      len = stage_n * FLOAT_LEN;

      if (stage_n && !(my_id & mask))
	{
	  source = dest = edge_node[edge];
	  type = 100001 + my_id + (num_nodes*edge);
	  if (my_id > dest)
	    {
	      csend(type, vals+segs[edge],len,dest,0);
	    }
	  else
	    {
	      type =  type - my_id + source;
	      crecv(type,work,len);
	      rvec_add(vals+segs[edge], work, stage_n);
	    }
	}
      mask <<= 1;
      mask += 1;
    }

  mask >>= 1;
  for (edge = 0; edge < level; edge++)
    {
      stage_n = (segs[level] - segs[level-1-edge]);
      len = stage_n * FLOAT_LEN;
      
      if (stage_n && !(my_id & mask))
	{
	  source = dest = edge_node[level-edge-1];
	  type = 10000001 + my_id + (num_nodes*edge);
	  if (my_id < dest)
	    {
	      csend(type,vals+segs[level-1-edge],len,dest,0);
	    }
	  else
	    {
	      type =  type - my_id + source;
	      crecv(type,vals+segs[level-1-edge],len);
	    }
	}
      mask >>= 1;
    }

/*
  for (edge=0;edge<10;edge++)
    {if (msg_id[edge] != -1) msgwait(msg_id[edge]);}
*/
}  



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
 
ii+1 entries in seg :: 0 .. ii

******************************************************************************/
void 
fast_staged_gs_ (register double *vals, register double *work, 
		 register int *level, register int *segs)
{
  fssgl_radd(vals, work, *level, segs);
}



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
 
ii+1 entries in seg :: 0 .. ii

******************************************************************************/
void 
fssgl_radd(register double *vals, register double *work, register int level, 
	   register int *segs)
{
  register int edge, dest, mask, stage_n, off, end;

#ifdef DEBUG
  if (level != i_log2_num_nodes)
    {error_msg_fatal("fssgl_radd() :: level != log_2(P)!!!");}

  for (edge = 0; edge < level; edge++)
    {
      if (!(segs[level] - segs[edge]))
	{ierror_msg_warning("fssgl_radd() :: stage_n=0 on",edge);}
    }
#endif

  /* all msgs are *NOT* the same length */
  /* fan in/out exchange ==> contention free on mesh */
  end = segs[level];
  for (mask=1, edge=0; edge<level; edge++, mask<<=1)
    {
      off = segs[edge];
      stage_n = (end - off);
      dest = my_id^mask;
      
      if (my_id > dest)
	{csend(my_id, vals+off,stage_n<<3,dest,0); break;}
      else
	{
	  crecv(dest,work,stage_n<<3);
	  if (stage_n<250)
	    {
	      for (dest=0; dest<stage_n; dest++)
		{vals[off+dest] += work[dest];}
	    }
	  else
	    {daxpy(stage_n,1.0,work,1,vals+off,1);}
	}
    }
 
  for (mask=num_nodes>>1, edge=0; edge<level; edge++, mask>>=1)
    {
      off = segs[level-1-edge];
      stage_n = (end - off);

      if (my_id%mask)
	{continue;}
      
      dest = my_id^mask;
      if (my_id < dest)
	{csend(my_id,vals+off,stage_n<<3,dest,0);}
      else
	{crecv(dest,vals+off,stage_n<<3);}
    }
}



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
 
ii+1 entries in seg :: 0 .. ii

******************************************************************************/
void 
calc_fssgl_radd(register double *vals, register double *work, register int level, 
	   register int *segs)
{
  register int edge, dest, mask, stage_n, off;

#ifdef DEBUG
  if (level != i_log2_num_nodes)
    {error_msg_fatal("fssgl_radd() :: level != log_2(P)!!!");}

  for (edge = 0; edge < level; edge++)
    {
      if (!(segs[level] - segs[edge]))
	{ierror_msg_warning("fssgl_radd() :: stage_n=0 on",edge);}
    }
#endif

  /* all msgs are *NOT* the same length */
  /* fan in/out exchange ==> contention free on mesh */
  for (mask=1, edge=0; edge<level; edge++, mask<<=1)
    {
      off = segs[edge];
      stage_n = (segs[level] - off);
      dest = my_id^mask;
      
      if (my_id > dest)
	{csend(my_id, vals+off,stage_n<<3,dest,0); break;}
      else
	{
	  crecv(dest,work,stage_n<<3);
	  for (dest=0; dest<stage_n; dest++)
	    {vals[off+dest] += work[dest];}
/*
	  if (stage_n>>5)
	    {daxpy(stage_n,1,work,1,vals+off,1);}
	  else
	    {
	    rvec_add(vals+off,work,stage_n); 
	    }
*/
	}
    }
  
  for (mask=num_nodes>>1, edge=0; edge<level; edge++, mask>>=1)
    {
      off = segs[level-1-edge];
      stage_n = (segs[level] - off);

      if (my_id%mask)
	{continue;}
      
      dest = my_id^mask;
      if (my_id < dest)
	{csend(my_id,vals+off,stage_n<<3,dest,0);}
      else
	{crecv(dest,vals+off,stage_n<<3);}
    }
}


void 
flush_msgs_ (void)
{
  register int edge;

  for (edge=0;edge<10;edge++)
    {if (msg_id[edge] != -1) msgwait(msg_id[edge]);}
}


/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
 
ii+1 entries in seg :: 0 .. ii
******************************************************************************/
void 
hmt_concat_ (double *vals, double *work, int *size)
{
  hmt_concat(vals, work, *size);
}



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
 
ii+1 entries in seg :: 0 .. ii

******************************************************************************/
void 
hmt_concat(double *vals, double *work, int size)
{
  int edge, type, dest, source, len, mask, ceil;
  long msg_id;
  int offset, stage_n;
  double *dptr;

  /* all msgs are *NOT* the same length */
  /* implement the mesh fan in/out exchange algorithm */
  mask = 0;
  stage_n = size;
  rvec_copy(work,vals,size);
  
  dptr  = work+size;
  for (edge = 0; edge < i_log2_num_nodes; edge++)
    {
      len = stage_n * FLOAT_LEN;

      if (!(my_id & mask))
	{
	  source = dest = edge_node[edge];
	  type = 100001 + my_id + (num_nodes*edge);
	  if (my_id > dest)
	    {
	      msg_id = isend(type, work, len,dest,0);
	      msgwait(msg_id);
	    }
	  else
	    {
	      type =  type - my_id + source;
	      msg_id = irecv(type, dptr,len);
	      msgwait(msg_id);
	    }
	}
      
#ifdef DEBUG_1      
      ierror_msg_warning_n0("stage_n = ",stage_n);
#endif

      dptr += stage_n;
      stage_n <<=1;
      mask <<= 1;
      mask += 1;
    }

  size = stage_n;
  stage_n >>=1;
  dptr -= stage_n;

  mask >>= 1;

  for (edge = 0; edge < i_log2_num_nodes; edge++)
    {
      len = (size-stage_n) * FLOAT_LEN;
      
      if (!(my_id & mask) && stage_n)
	{
	  source = dest = edge_node[i_log2_num_nodes-edge-1];
	  type = 10000001 + my_id + (num_nodes*edge);
	  if (my_id < dest)
	    {
	      msg_id = isend(type,dptr,len,dest,0);
	      msgwait(msg_id);
	    }
	  else
	    {
	      type =  type - my_id + source;
	      msg_id = irecv(type,dptr,len);
	      msgwait(msg_id);
	    }
	}

#ifdef DEBUG_1      
      ierror_msg_warning_n0("size-stage_n = ",size-stage_n);
#endif
      
      stage_n >>= 1;
      dptr -= stage_n;
      mask >>= 1;
    }
}




/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
gl_ivec_add(int *vals, int *work, int n)
{
  int edge, type, dest, source, len, mask, ceil;
  long msg_id;


  if (num_nodes<2)
    {return;}

  /* all msgs will be of the same length */
  len = n*INT_LEN;


  /* implement the mesh fan in/out exchange algorithm */
  if (my_id<floor_num_nodes)
    {
      mask = 0;
      for (edge = 0; edge < i_log2_num_nodes; edge++)
	{
	  if (!(my_id & mask))
	    {
	      source = dest = edge_node[edge];
	      type = 100001 + my_id + (num_nodes*edge);
	      if (my_id > dest)
		{
		  msg_id = isend(type,vals,len,dest,0);
		  msgwait(msg_id);
		}
	      else
		{
		  type =  type - my_id + source;
		  msg_id = irecv(type,work,len);
		  msgwait(msg_id);
		  ivec_add(vals, work, n);
		}
	    }
	  mask <<= 1;
	  mask += 1;
	}
    }

  if (my_id<floor_num_nodes)
    {
      mask >>= 1;
      for (edge = 0; edge < i_log2_num_nodes; edge++)
	{
	  if (!(my_id & mask))
	    {
	      source = dest = h2l_edge_node[edge];
	      type = 100001 + my_id + (num_nodes*edge);
	      if (my_id < dest)
		{
		  msg_id = isend(type,vals,len,dest,0);
		  msgwait(msg_id);
		}
	      else
		{
		  type =  type - my_id + source;
		  msg_id = irecv(type,work,len);
		  msgwait(msg_id);
		  ivec_copy(vals, work, n);
		}
	    }
	  mask >>= 1;
	}
    }
}  



/********************************ivec.c**************************************
Function ivec_copy()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_copy(register int *arg1, register int *arg2, register int n)
{
  while (n--)  {*arg1++ = *arg2++;} 
}



/********************************ivec.c**************************************
Function ivec_zero()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_zero(register int *arg1, register int n)
{
  while (n--)  {*arg1++ = 0;}
}



/********************************ivec.c**************************************
Function ivec_neg_one()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_neg_one(register int *arg1, register int n)
{
  while (n--)  {*arg1++ = -1;}
}



/********************************ivec.c**************************************
Function ivec_pos_one()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_pos_one(register int *arg1, register int n)
{
  while (n--)  {*arg1++ = 1;}
}



/********************************ivec.c**************************************
Function ivec_set()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_set(register int *arg1, register int arg2, register int n)
{
  while (n--)  {*arg1++ = arg2;}
}



/********************************ivec.c**************************************
Function ivec_cmp()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
int
ivec_cmp(register int *arg1, register int *arg2, register int n)
{
  while (n--)  {if (*arg1++ != *arg2++)  {return(FALSE);}}
  return(TRUE);
}



/********************************ivec.c**************************************
Function ivec_max()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_max(register int *arg1, register int *arg2, register int n)
{
  while (n--)  {*arg1 = MAX(*arg1,*arg2); arg1++; arg2++;}
}



/********************************ivec.c**************************************
Function ivec_min()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_min(register int *arg1, register int *arg2, register int n)
{
  while (n--)  {*(arg1) = MIN(*arg1,*arg2); arg1++; arg2++;}
}



/********************************ivec.c**************************************
Function ivec_mult()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_mult(register int *arg1, register int *arg2, register int n)
{
  while (n--)  {*arg1++ *= *arg2++;}
}



/********************************ivec.c**************************************
Function ivec_add()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_add(register int *arg1, register int *arg2, register int n)
{
  while (n--)  {*arg1++ += *arg2++;}
}



/********************************ivec.c**************************************
Function ivec_lxor()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_lxor(register int *arg1, register int *arg2, register int n)
{
  while (n--) {*arg1=((*arg1 || *arg2) && !(*arg1 && *arg2)) ; arg1++; arg2++;}
}



/********************************ivec.c**************************************
Function ivec_xor()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_xor(register int *arg1, register int *arg2, register int n)
{
  while (n--)  {*arg1++ ^= *arg2++;}
}



/********************************ivec.c**************************************
Function ivec_or()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_or(register int *arg1, register int *arg2, register int n)
{
  while (n--)  {*arg1++ |= *arg2++;}
}



/********************************ivec.c**************************************
Function ivec_lor()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_lor(register int *arg1, register int *arg2, register int n)
{
  while (n--)  {*arg1 = (*arg1 || *arg2); arg1++; arg2++;} 
}



/********************************ivec.c**************************************
Function ivec_or3()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_or3(register int *arg1, register int *arg2, register int *arg3, 
	 register int n)
{
  while (n--)  {*arg1++ = (*arg2++ | *arg3++);}
}



/********************************ivec.c**************************************
Function ivec_and()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_and(register int *arg1, register int *arg2, register int n)
{
  while (n--)  {*arg1++ &= *arg2++;}
}



/********************************ivec.c**************************************
Function ivec_land()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_land(register int *arg1, register int *arg2, register int n)
{
  while (n--) {*arg1++ = (*arg1 && *arg2); arg1++; arg2++;} 
}



/********************************ivec.c**************************************
Function ivec_and3()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_and3(register int *arg1, register int *arg2, register int *arg3, 
	  register int n)
{
  while (n--)  {*arg1++ = (*arg2++ & *arg3++);}
}



/********************************ivec.c**************************************
Function ivec_sum

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
int 
ivec_sum(register int *arg1, register int n)
{
  register int tmp = 0;


  while (n--) {tmp += *arg1++;}
  return(tmp);
}



/********************************ivec.c**************************************
Function ivec_prod

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
int 
ivec_prod(register int *arg1, register int n)
{
  register int tmp = 1;


  while (n--)  {tmp *= *arg1++;}
  return(tmp);
}



/********************************ivec.c**************************************
Function ivec_u_sum

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
int 
ivec_u_sum(register unsigned *arg1, register int n)
{
  register unsigned tmp = 0;


  while (n--)  {tmp += *arg1++;}
  return(tmp);
}



/********************************ivec.c**************************************
Function ivec_lb()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
int 
ivec_lb(register int *arg1, register int n)
{
  register int min = MAX_INT;


  while (n--)  {min = MIN(min,*arg1); arg1++;}
  return(min);
}



/********************************ivec.c**************************************
Function ivec_ub()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
int 
ivec_ub(register int *arg1, register int n)
{
  register int max = MIN_INT;


  while (n--)  {max = MAX(max,*arg1); arg1++;}
  return(max);
}



/********************************ivec.c**************************************
Function split_buf()

Input : 
Output: 
Return: 
Description: 

assumes that sizeof(int) == 4bytes!!!
*********************************ivec.c*************************************/
int
ivec_split_buf(int *buf1, int **buf2, register int size)
{
  *buf2 = (buf1 + (size>>3));
  return(size);
}



/********************************ivec.c**************************************
Function ivec_non_uniform()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
ivec_non_uniform(int *arg1, int *arg2, register int n, register int *arg3)
{
  register int i, j, type;


  /* LATER: if we're really motivated we can sort and then unsort */
  for (i=0;i<n;)
    {
      /* clump 'em for now */
      j=i+1;
      type = arg3[i];
      while ((j<n)&&(arg3[j]==type))
	{j++;}
      
      /* how many together */
      j -= i;

      /* call appropriate ivec function */
      if (type == GL_MAX)
	{ivec_max(arg1,arg2,j);}
      else if (type == GL_MIN)
	{ivec_min(arg1,arg2,j);}
      else if (type == GL_MULT)
	{ivec_mult(arg1,arg2,j);}
      else if (type == GL_ADD)
	{ivec_add(arg1,arg2,j);}
      else if (type == GL_B_XOR)
	{ivec_xor(arg1,arg2,j);}
      else if (type == GL_B_OR)
	{ivec_or(arg1,arg2,j);}
      else if (type == GL_B_AND)  
	{ivec_and(arg1,arg2,j);}
      else if (type == GL_L_XOR)
	{ivec_lxor(arg1,arg2,j);}
      else if (type == GL_L_OR)
	{ivec_lor(arg1,arg2,j);}
      else if (type == GL_L_AND)   
	{ivec_land(arg1,arg2,j);}
      else
	{error_msg_fatal("unrecognized type passed to ivec_non_uniform()!");}

      arg1+=j; arg2+=j; i+=j;
    }
}



/********************************ivec.c**************************************
Function ivec_addr()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
vfp ivec_fct_addr(register int type)
{
  if (type == NON_UNIFORM)
    {return(&ivec_non_uniform);}
  else if (type == GL_MAX)
    {return(&ivec_max);}
  else if (type == GL_MIN)
    {return(&ivec_min);}
  else if (type == GL_MULT)
    {return(&ivec_mult);}
  else if (type == GL_ADD)
    {return(&ivec_add);}
  else if (type == GL_B_XOR)
    {return(&ivec_xor);}
  else if (type == GL_B_OR)
    {return(&ivec_or);}
  else if (type == GL_B_AND)  
    {return(&ivec_and);}
  else if (type == GL_L_XOR)
    {return(&ivec_lxor);}
  else if (type == GL_L_OR)
    {return(&ivec_lor);}
  else if (type == GL_L_AND)   
    {return(&ivec_land);}
 
  /* catch all ... not good if we get here */
  return(NULL);
}


/********************************ivec.c**************************************
Function ct_bits()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
static
int 
ivec_ct_bits(register int *ptr, register int n)
{
  register int tmp=0;


  /* should expand to full 32 bit */
  while (n--)
    {
      if (*ptr&128) {tmp++;}
      if (*ptr&64)  {tmp++;}
      if (*ptr&32)  {tmp++;}
      if (*ptr&16)  {tmp++;}
      if (*ptr&8)   {tmp++;}
      if (*ptr&4)   {tmp++;}
      if (*ptr&2)   {tmp++;}
      if (*ptr&1)   {tmp++;}
      ptr++;
    }

  return(tmp);
}


/******************************************************************************
Function: my_sort().
Input : offset of list to be sorted, number of elements to be sorted.
Output: sorted list (in ascending order).
Return: none.
Description: stack based (nonrecursive) quicksort w/brute-shell bottom. 
******************************************************************************/
void
ivec_sort(register int *ar, register int size)
{
  register int *pi, *pj, temp;
  register int **top_a = offset_stack;
  register int *top_s = size_stack, *bottom_s = size_stack; 


  /* we're really interested in the offset of the last element */
  /* ==> length of the list is now size + 1                    */
  size--;

  /* do until we're done ... return when stack is exhausted */
  for (;;)
    {
      /* if list is large enough use quicksort partition exchange code */
      if (size > SORT_OPT)
	{	
	  /* start up pointer at element 1 and down at size     */  
	  pi = ar+1;
	  pj = ar+size;

	  /* find middle element in list and swap w/ element 1 */
	  SWAP(*(ar+(size>>1)),*pi)

	  pj = ar+size; 

	  /* order element 0,1,size-1 st {M,L,...,U} w/L<=M<=U */
	  /* note ==> pivot_value in index 0                   */
	  if (*pi > *pj)  
	    {SWAP(*pi,*pj)}
	  if (*ar > *pj) 
	    {SWAP(*ar,*pj)}
	  else if (*pi > *ar)   
	    {SWAP(*(ar),*(ar+1))}

	  /* partition about pivot_value ...  	                    */
	  /* note lists of length 2 are not guaranteed to be sorted */
	  for(;;)
	    {
	      /* walk up ... and down ... swap if equal to pivot! */
	      do pi++; while (*pi<*ar);
	      do pj--; while (*pj>*ar);

	      /* if we've crossed we're done */
	      if (pj<pi) break;

	      /* else swap */
	      SWAP(*pi,*pj)
	    }

	  /* place pivot_value in it's correct location */
	  SWAP(*ar,*pj)

	  /* test stack_size to see if we've exhausted our stack */
	  if (top_s-bottom_s >= MAX_STACK)
	    {
	      printf("\nSTACK EXHAUSTED!!!\n");
	      exit(1);
	    }

	  /* push right hand child iff length > 1 */
	  if (*top_s = size-(pi-ar))
	    {
	      *(top_a++) = pi;
	      size -= *top_s+2;  
	      top_s++;
	    }
	  /* set up for next loop iff there is something to do */
	  else if (size -= *top_s+2) 
	    {;}
	  /* might as well pop - note NR_OPT >=2 ==> we're ok! */
	  else
	    {
	      ar = *(--top_a);
	      size = *(--top_s);
	    }
	}

      /* else sort small list directly then pop another off stack */
      else
	{
	  /* insertion sort for bottom */
          for (pj=ar+1;pj<=ar+size;pj++)
            {
              temp = *pj;
              for (pi=pj-1;pi>=ar;pi--)
                {
                  if (*pi <= temp) break;
                  *(pi+1)=*pi;
                }
              *(pi+1)=temp;
	    }

	  /* check to see if stack is exhausted ==> DONE */
	  if (top_s==bottom_s) return;
	  
	  /* else pop another list from the stack */
	  ar = *(--top_a);
	  size = *(--top_s);
	}
    }
}



/******************************************************************************
Function: my_sort().
Input : offset of list to be sorted, number of elements to be sorted.
Output: sorted list (in ascending order).
Return: none.
Description: stack based (nonrecursive) quicksort w/brute-shell bottom. 
******************************************************************************/
void
ivec_sort_companion(register int *ar, register int *ar2, register int size)
{
  register int *pi, *pj, temp, temp2;
  register int **top_a = offset_stack;
  register int *top_s = size_stack, *bottom_s = size_stack; 
  register int *pi2, *pj2;
  register int mid;


  /* we're really interested in the offset of the last element */
  /* ==> length of the list is now size + 1                    */
  size--;

  /* do until we're done ... return when stack is exhausted */
  for (;;)
    {
      /* if list is large enough use quicksort partition exchange code */
      if (size > SORT_OPT)
	{	
	  /* start up pointer at element 1 and down at size     */  
	  mid = size>>1;
	  pi = ar+1;
	  pj = ar+mid;
	  pi2 = ar2+1;
	  pj2 = ar2+mid;

	  /* find middle element in list and swap w/ element 1 */
	  SWAP(*pi,*pj)
	  SWAP(*pi2,*pj2)

	  /* order element 0,1,size-1 st {M,L,...,U} w/L<=M<=U */
	  /* note ==> pivot_value in index 0                   */
	  pj = ar+size;
	  pj2 = ar2+size;
	  if (*pi > *pj)  
	    {SWAP(*pi,*pj) SWAP(*pi2,*pj2)}
	  if (*ar > *pj) 
	    {SWAP(*ar,*pj) SWAP(*ar2,*pj2)}
	  else if (*pi > *ar)   
	    {SWAP(*(ar),*(ar+1)) SWAP(*(ar2),*(ar2+1))}

	  /* partition about pivot_value ...  	                    */
	  /* note lists of length 2 are not guaranteed to be sorted */
	  for(;;)
	    {
	      /* walk up ... and down ... swap if equal to pivot! */
	      do {pi++; pi2++;} while (*pi<*ar);
	      do {pj--; pj2--;} while (*pj>*ar);

	      /* if we've crossed we're done */
	      if (pj<pi) break;

	      /* else swap */
	      SWAP(*pi,*pj)
	      SWAP(*pi2,*pj2)
	    }

	  /* place pivot_value in it's correct location */
	  SWAP(*ar,*pj)
	  SWAP(*ar2,*pj2)

	  /* test stack_size to see if we've exhausted our stack */
	  if (top_s-bottom_s >= MAX_STACK)
	    {
	      printf("\nSTACK EXHAUSTED!!!\n");
	      exit(1);
	    }

	  /* push right hand child iff length > 1 */
	  if (*top_s = size-(pi-ar))
	    {
	      *(top_a++) = pi;
	      *(top_a++) = pi2;
	      size -= *top_s+2;  
	      top_s++;
	    }
	  /* set up for next loop iff there is something to do */
	  else if (size -= *top_s+2) 
	    {;}
	  /* might as well pop - note NR_OPT >=2 ==> we're ok! */
	  else
	    {
	      ar2 = *(--top_a);
	      ar  = *(--top_a);
	      size = *(--top_s);
	    }
	}

      /* else sort small list directly then pop another off stack */
      else
	{
	  /* insertion sort for bottom */
          for (pj=ar+1, pj2=ar2+1;pj<=ar+size;pj++,pj2++)
            {
              temp = *pj;
              temp2 = *pj2;
              for (pi=pj-1,pi2=pj2-1;pi>=ar;pi--,pi2--)
                {
                  if (*pi <= temp) break;
                  *(pi+1)=*pi;
                  *(pi2+1)=*pi2;
                }
              *(pi+1)=temp;
              *(pi2+1)=temp2;
	    }

	  /* check to see if stack is exhausted ==> DONE */
	  if (top_s==bottom_s) return;
	  
	  /* else pop another list from the stack */
	  ar2 = *(--top_a);
	  ar  = *(--top_a);
	  size = *(--top_s);
	}
    }
}



/******************************************************************************
Function: my_sort().
Input : offset of list to be sorted, number of elements to be sorted.
Output: sorted list (in ascending order).
Return: none.
Description: stack based (nonrecursive) quicksort w/brute-shell bottom. 
******************************************************************************/
void
SMI_sort(void *ar1, void *ar2, int size, int type)
{
  if (type == INTEGER)
    {
      if (ar2)
	{ivec_sort_companion((int *)ar1,(int *)ar2,size);}
      else
	{ivec_sort(ar1,size);}
    }
  else
    {
      error_msg_fatal("SMI_sort only does INTEGER!");
    }
/*
  if (type == REAL)
    {
      if (ar2)
	{rvec_sort_companion(ar2,ar1,size);}
      else
	{rvec_sort(ar1,size);}
    }
*/
}



/********************************ivec.c**************************************
Function ivec_linear_search()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
int
ivec_linear_search(register int item, register int *list, register int n)
{
  register int tmp = n-1;

  while (n--)  {if (*list++ == item) {return(tmp-n);}}
  return(-1);
}




/********************************ivec.c**************************************
Function ivec_binary_search()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
int
ivec_binary_search(register int item, register int *list, register int rh)
{
  register int mid, lh=0;

  rh--;
  while (lh<=rh)
    {
      mid = (lh+rh)>>1;
      if (list[mid] == item) 
	{return(mid);}
      if (list[mid] < item)  
	{rh = mid-1;}
      else 
	{lh = mid+1;}
    }
  return(-1);
}






/********************************ivec.c**************************************
Function ivec_add()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
rvec_add(register double *arg1, register double *arg2, register int n)
{
  while (n--)  {*arg1++ += *arg2++;}
}



/********************************ivec.c**************************************
Function ivec_copy()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
rvec_copy(register double *arg1, register double *arg2, register int n)
{
  while (n--)  {*arg1++ = *arg2++;}
}



/********************************ivec.c**************************************
Function ivec_zero()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
void 
rvec_zero(register double *arg1, register int n)
{
  while (n--)  {*arg1++ = 0.0;}
}



/******************************************************************************
Function: perm_malloc()


add ability to pass in later ...

Space to be passed later must be FLOAT aligned!!!
******************************************************************************/
void 
perm_init(void)
{
  perm_req = 0;
  num_perm_req = 0;
  num_perm_frees = 0;

#ifdef MYMALLOC
  perm_top = perm_buf;
#endif
}



/******************************************************************************
Function: perm_malloc()

******************************************************************************/
void *
perm_malloc(size_t size)
{
  void *tmp;
  #ifdef MYMALLOC  
  double *space;
  int num_blocks;
  #endif

  
#ifdef DEBUG
      gsync();
      ierror_msg_warning("perm_malloc :: malloc!", size);
      gsync();
#endif

  if (!size)
    {
      error_msg_warning("perm_malloc() :: asking for no space?");
      return(NULL);
    }
     
  #ifdef MYMALLOC
  if (size%sizeof(double))
    {num_blocks = size/sizeof(double) + 1;}
  else
    {num_blocks = size/sizeof(double);}
  if (num_blocks < (PERM_MALLOC_BUF/sizeof(double) - (perm_top - perm_buf)))
    {
      space = perm_top;
      perm_top += num_blocks;
      perm_req+=size;
      num_perm_req++;
      return(space);
    }
  else
    {error_msg_fatal("perm_malloc() :: not enough space to satisfy request");}
  #endif


  tmp = (void *) malloc(size);

#ifdef DEBUG
      gsync();
      ierror_msg_warning("bss_malloc :: malloc!", (int) tmp);
      gsync();
#endif

  if (!tmp&&size)
    {
      printf("PERM ERROR :: %d:: %d,%d\n",my_id,(int) tmp,size);
      error_msg_fatal("Not that much memory left in PERM Pool!?!");
    }

  perm_req+=size;
  num_perm_req++;
  return(tmp);
}


/******************************************************************************
Function: perm_malloc_rem()


******************************************************************************/
int
perm_rem(void)
{
  error_msg_warning_n0("perm_rem() has been disabled!!!");
  return(MAX_INT);
}



/******************************************************************************
Function: perm_malloc()

******************************************************************************/
void 
perm_free(void *ptr)
{
  num_perm_frees--;  

  if (ptr)
    {
      #ifdef MYMALLOC
      return;
      #endif

      free((void *) ptr);
    }
  else
    {error_msg_warning("perm_free() :: nothing to free!!!");}
}



/******************************************************************************
Function: bss_malloc()

******************************************************************************/
void 
perm_stats(void)
{
  int min, max, ave, work;
  
  min = max = ave = perm_req;
#ifdef NXSRC
  gisum(&ave,1,&work);
  ave /= num_nodes;

  gilow(&min,1,&work);
  gihigh(&max,1,&work);
#endif
#ifdef MPISRC
  MPI_Allreduce (&ave, &work, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  ave = work/num_nodes;

/* Maybe needs a synchronization here to ensure work is not corrupted */
  MPI_Allreduce (&min, &work, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  min = work;
/* Maybe needs a synchronization here to ensure work is not corrupted */
  MPI_Allreduce (&max, &work, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  max = work;
#endif  

  if (!my_id)
    {
      printf("\n");
      printf("perm_malloc stats:\n");
      printf("perm_req      = %d\n",perm_req);
      printf("perm_req min  = %d\n",min);
      printf("perm_req ave  = %d\n",ave);
      printf("perm_req max  = %d\n",max);
    }

  min = max = ave = num_perm_req;
#ifdef NXSRC
  gisum(&ave,1,&work);
  gilow(&min,1,&work);
  gihigh(&max,1,&work);
#endif
#ifdef MPISRC
  MPI_Allreduce (&ave, &work, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  ave = work;
/* Maybe needs a synchronization here to ensure work is not corrupted */
  MPI_Allreduce (&min, &work, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  min = work;
/* Maybe needs a synchronization here to ensure work is not corrupted */
  MPI_Allreduce (&max, &work, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  max = work;
#endif
  if (!my_id)
    {
      printf("num_frees     = %d\n",num_perm_frees);
      printf("num_calls     = %d\n",num_perm_req);
      printf("num_calls min = %d\n",min);
      printf("num_calls ave = %.0f\n",ave/(1.0*num_nodes));
      printf("num_calls max = %d\n",max);
      
      printf("\n\n");
    }
  fflush(stdout);
}



/******************************************************************************
Function: bss_malloc()


add ability to pass in later ...

Space to be passed later must be FLOAT aligned!!!
******************************************************************************/
void 
bss_init(void)
{
  bss_req = 0;
  num_bss_req = 0;
  num_bss_frees = 0;

#ifdef MYMALLOC
  bss_top = bss_buf;
#endif
}



/******************************************************************************
Function: bss_malloc()

******************************************************************************/
void *
bss_malloc(size_t size)
{
  void *tmp;
  int max, work;
  #ifdef MYMALLOC  
  double *space;
  int num_blocks;
  #endif
 
#ifdef DEBUG
      gsync();
      ierror_msg_warning("bss_malloc :: malloc!", size);
      max = (int) size;
      gihigh((long)&max, 1, (long) &work);
      ierror_msg_warning("bss_malloc :: malloc!", max);
      gsync();
#endif
 

  #ifdef MYMALLOC
  if (size%sizeof(double))
    {num_blocks = size/sizeof(double) + 1;}
  else
    {num_blocks = size/sizeof(double);}
  if (num_blocks < (BSS_MALLOC_BUF/sizeof(double) - (bss_top - bss_buf)))
    {
      space = bss_top;
      bss_top += num_blocks;
      bss_req+=size;
      num_bss_req++;
      return(space);
    }
  else
    {error_msg_fatal("bss_malloc() :: not enough space to satisfy request");}
  #endif


/*     tmp = (void *) calloc(1,size);  */
  tmp = (void *) malloc(size); 
  
#ifdef DEBUG
      gsync();
      ierror_msg_warning("bss_malloc :: malloc!", (int) tmp);
      gsync();
#endif

  printf("%d:: %d,%d\n",my_id,(int) tmp,size);
  
  if (!tmp&&size)
    {
      printf("ERROR :: %d :: (%d,%d)\n",my_id,(int) tmp,size);
      error_msg_fatal("Not that much memory left in BSS Pool!?!");
    }

  bss_req+=size;
  num_bss_req++;
  return(tmp);
}


/******************************************************************************
Function: bss_malloc_rem()


******************************************************************************/
int
bss_rem(void)
{
  error_msg_warning_n0("bss_rem() has been disabled!!!");
  return(MAX_INT);
}



/******************************************************************************
Function: bss_malloc()

******************************************************************************/
void 
bss_free(void *ptr)
{
  num_bss_frees--;
  if (ptr)
    {
      #ifdef MYMALLOC
      return;
      #endif
      
      free((void *) ptr);
    }
  else
    {error_msg_warning("bss_free() :: nothing to free!!!");}
}



/******************************************************************************
Function: bss_malloc()

******************************************************************************/
void 
bss_stats(void)
{
  int min, max, ave, work;
  
  min = max = ave = bss_req;
#ifdef NXSRC
  gisum(&ave,1,&work);
  ave /= num_nodes;

  gilow(&min,1,&work);
  gihigh(&max,1,&work);
#endif
#ifdef MPISRC
  MPI_Allreduce (&ave, &work, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  ave = work/num_nodes;

/* Maybe needs a synchronization here to ensure work is not corrupted */
  MPI_Allreduce (&min, &work, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  min = work;
/* Maybe needs a synchronization here to ensure work is not corrupted */
  MPI_Allreduce (&max, &work, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  max = work;
#endif  

  if (!my_id)
    {
      printf("\n");
      printf("bss_malloc stats:\n");
      printf("bss_req      = %d\n",bss_req);
      printf("bss_req min  = %d\n",min);
      printf("bss_req ave  = %d\n",ave);
      printf("bss_req max  = %d\n",max);
    }

  min = max = ave = num_bss_req;
#ifdef NXSRC
  gisum(&ave,1,&work);
  gilow(&min,1,&work);
  gihigh(&max,1,&work);
#endif
#ifdef MPISRC
  MPI_Allreduce (&ave, &work, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  ave = work;
/* Maybe needs a synchronization here to ensure work is not corrupted */
  MPI_Allreduce (&min, &work, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  min = work;
/* Maybe needs a synchronization here to ensure work is not corrupted */
  MPI_Allreduce (&max, &work, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  max = work;
#endif

  if (!my_id)
    {
      printf("num_frees     = %d\n",num_bss_frees);
      printf("num_calls     = %d\n",num_bss_req);
      printf("num_calls min = %d\n",min);
      printf("num_calls ave = %.0f\n",ave/(1.0*num_nodes));
      printf("num_calls max = %d\n",max);
      printf("\n\n");

    }
  fflush(stdout);
}




/********************************driver.c**************************************
Function error_msg_fatal()

Input : pointer to error message.
Output: prints message to stdio, ofp.
Return: none
Description: prints string passed in and initiated program termination.
*********************************driver.c*************************************/
void 
error_msg_fatal(char *msg)
{
  /* print error message along w/node identifier */
  printf("\nFATAL :: P# %d :: %s\n", my_id,msg);
  fflush(stdout);

  /* exit program */
#ifdef MPISRC
/* Try with MPI_Finalize() as well _only_ if all procs call this routine */
/* Choose a more meaningful error code than -12 */
  MPI_Abort(MPI_COMM_WORLD, -12);
#endif
#ifdef NXSRC
  abort();
#endif
}



/********************************driver.c**************************************
Function error_msg_warning()

Input : pointer to error message.
Output: prints message to stdio.
Return: none.
Description: prints warning to stdio. Program continues ...
*********************************driver.c*************************************/
void 
error_msg_warning(char *msg)
{
#ifdef INFO
  /* print error message along w/node identifier */
  printf("\nWARNING :: P# %d :: %s\n", my_id,msg);
  fflush(stdout);

  /* do not exit program */
#endif
}



/********************************driver.c**************************************
Function error_msg_warning_n0()

Input : pointer to error message.
Output: prints message to stdio.
Return: none.
Description: prints warning to stdio. Program continues ...
*********************************driver.c*************************************/
void 
error_msg_warning_n0(char *msg)
{
#ifdef INFO
  /* if master node (node 0) then print error message */
  if (!my_id)
    {
      printf("\nMASTER WARNING :: P# %d :: %s\n",my_id,msg);
      fflush(stdout);
    }

  /* do not exit program */
#endif
}



/********************************driver.c**************************************
Function ierror_msg_warning_n0()

Input : pointer to error message.
Output: prints message to stdio.
Return: none.
Description: prints warning to stdio. Program continues ...
*********************************driver.c*************************************/
void 
ierror_msg_warning_n0(char *msg, int val)
{
#ifdef INFO
  /* if master node (node 0) then print error message */
  if (!my_id)
    {
      printf("\nMASTER WARNING :: P# %d :: %s :: val=%d\n",my_id,msg,val);
      fflush(stdout);
    }

  /* do not exit program */
#endif
}


void 
ierror_msg_warning(char *msg, int val)
{
#ifdef INFO
  printf("\nMASTER WARNING :: P# %d :: %s :: val=%d\n",my_id,msg,val);
  fflush(stdout);

  /* do not exit program */
#endif
}




/********************************driver.c**************************************
Function ferror_msg_warning_n0()

Input : pointer to error message.
Output: prints message to stdio.
Return: none.
Description: prints warning to stdio. Program continues ...
*********************************driver.c*************************************/
void 
ferror_msg_warning_n0(char *msg, float val)
{
#ifdef INFO
  /* if master node (node 0) then print error message */
  if (!my_id)
    {
      printf("\nMASTER WARNING :: P# %d :: %s :: val=%f\n",my_id,msg,val);
      fflush(stdout);
    }

  /* do not exit program */
#endif
}


/********************************driver.c**************************************
Function bm_to_proc

Input : 
Output: 
Return: 
Description: 
*********************************driver.c*************************************/
void 
bm_to_proc(register char *ptr, int p_mask, register int *msg_list)
{
  register int i, tmp;

  if (msg_list)
    {
      /* low to high */
      ptr+=(p_mask-1);
      for (i=p_mask-1;i>=0;i--)
	{
	  tmp = BYTE*(p_mask-i-1);
	  if (*ptr&BIT_0) 
	    {*msg_list = tmp; msg_list++;}
	  if (*ptr&BIT_1) 
	    {*msg_list = tmp+1; msg_list++;}
	  if (*ptr&BIT_2) 
	    {*msg_list = tmp+2; msg_list++;}
	  if (*ptr&BIT_3) 
	    {*msg_list = tmp+3; msg_list++;}
	  if (*ptr&BIT_4)
	    {*msg_list = tmp+4; msg_list++;}
	  if (*ptr&BIT_5)
	    {*msg_list = tmp+5; msg_list++;}
	  if (*ptr&BIT_6)
	    {*msg_list = tmp+6; msg_list++;}
	  if (*ptr&BIT_7) 
	    {*msg_list = tmp+7; msg_list++;}
	  ptr --;
	}

      /* high to low */
      /*
      for (i=0;i<p_mask;i++)
	{
	  tmp = BYTE*(p_mask-i-1);
	  if (*ptr&128) 
	    {*msg_list = tmp+7; msg_list++;}
	  if (*ptr&64)
	    {*msg_list = tmp+6; msg_list++;}
	  if (*ptr&32)
	    {*msg_list = tmp+5; msg_list++;}
	  if (*ptr&16)
	    {*msg_list = tmp+4; msg_list++;}
	  if (*ptr&8) 
	    {*msg_list = tmp+3; msg_list++;}
	  if (*ptr&4) 
	    {*msg_list = tmp+2; msg_list++;}
	  if (*ptr&2) 
	    {*msg_list = tmp+1; msg_list++;}
	  if (*ptr&1) 
	    {*msg_list = tmp; msg_list++;}
	  ptr ++;
	}
      */

    }
}



/********************************driver.c**************************************
Function ct_bits()

Input : 
Output: 
Return: 
Description: 
*********************************driver.c*************************************/
int 
ct_bits(register char *ptr, int n)
{
  register int i, tmp=0;


  for(i=0;i<n;i++)
    {
      if (*ptr&128) {tmp++;}
      if (*ptr&64)  {tmp++;}
      if (*ptr&32)  {tmp++;}
      if (*ptr&16)  {tmp++;}
      if (*ptr&8)   {tmp++;}
      if (*ptr&4)   {tmp++;}
      if (*ptr&2)   {tmp++;}
      if (*ptr&1)   {tmp++;}
      ptr++;
    }

  return(tmp);
}



/********************************driver.c**************************************
Function len_buf()

Input : 
Output: 
Return: 
Description: 
*********************************driver.c*************************************/
int
div_ceil(register int numin, register int denom)
{
  register int rt_val = 0;

  if ((numin<0)||(denom<=0))
    {error_msg_fatal("div_ceil() takes numin>0 denon>=0");}

  /* if integer division remainder then increment */
  rt_val = numin/denom;
  if (numin%denom) 
    {rt_val++;}
  
  return(rt_val);
}



/********************************driver.c**************************************
Function set_bit_mask()

Input : 
Output: 
Return: 
Description: 
*********************************driver.c*************************************/
void
set_bit_mask(register int *bm, int len, int val)
{
  register int i, offset;
  register char mask = 1;
  char *cptr;


  if (len_bit_mask(val)>len)
    {error_msg_fatal("The Bit Mask Isn't That Large!");}

  cptr = (char *) bm;

  offset = len/INT_LEN;
  for (i=0;i<offset;i++)
    {*bm=0; bm++;}

  offset = val%BYTE;
  for (i=0;i<offset;i++)
    {mask <<= 1;}

  offset = len - val/BYTE - 1;
  cptr[offset] = mask;
}


/********************************driver.c**************************************
Function len_bit_mask()

Input : 
Output: 
Return: 
Description: 
*********************************driver.c*************************************/
int
len_bit_mask(register int num_items)
{
  register int rt_val, tmp;

  if (num_items<0)
    {error_msg_fatal("Value Sent To len_bit_mask() Must be >= 0!");}

  /* mod BYTE ceiling function */
  rt_val = num_items/BYTE;
  if (num_items%BYTE) 
    {rt_val++;}
  
  /* make mults of sizeof int */
  if (tmp=rt_val%INT_LEN) 
    {return(rt_val+=(INT_LEN-tmp));}

  return(rt_val);
}



/********************************driver.c**************************************
Function len_buf()

Input : 
Output: 
Return: 
Description: 
*********************************driver.c*************************************/
int
len_buf(int item_size, int num_items)
{
  register int rt_val, tmp;

  rt_val = item_size * num_items;

  /*  byte align for now ... consider page later */
  if (tmp = (rt_val%FLOAT_LEN))
    {rt_val += (FLOAT_LEN - tmp);}

  return(rt_val);
}


int 
gs_init_ (int *elms, int *nel, int *level)
{
  gs_handle =  gs_init(elms, *nel, *level);
  return((int) gs_handle);
}



void
gs_gop_ (int *gs, double *vals)
{
  char *op = "+";

/*
  if ((gs_id *) *gs != gs_handle)
    {error_msg_warning_n0("hmmm ... gs_gop_ ()");}
*/
  
  gs_gop((gs_id *) *gs, vals, op);

/*
  gs_gop(gs_handle, vals, op);
*/
}



void
gs_free_ (int *gs)
{
/*
  if ((gs_id *) *gs != gs_handle)
    {error_msg_warning_n0("hmmm ... gs_free_ ()");}
*/

  gs_free((gs_id *) *gs);

/*
  gs_free(gs_handle);
*/
}



/******************************************************************************
Function: gs_init()

Input : 

Output: 

RETURN: 

Description:  
******************************************************************************/
gs_id *
gs_init(register int *elms, int nel, int level)
{
  register gs_id *gs;

  if (!p_init) {c_init();}

  /* determines if we have enough dynamic/semi-static memory */
  /* checks input, allocs and sets gd_id template            */
  gs = gsi_check_args(elms,nel,level);

  /* only bit mask version up and working for the moment    */
  /* LATER :: get int list version working for sparse pblms */
  gsi_via_bit_mask(gs);

  #ifdef INFO
/*  gs_print_template(gs,my_id);  
  bss_stats();
  perm_stats();
*/
  #endif

  /* clean up */
  bss_free((void *) gs->local_elms);
  bss_free((void *) gs->companion);
  bss_free((void *) gs->elms);
  bss_free((void *) gs->ngh_buf);
  gs->local_elms = gs->companion = gs->elms = gs->ngh_buf = NULL;

  /* print out P0's template as well as malloc stats */
  #ifdef INFO
/*   gs_print_template(gs,my_id); */
  bss_stats();
  perm_stats();

/*
  printf("%d :: gs = %d\n",my_id,(int) gs);
  error_msg_warning_n0("gs_init done!!!");
*/
  
  #endif

  return(gs);
}



/******************************************************************************
Function: gsi_check_args()

Input : 
Output: 
Return: 
Description: 

elm list must >= 0!!!
elm repeats allowed
******************************************************************************/
static
gs_id *
gsi_check_args(int *in_elms, int nel, int level)
{
  register int i=0, j, k;
  register gs_id *gs;
  int t2;
  int *vals, *work, *oprs;
  int *companion, *elms, *unique;
  int num_local=0, *num_to_reduce, **local_reduce;
  int *iptr;
  

  if (!in_elms)
    {error_msg_fatal("elms point to nothing!!!");}

  if (nel < 0)
    {error_msg_fatal("can't have fewer than 0 elms!!!");}

  if (nel == 0)
    {error_msg_warning("I don't have any elements!!!");}


  /* get space for gs template */
  gs = (gs_id *) perm_malloc(sizeof(gs_id));
  gs->id = ++num_gs_ids;

  /* default parameters */
  gs->repeats = FALSE;
  gs->ordered = TRUE;
  gs->positive = TRUE;

  /* copy over in_elms list and create inverse map */
  elms = (int *) bss_malloc((nel+1)*INT_LEN);
  companion = (int *) bss_malloc(nel*INT_LEN);
  for (i=0;i<nel;i++)
    {companion[i] = i; elms[i] = in_elms[i];}
  elms[i] = MIN_INT;
  SMI_sort((void *)elms, (void *)companion, nel, INTEGER);

  /* determine number of unique elements, check pd */
  i=0;
  while(i<nel)
    {
      j=i+1;
      t2 = elms[i];
      
      if (t2 <= 0)
	{gs->positive = FALSE;}
     
      /* clump 'em for now */ 
      while ((j<nel)&&(elms[j]==t2))
	{j++;}
      
      /* how many together and num local */
      j -= i;
      if (j>1)
	{
	  num_local++;
	  gs->repeats += (j-1);
	}
      i+=j;
    }

  /* two extras for sentinal and partition ... */
  gs->num_local = num_local;
  num_local+=2;
  gs->local_reduce=local_reduce=(int **)perm_malloc(num_local*INT_PTR_LEN);
  gs->num_local_reduce=num_to_reduce=(int *) perm_malloc(num_local*INT_LEN);
  gs->nel = (nel - gs->repeats);
  unique = (int *) bss_malloc(gs->nel*INT_LEN);

  /* compess map as well as keep track of local ops */
  num_local = 0;
  for (i=0, j=0;i<gs->nel;i++)
    {
      t2 = unique[i] = elms[j];
      companion[i] = companion[j];
     
      k=j;
      while ((j<nel)&&(elms[j]==t2))
	{j++;}

      if ((j-k)>1)
	{
	  num_to_reduce[num_local] = (j-k);
	  iptr = local_reduce[num_local] = (int *)perm_malloc((j-k+1)*INT_LEN);
	  while (k<j)
	    {*(iptr++) = companion[k++];}
	  *iptr = -1;
	  num_local++;
	}
    }
  num_to_reduce[num_local] = 0;
  local_reduce[num_local] = NULL;
  if (num_local!=gs->num_local)
    {error_msg_fatal("compression of maps wrong!");}

  num_to_reduce[++num_local] = 0;
  local_reduce[num_local] = NULL;

  gs->elms = unique; 
  gs->nel_total = nel;
  gs->local_elms = elms;
  gs->companion = companion;


  /* gather some global problem info */
  /* buffer included from c_init.e   */
  /* min size check done in c_init.c */
  vals = (int *) boot_msg_buf;
  work = (vals + 11);
  oprs = (work + 11);


  /* load 'em up */
  /* note one extra to hold NON_UNIFORM flag!!! */
  oprs[0] = NON_UNIFORM;
  vals[0] = perm_rem(); 
  oprs[1] = GL_MIN; 
  vals[1] = bss_rem();  
  oprs[2] = GL_MIN;
  vals[2] = nel;
  oprs[3] = GL_MIN;
  vals[3] = ivec_lb(elms,nel); 
  oprs[4] = GL_MIN;
  vals[4] = nel;
  oprs[5] = GL_MAX;
  vals[5] = ivec_ub(elms,nel); 
  oprs[6] = GL_MAX;
  vals[6] = nel;  
  oprs[7] = GL_ADD;
  vals[7] = level;
  oprs[8] = GL_B_AND;
  vals[8] = num_gs_ids;
  oprs[9] = GL_B_AND;
  vals[9] = gs->ordered;
  oprs[10] = GL_B_AND;
  vals[10] = gs->positive;
  oprs[11] = GL_B_AND;

  /* GLOBAL: send 'em out */
  giop(vals,work,11,oprs);

  /* unpack 'em */
  /* ordered ? */
  if (!vals[9])
    {error_msg_warning_n0("System not ordered!");}

  /* positive ? */
  if (!vals[10])
    {error_msg_warning_n0("System not positive!");}
  if (vals[3] <= 0)
    {error_msg_warning_n0("System not positive!");}

  /* check gs template count */
  if (vals[8] != num_gs_ids)
    {error_msg_fatal("num_gs_ids mismatch!!!");}

  /* check all have same level threshold */
  if (level != vals[7])
    {error_msg_fatal("all must call gs_init w/same level!!!");}

  /* check buffer space existence and size */
/*
  if ((vals[0] < (2*N_HALF)) || (vals[1] < (2*N_HALF)))
    {rt_val += NOT_ENOUGH_MEM;}
*/

  /* if any error nuke 'em */
/* 
   if (rt_val)
    {return(rt_val);}
*/

  /*  gs->nel = nel; */
  gs->nel_max = vals[4];
  gs->nel_min = vals[2];
  gs->nel_sum = vals[6];
  gs->negl =  (vals[5] - vals[3]) + 1;
  gs->gl_max = vals[5];
  gs->gl_min = vals[3];
  gs->vals = NULL;
  gs->gl_bss_min = vals[0];
  gs->gl_perm_min = vals[1];

  /* LATER :: add level == -1 -> program selects level */
  if (vals[7]<0)
    {vals[7]=0;}
  if (vals[7]>num_nodes)
    {vals[7]=num_nodes;}
  gs->level = vals[7];

  gs->num_pairs = 0;
  gs->max_pairs = 0;
  gs->pair_list = NULL;  
  gs->msg_sizes = NULL;
  gs->node_list = NULL;
 

  /*  
  for (i=0;i<(nel-1);i++)
    {
      if (elms[i] > elms[i+1]) 
	{
	  gs->ordered = FALSE;
	  error_msg_fatal("local elms  must be asc ordered for now!");
	}
      else if ((elms[i] <= 0) || (elms[i+1] <= 0))
	{error_msg_fatal("local elms  must be all greater than zero!");}
      if (elms[i] == elms[i+1]) 
	{
	  gs->repeats = TRUE;
	  error_msg_warning("No repeats allowed!");
	}
    }
  */

  /* LATER: do crystal accumulator init */
  gs->max_left_over = 0;
  gs->in_num    = NULL;
  gs->out_num   = NULL;
  gs->in_list   = NULL;
  gs->out_list  = NULL;


  /* ngh_buf is data drives the entire process */
  gs->ngh_buf_sz = 0;
  gs->ngh_buf    = NULL;
  gs->neighbors = NULL;


  /* LATER: right now accept only pos ordered (asc) lists w/no repeats!!!   */
  /* to fix: copy elms w/shift, compression and sorting -> {0, ..., nel'-1} */
  /* after gs_gop the reverse compression and sorting, no need to shift     */

  /* Relation: Thresh =  (negl * (ceil(#P/32) -1))        */
  /* if (C_IDA * #S(hared)) < Thresh then int list alg    */
  /* else bit_mask alg                                    */
  /* three ideas on getting #S(hared)                     */
  /* 1). let #S = (nel_sum-negl). Fine if all glnums used */  
  /* 2). let #S = negl * %shared up to glnum(i) < negl_max*/  
  /* 3). do addition on all glnums vector. Exact but EXP  */  
  /* (nel_sum > (negl*(1-10*(div_ceil(num_nodes,32)-1)))) */

  /* LATER: add this and list alg */
  /*
  if (gs->nel_sum<(gs->negl*(1.0+0.2*(div_ceil(num_nodes,32)*4-1))))
      {return(INT_LIST);}
  */

  /* only bit routine working now! */
  /*  error_msg_warning_n0("Step ... finished gs_init()!"); */

  return(gs);
}



/******************************************************************************
Function: gsi_via_bit_mask()

Input : 
Output: 
Return: 
Description: 


******************************************************************************/
static
void
gsi_via_bit_mask(gs_id *gs)
{
  register int i, nel, *elms;
  int t1,t2,op;
  int **reduce;
  int *map;


  /* totally local removes ... ct_bits == 0 */
  get_ngh_buf(gs);
  
#ifdef DEBUG
  printf("gsi_via_bit_mask() :: gs->level = %d\n",gs->level);
  error_msg_warning("Steping into set_pairwise()"); 
#endif


  
  if (gs->level)
    {set_pairwise(gs);}

#ifdef DEBUG
  error_msg_warning("Steping out of  set_pairwise()"); 
#endif

  t1 = 0;
  elms = gs->elms;
  nel = gs->nel;
  for (i=0;i<nel;i++)
    {
      if (!(elms[i] & BIT_31))
	{t1++;}
    }

  op = GL_MAX;
  giop(&t1,&t2,1,&op);
  
  if (gs->max_left_over = t1)
    {set_tree(gs);}

  /* reset tmp->elms to positive */
  for (i=0;i<nel;i++)
    {
      if (*elms & BIT_31)
	{*elms ^= BIT_31;}
      elms++;
    }

  /* remap pairwise */
  map = gs->companion;
  elms = gs->pw_elm_list;
  for (i=0; i<gs->len_pw_list; i++)
    {elms[i] = map[elms[i]];}
      
  
  /* LATER ::: remap tree */

  /* intersection local and pairwise/tree? */
  gs->num_local_total = gs->num_local;
  if (gs->num_local == 0)
    {
      gs->local_strength = NONE;
      gs->num_local_gop = 0;
      gs->num_local_total =  0;
      gs->gop_local_reduce = gs->local_reduce;
      gs->num_gop_local_reduce = gs->num_local_reduce;
    }
  else
    {
      reduce = gs->local_reduce;  
      for (i=0, t1=0; i<gs->num_local; i++, reduce++)
	{
	  /* LATER :: also must search tree list when done!!! */
	  if (ivec_linear_search(**reduce,gs->pw_elm_list,gs->len_pw_list)>=0)
	    {t1++; gs->num_local_reduce[i] *= -1;}
	}

      if (!t1)
	{gs->local_strength = FULL;}
      else
	{
	  gs->local_strength = PARTIAL;
	  SMI_sort((void *)gs->num_local_reduce, (void *)gs->local_reduce, 
		           gs->num_local + 1, INTEGER);

	  gs->num_local_gop = t1;
	  gs->num_local_total =  gs->num_local;
	  gs->num_local    -= t1;
	  gs->gop_local_reduce = gs->local_reduce;
	  gs->num_gop_local_reduce = gs->num_local_reduce;

	  
	  for (i=0; i<t1; i++)
	    {
	      gs->num_gop_local_reduce[i] *= -1;
	      gs->local_reduce++;
	      gs->num_local_reduce++;
	    }
	  gs->local_reduce++;
	  gs->num_local_reduce++;
	}
    }
#ifdef DEBUG
  gsync();
  error_msg_warning_n0("Done via_pairwise()"); 
#endif
}



/******************************************************************************
Function: get_ngh_buf()

Input : 
Output: 
Return: 
Description: 


******************************************************************************/
static
void
get_ngh_buf(gs_id *gs)
{
  register int i, j;
  int p_mask_size, ngh_buf_size, buf_size;
  int *p_mask, *sh_proc_mask;
  int *ngh_buf, *buf1, *buf2;
  int offset, per_load, num_loads, or_ct, start, end;
  int *msg_list;
  int *ptr1, *ptr2, i_start, negl, nel, *elms;
  int oper[1];


  /* to make life easier and faster */
  nel  = gs->nel;
  negl = gs->negl;
  elms = gs->elms;

  /* det #bytes needed for processor bit masks and init */
  p_mask_size   = len_bit_mask(num_nodes);
  p_mask        = (int *) bss_malloc(p_mask_size);
  set_bit_mask(p_mask,p_mask_size,my_id);

  /* allocate space for additional masks and init */
  gs->neighbors = sh_proc_mask = (int *) perm_malloc(p_mask_size);

  gs->ngh_buf_sz = ngh_buf_size = p_mask_size*nel;
  gs->ngh_buf = ngh_buf = (int *) bss_malloc(ngh_buf_size);

  /* give remaining amount of globally available bss space to msg buf */
  /* determine min msg buffer size across all nodes */
  /*  
  buf1 = (int *) bss_malloc(MAX_MSG_BUF<<1); 
  buf_size = ivec_split_buf(buf1,&buf2,MAX_MSG_BUF<<1);
  */
  buf1 = (int *) bss_malloc(MAX_MSG_BUF);
  buf2 = (int *) bss_malloc(MAX_MSG_BUF);
  buf_size = MAX_MSG_BUF;


  /* do we have the space to do the job */
  if ((!p_mask)||(!sh_proc_mask)||(!buf1)||(buf_size<p_mask_size))
    {error_msg_fatal("get_ngh_buf() :: Malloc Failure Or Buf Too Sm!");}
  if (!ngh_buf && nel)
    {error_msg_fatal("get_ngh_buf() :: Ngh_buf Malloc Failure!");}

  /* are we about to overwork ourselves? */
  if (buf_size < N_HALF)
    {error_msg_warning_n0("get_ngh_buf() :: Msg Buf smaller than N_HALF!");}

  /* convert buf sizes from bytes to bytes/INT_LEN */
  p_mask_size>>=2; ngh_buf_size>>=2; buf_size>>=2;
  /* p_mask_size/=INT_LEN; ngh_buf_size/=INT_LEN; buf_size/=INT_LEN; */
  gs->mask_sz = p_mask_size;

  /* init buffers and det number of gior exchanges */
  ivec_zero(sh_proc_mask,p_mask_size);
  ivec_zero(ngh_buf,ngh_buf_size);

  if (buf_size>=(p_mask_size*negl))
    {
      buf_size = p_mask_size*negl;
      per_load = negl;
      gs->num_loads = num_loads = 1;
    }
  else
    {
      per_load = buf_size/p_mask_size;
      gs->num_loads = num_loads = div_ceil(negl,per_load);
    }

  if (num_loads>1)
   {error_msg_warning_n0("gsi_via_bit_mask() :: Switch to int list setup?!?");}

  /* LATER test ivec_or vs.() ivec_copy() */
  /* ivec_or  (buf1+(elms[i]-start)*p_mask_size,p_mask,p_mask_size); i++;} */
  oper[0] = GL_B_OR;
  ptr1 = ngh_buf;
  ptr2 = elms;
  end = gs->gl_min;
  for (or_ct=0,i=0;or_ct<num_loads;or_ct++)
    {
      ivec_zero(buf1,buf_size);
      start = end;
      end+=per_load;
      i_start = i;

      /* load msg buffer */
      while ((i<nel)&&((offset=*ptr2)<end))
	{
	  offset = (offset-start)*p_mask_size;
	  ivec_copy(buf1+offset,p_mask,p_mask_size); 
	  i++; ptr2++;
	}

      /* GLOBAL: pass buffer */
      oper[0] = GL_B_OR;
      giop(buf1, buf2, buf_size, oper);



/*      gior((long *) buf1, buf_size, (long *) buf2); */


      ptr2=(elms+i_start);
      /* unload buffer into ngh_buf */
      for(j=i_start;j<i;j++) 
	{
	  offset = *ptr2;
	  offset = (offset-start)*p_mask_size;
	  ivec_copy(ptr1,buf1+offset,p_mask_size);
	  ptr1+=p_mask_size; ptr2++;
	}
    }

  /* unset processor mask bit and collect */
  ivec_zero(sh_proc_mask, p_mask_size);
  for (i=0,buf2=ngh_buf;i<nel;i++,buf2+=p_mask_size)
    {
      ivec_xor(buf2,p_mask,p_mask_size);
      ivec_or(sh_proc_mask,buf2,p_mask_size);

      /* elms not shared are complete */
      if (!ct_bits((char *) buf2, p_mask_size*INT_LEN))
	{elms[i] |= BIT_31;}
    }
  
  gs->num_nghs = ct_bits((char *)sh_proc_mask,p_mask_size*INT_LEN);

  bss_free((void *)p_mask);
  bss_free((void *)buf1);
  bss_free((void *)buf2);

#ifdef DEBUG
  error_msg_warning("returning from get_ngh_buf!");
  gsync();
#endif
}



/******************************************************************************
Function: pairwise_init()

Input : 
Output: 
Return: 
Description: 

if an element is shared by fewer that level# of nodes do pairwise exch 
******************************************************************************/
static
void
set_pairwise(gs_id *gs)
{
  register int i, j;
  int p_mask_size;
  int *p_mask, *sh_proc_mask, *tmp_proc_mask;
  int *ngh_buf, *buf2;
  int offset;
  int *msg_list, *msg_size, **msg_nodes, nprs;
  int *pairwise_elm_list, len_pair_list=0;
  int *iptr, t1, i_start, nel, *elms;
  int ct, level;


  /* to make life easier and faster */
  nel  = gs->nel;
  elms = gs->elms;
  ngh_buf = gs->ngh_buf;
  level = gs->level;


#ifdef DEBUG
      error_msg_warning("set_pairwise() :: asgns!");
#endif

  p_mask_size   = len_bit_mask(num_nodes);

#ifdef DEBUG
      ierror_msg_warning("set_pairwise() :: len_bit_mask!",p_mask_size);
#endif

  p_mask        = (int *) bss_malloc(p_mask_size);
  sh_proc_mask  = (int *) bss_malloc(p_mask_size);
  tmp_proc_mask = (int *) bss_malloc(p_mask_size);

#ifdef DEBUG
      error_msg_warning("set_pairwise() :: bss_malloc!");
#endif

  if ((!p_mask)||(!sh_proc_mask)||(!tmp_proc_mask))
    {error_msg_fatal("pairwise_init() :: malloc failure for masks!");}
  set_bit_mask(p_mask,p_mask_size,my_id);

#ifdef DEBUG
      error_msg_warning("set_pairwise() :: set_bit_mask!");
#endif

  p_mask_size >>= 2;
  /* p_mask_size /= INT_LEN; */

  /* find cover st all elms under cover have fewer than level nghs */
  ivec_zero(sh_proc_mask, p_mask_size);

#ifdef DEBUG
      error_msg_warning("set_pairwise() :: ivec_zero!");
#endif

  for (i=0,buf2=ngh_buf;i<nel;i++,buf2+=p_mask_size)
    {
      if (ct = ct_bits((char *)buf2,p_mask_size*INT_LEN))
	{
	  if (ct <= level)
	    {
	      ivec_or(sh_proc_mask,buf2,p_mask_size);
	      len_pair_list++;
	    }
	}
    }
	  
  gs->len_pw_list=len_pair_list;
  gs->pw_elm_list=pairwise_elm_list=perm_malloc((len_pair_list+1)*INT_LEN);

#ifdef DEBUG
      error_msg_warning("set_pairwise() :: perm_malloc!");
#endif

  /* how many processors (nghs) do we have to exchange with? */
  nprs=gs->num_pairs=ct_bits((char *)sh_proc_mask,p_mask_size*INT_LEN);
  
#ifdef DEBUG
  if (!nprs)
    {
      
      error_msg_warning("set_pairwise() :: nprs is zero!");
      /* must take into account giop and frees!!! */
/*
      bss_free((void *)p_mask);
      bss_free((void *)sh_proc_mask);
      bss_free((void *)tmp_proc_mask);
      return;
*/
    }
#endif

  /* allocate space for gs_gop() info */
  gs->pair_list = msg_list = (int *)  perm_malloc(INT_LEN*nprs);
  gs->msg_sizes = msg_size  = (int *)  perm_malloc(INT_LEN*nprs);
  gs->node_list = msg_nodes = (int **) perm_malloc(INT_PTR_LEN*(nprs+1));

  /* do we have the space to do the job - poss. NO nghs! */
  if ((nprs)&&((!msg_list)||(!msg_size)||(!msg_nodes)))
    {error_msg_fatal("set_pairwise :: Malloc Failure for info ptrs!");}

  /* init msg_size list */
  ivec_zero(msg_size,nprs);  

  /* expand from bit mask list to int list */
  bm_to_proc((char *)sh_proc_mask,p_mask_size*INT_LEN,msg_list);
  
  /* who are we going to take care of ... make neg to indicate */
  for (i=0,j=0,buf2=ngh_buf;i<nel;i++,buf2+=p_mask_size)
    {
      if (ct = ct_bits((char *)buf2,p_mask_size*INT_LEN))
	{
	  if (ct <= level)
	    {
	      pairwise_elm_list[j++] = i;
	      elms[i] |= BIT_31;
	    }
	}
    }
  pairwise_elm_list[j] = -1;
  if (j!=len_pair_list)
    {error_msg_fatal("oops ... bad paiwise list in set_pairwise!");}

#ifdef NXSRC
  gs->msg_ids_out = (int *)  perm_malloc(INT_LEN*(nprs+1));
  ivec_zero(gs->msg_ids_out,nprs);
  gs->msg_ids_out[nprs] = -1;
  gs->msg_ids_in = (int *)  perm_malloc(INT_LEN*(nprs+1));
  ivec_zero(gs->msg_ids_in,nprs);
  gs->msg_ids_in[nprs] = -1;
  gs->pw_vals = (double *) perm_malloc(sizeof(double)*len_pair_list);
#endif
#ifdef MPISRC
  gs->msg_ids_out = (MPI_Request *)  perm_malloc(sizeof(MPI_Request)*(nprs+1));
/* no need to initialize the rest of the entries to zero */
  gs->msg_ids_out[nprs] = MPI_REQUEST_NULL;
  gs->msg_ids_in = (MPI_Request *)  perm_malloc(sizeof(MPI_Request)*(nprs+1));
/* no need to initialize the rest of the entries to zero */
  gs->msg_ids_in[nprs] = MPI_REQUEST_NULL;
  gs->pw_vals = (double *) perm_malloc(sizeof(double)*len_pair_list);
#endif

  /* cover method not working for obv. reasons */
  /*
      if (!ct_bits((char *) buf2,p_mask_size*INT_LEN))
	{elms[i] |= BIT_31;}
      else
	{
	  ivec_or3(tmp_proc_mask, buf2, sh_proc_mask, p_mask_size);
	  if (ivec_cmp(tmp_proc_mask, sh_proc_mask, p_mask_size))
	    {elms[i] |= BIT_31;}
	}
    }
  */

  i_start = 0;
  for (i=0;i<nprs;i++)
    {
      set_bit_mask(p_mask,p_mask_size*INT_LEN,msg_list[i]);
      buf2 = ngh_buf;
      ct =0;
      for (j=0;j<nel;j++)
	{
	  if (elms[j] & BIT_31)
	    {
	      ivec_and3(tmp_proc_mask,p_mask,buf2,p_mask_size);
	      if (ct_bits((char *)tmp_proc_mask,p_mask_size*INT_LEN))
		{ct++;}
	    }
	  buf2+=p_mask_size;
	}

      /* set number of messages */
      msg_size[i] = ct;
      i_start = MAX(i_start,ct);

      /*space to hold nodes in message to first neighbor */
      msg_nodes[i] = iptr = (int *) perm_malloc(INT_LEN*(ct+1));

      /* do we have the space to do the job */
      /*      gs_print_template(gs->id,0); */
      if (!iptr)
	{error_msg_fatal("Malloc Failure for iptr In gs_init()!");}

      buf2 = ngh_buf;
      t1 =0;
      for (j=0;j<nel;j++)
	{
	  if (elms[j] & BIT_31)
	    {
	      ivec_and3(tmp_proc_mask,p_mask,buf2,p_mask_size);
	      if (ct_bits((char *)tmp_proc_mask,p_mask_size*INT_LEN))
		/*		{*(iptr+t1) = j; t1++; ct--;} */
		{
       *(iptr+t1)=ivec_linear_search(j,pairwise_elm_list,len_pair_list); 
		  t1++; 
		  ct--;
		} 
	    }
	  buf2+=p_mask_size;
	}
      *(iptr+t1) = -1;

      if (ct)
	{printf("\n Count should be zero\n"); fflush(stdout);}
    }
  msg_nodes[nprs] = NULL;

  /* build custom exch when we determine exactly what info gop() needs */
  /* LATER: more info to assist in tuning */
  /* LATER: who knows what info to collect so ... */
  /*  
  gs->out = (FLOAT *) perm_malloc(sizeof(double)*len_pair_list);
  gs->in  = (FLOAT *) perm_malloc(sizeof(double)*len_pair_list);
  */

  t1 = GL_MAX;
  giop(&i_start,&offset,1,&t1);
  gs->max_pairs = i_start;
  gs->msg_total = ivec_sum(gs->msg_sizes,nprs);
  gs->out = (double *) perm_malloc(sizeof(double)*gs->msg_total);
  gs->in  = (double *) perm_malloc(sizeof(double)*gs->msg_total);

  /* reset malloc pool */
  bss_free((void *)p_mask);
  bss_free((void *)sh_proc_mask);
  bss_free((void *)tmp_proc_mask);
}    



/******************************************************************************
Function: gsi_via_int_list()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
void
gsi_via_int_list(gs_id *gs)
{

  /* LATER: for P large the bit masks -> too many passes */
  /* LATER: strategy: do gsum w/1 in position i in negl if owner */
  /* LATER: then sum of entire vector 1 ... negl determines min buf len */
  /* LATER: So choose min from this or mask method */
}



/******************************************************************************
Function: gather_scatter

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
gs_gop(register gs_id *gs, register double *vals, register char *operation)
{
  if (!gs) {error_msg_fatal("gs_gop() passed NULL gs handle!!!");}

  /* local only operations!!! */
  if (gs->num_local)
    {gs_gop_local(gs,vals,operation);}

  /* if intersection tree/pairwise and local isn't empty */
  if (gs->num_local_gop)
    {
      gs_gop_local_in(gs,vals,operation);

      /* pairwise */
      if (gs->num_pairs)
	{gs_gop_pairwise(gs,vals,operation);}
      
      /* tree */
      if (gs->max_left_over)
	{gs_gop_tree(gs,vals,operation);}
  
      gs_gop_local_out(gs,vals,operation);
    }
  /* if intersection tree/pairwise and local is empty */
  else
    {
      /* pairwise */
      if (gs->num_pairs)
	{gs_gop_pairwise(gs,vals,operation);}
      
      /* tree */
      if (gs->max_left_over)
	{gs_gop_tree(gs,vals,operation);}
    }
}



/******************************************************************************
Function: gather_scatter

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
void
gs_gop_local(register gs_id *gs, register double *vals, char *operation)
{
  register int *num, *map, **reduce;
  register double tmp;


  num    = gs->num_local_reduce;  
  reduce = gs->local_reduce;  
  while (map = *reduce)
    {
      /* wall */
      if (*num == 2)
	{
	  num ++; reduce++;
	  vals[map[1]] = vals[map[0]] += vals[map[1]];
	}
      /* corner shared by three elements */
      else if (*num == 3)
	{
	  num ++; reduce++;
	  vals[map[2]]=vals[map[1]]=vals[map[0]]+=(vals[map[1]]+vals[map[2]]);
	}
      /* corner shared by four elements */
      else if (*num == 4)
	{
	  num ++; reduce++;
	  vals[map[1]]=vals[map[2]]=vals[map[3]]=vals[map[0]] += 
	                         (vals[map[1]] + vals[map[2]] + vals[map[3]]);
	}
      /* general case ... odd geoms ... 3D*/
      else
	{
	  num ++;
 	  tmp = 0.0;
	  while (*map >= 0)
	    {tmp += *(vals + *map++);}

	  map = *reduce++;
	  while (*map >= 0)
	    {*(vals + *map++) = tmp;}
	}
    }
}



/******************************************************************************
Function: gather_scatter

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
void
gs_gop_local_in(register gs_id *gs, register double *vals, char *operation)
{
  register int *num, *map, **reduce;
  register double *base;


  num    = gs->num_gop_local_reduce;  
  reduce = gs->gop_local_reduce;  
  while (map = *reduce++)
    {
      /* wall */
      if (*num == 2)
	{
	  num ++;
	  vals[map[0]] += vals[map[1]];
	}
      /* corner shared by three elements */
      else if (*num == 3)
	{
	  num ++;
	  vals[map[0]] += (vals[map[1]] + vals[map[2]]);
	}
      /* corner shared by four elements */
      else if (*num == 4)
	{
	  num ++;
	  vals[map[0]] += (vals[map[1]] + vals[map[2]] + vals[map[3]]);
	}
      /* general case ... odd geoms ... 3D*/
      else
	{
	  num++;
	  base = vals + *map++;
	  while (*map >= 0)
	    {*base += *(vals + *map++);}
	}
    }
}



/******************************************************************************
Function: gather_scatter

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
void
gs_gop_local_out(register gs_id *gs, register double *vals, char *operation)
{
  register int *num, *map, **reduce;
  register double tmp;


  num    = gs->num_gop_local_reduce;  
  reduce = gs->gop_local_reduce;  
  while (map = *reduce++)
    {
      /* wall */
      if (*num == 2)
	{
	  num ++;
	  vals[map[1]] = vals[map[0]];
	}
      /* corner shared by three elements */
      else if (*num == 3)
	{
	  num ++;
	  vals[map[2]] = vals[map[1]] = vals[map[0]];
	}
      /* corner shared by four elements */
      else if (*num == 4)
	{
	  num ++;
	  vals[map[3]] = vals[map[2]] = vals[map[1]] = vals[map[0]];
	}
      /* general case ... odd geoms ... 3D*/
      else
	{
	  num++;
	  tmp = *(vals + *map++);
	  while (*map >= 0)
	    {*(vals + *map++) = tmp;}
	}
    }
}



/******************************************************************************
Function: gather_scatter

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
void
gs_gop_tree(gs_id *gs, double *vals, char *operation)
{
  error_msg_fatal("gs_gop_tree() :: shouldn't be calling this routine!!!");
}



/******************************************************************************
Function: set_tree()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
void
set_tree(gs_id *gs)
{
  error_msg_fatal("tree not ready yet!!!");
}



/******************************************************************************
Function: level_best_guess()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
int 
level_best_guess(void)
{
  /* full pairwise for now */
  return(num_nodes);
}



/******************************************************************************
Function: gs_print_template()

Input : 

Output: 

Return: 

Description:  
******************************************************************************/
void
gs_print_template(register gs_id* gs, int who)
{
  register int i, j, k, *iptr, *iptr2;

  if (!num_gs_ids) return;

  for (i=0;i<num_nodes;i++)
    {
      gsync();
      
  if ((my_id == i)&&(my_id == who))
    {
      printf("\n\nP#%d's GS#%d template:\n", my_id, gs->id);
      printf("id=%d\n",          gs->id);
      printf("nel(unique)=%d\n", gs->nel);
      printf("nel_max=%d\n",     gs->nel_max);
      printf("nel_min=%d\n",     gs->nel_min);
      printf("nel_sum=%d\n",     gs->nel_sum);
      printf("negl=%d\n",        gs->negl);
      printf("gl_max=%d\n",      gs->gl_max);
      printf("gl_min=%d\n",      gs->gl_min);
      printf("elms ordered=%d\n",gs->ordered);
      printf("repeats=%d\n",     gs->repeats);
      printf("positive=%d\n",    gs->positive);
      printf("elms=%d\n",        gs->elms);
      printf("elms(total)=%d\n", gs->local_elms);
      printf("vals=%d\n",        gs->vals);
      printf("gl_bss_min=%d\n",  gs->gl_bss_min);
      printf("gl_perm_min=%d\n", gs->gl_perm_min);
      printf("level=%d\n",       gs->level);
      printf("proc_mask_sz=%d\n",gs->mask_sz);
      printf("sh_proc_mask=%d\n",gs->neighbors);
      printf("ngh_buf_size=%d\n",gs->ngh_buf_sz);
      printf("ngh_buf=%d\n",     gs->ngh_buf);
      printf("num_nghs=%d\n",    gs->num_nghs);

      /* pairwise exchange information */
      printf("\nPaiwise Info:\n");
      printf("num_pairs=%d\n",   gs->num_pairs);
      printf("pair_list=%d\n",   gs->pair_list);  
      printf("msg_sizes=%d\n",   gs->msg_sizes);
      printf("node_list=%d\n",   gs->node_list);
      printf("max_pairs=%d\n",   gs->max_pairs);
      printf("len_pw_list=%d\n", gs->len_pw_list);
      printf("pw_elm_list=%d\n", gs->pw_elm_list);

      printf("pw_elm_list: ");
      if (iptr = gs->pw_elm_list)
	{
	  for (j=0;j<gs->len_pw_list;j++)
	    {printf("%d ", *iptr); iptr++;}
	}
      printf("\n");
    
      
      printf("processor_list: ");
      if (iptr = gs->pair_list)
	{
	  for (j=0;j<gs->num_pairs;j++)
	    {printf("%d ", *iptr); iptr++;}
	}
      printf("\n");

      printf("size_list: ");
      if (iptr = gs->msg_sizes)
	{
	  for (j=0;j<gs->num_pairs;j++)
	    {printf("%d ", *iptr); iptr++;}
	}
      printf("\n");

      if (iptr = gs->pair_list)
	{
	  for (j=0;j<gs->num_pairs;j++)
	    {
	      printf("node_list %d: ", *iptr);
	      if (iptr2 = (gs->node_list)[j])
		{
		  for (k=0;k<(gs->msg_sizes)[j];k++)
		    {printf("%d ", *iptr2); iptr2++;}
		}
	      iptr++;
	      printf("\n");
	    }
	}
      printf("\n");
      
      printf("elm_list(U): ");
      if (iptr = gs->elms)
	{
	  for (j=0;j<gs->nel;j++)
	    {printf("%d ", *iptr); iptr++;}
	}
      printf("\n");
      printf("\n");
      
      printf("elm_list(T): ");
      if (iptr = gs->local_elms)
	{
	  for (j=0;j<gs->nel_total;j++)
	    {printf("%d ", *iptr); iptr++;}
	}
      printf("\n");
      printf("\n");
      
      printf("map_list(T): ");
      if (iptr = gs->companion)
	{
	  for (j=0;j<gs->nel;j++)
	    {printf("%d ", *iptr); iptr++;}
	}
      printf("\n");
      printf("\n");
      

      /* local exchange information */
      printf("\nLocal Info:\n");
      printf("local_strength=%d\n",   gs->local_strength);  
      printf("num_local_total=%d\n",  gs->num_local_total);
      printf("num_local=%d\n",        gs->num_local);
      printf("num_local_gop=%d\n",    gs->num_local_gop);
      printf("num_local_reduce=%d\n", gs->num_local_reduce);
      printf("local_reduce=%d\n",     gs->local_reduce);
      printf("num_gop_local_reduce=%d\n", gs->num_gop_local_reduce);
      printf("gop_local_reduce=%d\n",     gs->gop_local_reduce);
      printf("\n");

      for (j=0;j<gs->num_local;j++)
	{
	  printf("local reduce_list %d: ", j);
	  if (iptr2 = (gs->local_reduce)[j])
	    {
	      if ((gs->num_local_reduce)[j] <= 0)
		{printf("oops");}
	  
	      for (k=0;k<(gs->num_local_reduce)[j];k++)
		{printf("%d ", *iptr2); iptr2++;}
	    }
	  printf("\n");
	}
      
      printf("\n");
      printf("\n");
      
      for (j=0;j<gs->num_local_gop;j++)
	{
	  printf("gop reduce_list %d: ", j);
	  iptr2 = (gs->gop_local_reduce)[j];
		
	  if ((gs->num_gop_local_reduce)[j] <= 0)
	    {printf("oops");}
	  

	  for (k=0;k<(gs->num_gop_local_reduce)[j];k++)
	    {printf("%d ", *iptr2); iptr2++;}
	  printf("\n");
	}
      printf("\n");
      printf("\n");

      /* crystal router information */
      printf("\n\n");
      printf("Tree Info:\n");
      printf("max_left_over=%d\n",   gs->max_left_over);
      printf("num_in_list=%d\n",     gs->in_num);  
      printf("in_list=%d\n",         gs->in_list);  
      printf("num_out_list=%d\n",    gs->out_num);  
      printf("out_list=%d\n",        gs->out_list);  

      printf("\n\n");
    }
    }
  gsync();
  fflush(stdout);
}



/******************************************************************************
Function: gs_free()

Input : 

Output: 

Return: 

Description:  
  if (gs->sss) {perm_free((void*) gs->sss);}
******************************************************************************/
void
gs_free(register gs_id *gs)
{
  register int i;
  
  #ifdef INFO
  perm_stats();
  bss_stats();
  #endif
  

  if (!gs) 
    {
#ifdef WARNING      
      error_msg_warning("NULL ptr passed to gs_free()");
#endif
      return;
    }
      

  if (gs->ngh_buf) {bss_free((void*) gs->ngh_buf);}
  if (gs->elms) {bss_free((void*) gs->elms);}
  if (gs->local_elms) {bss_free((void*) gs->local_elms);}
  if (gs->companion) {bss_free((void*) gs->companion);}


  if (gs->vals) {perm_free((void*) gs->vals);}
  if (gs->neighbors) {perm_free((void*) gs->neighbors);}
  if (gs->pw_elm_list) {perm_free((void*) gs->pw_elm_list);}
  if (gs->pw_vals) {perm_free((void*) gs->pw_vals);}
  if (gs->msg_ids_in) {perm_free((void*) gs->msg_ids_in);}  
  if (gs->msg_ids_out) {perm_free((void*) gs->msg_ids_out);}
  if (gs->in) {perm_free((void*) gs->in);}
  if (gs->out) {perm_free((void*) gs->out);}

  /* pairwise info */
  for (i=0;i<gs->num_pairs;i++)
    {if (gs->node_list[i]) {perm_free((void*) gs->node_list[i]);}}
  if (gs->node_list) {perm_free((void*) gs->node_list);}
  if (gs->msg_sizes) {perm_free((void*) gs->msg_sizes);}
  if (gs->pair_list) {perm_free((void*) gs->pair_list);}

  /* local info */
  for (i=0;i<gs->num_local_total+1;i++)
    {
      if (gs->gop_local_reduce[i]) 
	{perm_free((void*) gs->gop_local_reduce[i]);}
    }
  if (gs->gop_local_reduce) {perm_free((void*) gs->gop_local_reduce);}
  if (gs->num_gop_local_reduce) {perm_free((void*) gs->num_gop_local_reduce);}

  /* LATER :: remove tree info */
  perm_free((void *) gs);

  #ifdef INFO
  perm_stats();
  bss_stats();
  #endif
}



/******************************************************************************
Function: gather_scatter

Input : 
Output: 
Return: 
Description: 

NOT FUNCTIONAL - FROM OLD VERSION!!!
***********************************************************************/
int 
gs_dump_ngh(gs_id *id, int loc_num, int *num, int *ngh_list)
{
  register int size, *ngh_buf;

  return(0);

  size = id->mask_sz;
  ngh_buf = id->ngh_buf;
  ngh_buf += (size*loc_num);

  size *= INT_LEN;
  *num = ct_bits((char *) ngh_buf, size);
  bm_to_proc((char *) ngh_buf, size, ngh_list);  
  /* ... */
}




/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
giop(int *vals, int *work, int n, int *oprs)
{
  int edge, type, dest, source, len;
  long msg_id;
  vfp fp;


  if (!(num_nodes>>1)) return;

  /* set function pointer to appropriate vector operation */
  type = oprs[0];
  if ((fp = (vfp) ivec_fct_addr(type)) == NULL)
    {error_msg_fatal("Only Not A Recognized giop() operation!");}
  
  /* advance to list of n operations for custom */
  if (type == NON_UNIFORM)
    {oprs = &oprs[1];}
  else
    {oprs = NULL;}

  /* all msgs will be of the same length */
  len = n*INT_LEN;

  /* if not a hypercube must colapse partial dim */
  if (edge_not_pow_2)
    {
      if (my_id >= floor_num_nodes)
	{
	  type = 10001 + my_id;
	  msg_id = isend(type,vals,len,edge_not_pow_2,0);
	  msgwait(msg_id);
	}
      else if (my_id < floor_num_nodes)
	{
	  type = 10001 + edge_not_pow_2;
	  msg_id = irecv(type,work,len);
	  msgwait(msg_id);
	  (*fp)(vals, work, n, oprs);
	}
    }

  /* implement the hypercube interleaved exchange algorithm */
  /* LATER: code and test the Interleaved algorithm */
  if (my_id<floor_num_nodes)
    {
      for (edge=0;edge<i_log2_num_nodes;edge++)
	{
	  source = dest = edge_node[edge];
	  type = 100001 + my_id + (num_nodes*edge);
	  msg_id = isend(type,vals,len,dest,0);
	  msgwait(msg_id);
	  /*	  msgignore(msg_id); */
	  /*LATER: which is faster? 2 msgwaits or msgmerge;msgwait */
	
	  type =  type - my_id + source;
	  msg_id = irecv(type,work,len);
	  msgwait(msg_id);
	  (*fp)(vals, work, n, oprs);
	}
    }

  /* if not a hypercube must expand to partial dim */
  if (edge_not_pow_2)
    {
      if (my_id >= floor_num_nodes)
	{
	  type = 1000001 + edge_not_pow_2;
	  msg_id = irecv(type,work,len);
	  msgwait(msg_id);
	  ivec_copy(vals, work, n);
	}
      else if (my_id < floor_num_nodes)
	{
	  type = 1000001 + my_id;
	  msg_id = isend(type,vals,len,edge_not_pow_2,0);
	  msgwait(msg_id);
	}
    }
}  


/******************************************************************************
Function: giop
Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
mpi_giop(int *vals, int *work, int n, int *oprs)
{
  int edge, type, dest, source, len;
  long msg_id;
  vfp fp;
#ifdef MPISRC
  MPI_Status status;
#endif


  if (num_nodes<2)
    {return;}


  /* set function pointer to appropriate vector operation */
  type = oprs[0];
  if ((fp = (vfp) ivec_fct_addr(type)) == NULL)
    {error_msg_fatal("Only Not A Recognized giop() operation!");}
  
  /* advance to list of n operations for custom */
  if (type == NON_UNIFORM)
    {oprs = &oprs[1];}
  else
    {oprs = NULL;}

  /* all msgs will be of the same length */
  len = n*INT_LEN;

  /* if not a hypercube must colapse partial dim */
  if (edge_not_pow_2)
    {
      if (my_id >= floor_num_nodes)
	{
	  type = 10001 + my_id;
#ifdef NXSRC
	  msg_id = isend(type,vals,len,edge_not_pow_2,0);
	  msgwait(msg_id);
#endif
#ifdef MPISRC
/* Use MPI_Send instead of MPI_Isend as there is no useful work inbetween */
/* Use n and MPI_INT instead of len and MPI_BYTE */
	  MPI_Send(vals, n, MPI_INT, edge_not_pow_2, type, MPI_COMM_WORLD);
#endif
	}
      else if (my_id < floor_num_nodes)
	{
	  type = 10001 + edge_not_pow_2;
#ifdef NXSRC
	  msg_id = irecv(type,work,len);
	  msgwait(msg_id);
#endif
#ifdef MPISRC
/* Use MPI_Recv instead of MPI_Irecv as there is no useful work inbetween */
/* Use n and MPI_INT instead of len and MPI_BYTE */
/* Could we use edge_not_pow_2 instead of MPI_ANY_SOURCE? */
	  MPI_Recv(work, n, MPI_INT, MPI_ANY_SOURCE, type, MPI_COMM_WORLD, 
		   &status);
#endif
/* Could have used hrecv() for this in NX? */ 
	  (*fp)(vals, work, n, oprs);
	}
    }

  /* implement the hypercube interleaved exchange algorithm */
  /* LATER: code and test the fan-in/out algorithm */
  if (my_id<floor_num_nodes)
    {
      for (edge=0;edge<i_log2_num_nodes;edge++)
	{
	  source = dest = edge_node[edge];
	  type = 100001 + my_id + (num_nodes*edge);
#ifdef NXSRC
	  msg_id = isend(type,vals,len,dest,0);
	  msgwait(msg_id);
#endif
#ifdef MPISRC
/* Use MPI_Send instead of MPI_Isend as there is no useful work inbetween */
/* Use n and MPI_INT instead of len and MPI_BYTE */
	  MPI_Send(vals, n, MPI_INT, dest, type, MPI_COMM_WORLD);
#endif
	  /*	  msgignore(msg_id); */
	  /*LATER: which is faster? 2 msgwaits or msgmerge;msgwait */
	
	  type =  type - my_id + source;
#ifdef NXSRC
	  msg_id = irecv(type,work,len);
	  msgwait(msg_id);
#endif
#ifdef MPISRC
/* Use MPI_Recv instead of MPI_Irecv as there is no useful work inbetween */
/* Use n and MPI_INT instead of len and MPI_BYTE */
/* Could we use something else instead of MPI_ANY_SOURCE? */
	  MPI_Recv(work, n, MPI_INT, MPI_ANY_SOURCE, type, MPI_COMM_WORLD, 
		   &status);
#endif
/* Could have used hrecv() for this in NX? */ 
	  (*fp)(vals, work, n, oprs);
	}
    }

  /* if not a hypercube must expand to partial dim */
  if (edge_not_pow_2)
    {
      if (my_id >= floor_num_nodes)
	{
	  type = 1000001 + edge_not_pow_2;
#ifdef NXSRC
	  msg_id = irecv(type,work,len);
	  msgwait(msg_id);
#endif
#ifdef MPISRC
/* Use MPI_Recv instead of MPI_Irecv as there is no useful work inbetween */
/* Use n and MPI_INT instead of len and MPI_BYTE */
/* Could we use something else instead of MPI_ANY_SOURCE? */
	  MPI_Recv(work, n, MPI_INT, MPI_ANY_SOURCE, type, MPI_COMM_WORLD, 
		   &status);
#endif
	  ivec_copy(vals, work, n);
	}
      else if (my_id < floor_num_nodes)
	{
	  type = 1000001 + my_id;
#ifdef NXSRC
	  msg_id = isend(type,vals,len,edge_not_pow_2,0);
	  msgwait(msg_id);
#endif
#ifdef MPISRC
/* Use MPI_Send instead of MPI_Isend as there is no useful work inbetween */
/* Use n and MPI_INT instead of len and MPI_BYTE */
	  MPI_Send(vals, n, MPI_INT, edge_not_pow_2, type, MPI_COMM_WORLD);
#endif
	}
    }
}  


/******************************************************************************
Function: gather_scatter

VERSION 3 :: 

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
void
gs_gop_pairwise(register gs_id *gs, register double *in_vals, char *operation)
{
  register double *dptr1, *dptr2, *dptr3, *in1, *in2;
  register int *iptr, *msg_list, *msg_size, **msg_nodes;
  register int *pw, *list, *size, **nodes;
  register int *msg_ids_in, *msg_ids_out, *ids_in, *ids_out;


  /* strip and load registers */
  msg_list =list         = gs->pair_list;
  msg_size =size         = gs->msg_sizes;
  msg_nodes=nodes        = gs->node_list;
  iptr=pw                = gs->pw_elm_list;  
  dptr1=dptr3            = gs->pw_vals;
  msg_ids_in  = ids_in   = gs->msg_ids_in;
  msg_ids_out = ids_out  = gs->msg_ids_out;
  dptr2                  = gs->out;
  in1=in2                = gs->in;

  /* load gs values into in out gs buffers */  
  while (*iptr >= 0)
    {*dptr3++ = *(in_vals + *iptr++);}

  /* load out buffers and post the sends */
  while (iptr = *msg_nodes++)
    {
      dptr3 = dptr2;
      while (*iptr >= 0)
	{*dptr2++ = *(dptr1 + *iptr++);}
      *msg_ids_out++ = isend(1001+my_id,dptr3,*(msg_size++)<<3,*msg_list++,0);
    }

  /* post the receives */
  do 
    {
      *msg_ids_in++ = irecv(1001 + *list++,in1,*size<<3);
      in1 += *size++;
    }
  while (*msg_ids_in >= 0);

  /* process the received data */
  while (iptr = *nodes++)
    {
      msgwait(*ids_in++);
      while (*iptr >= 0)
	{*(dptr1 + *iptr++) += *in2++;}
    }

  /* replace vals */
  while (*pw >= 0)
    {*(in_vals + *pw++) = *dptr1++;}

  /* clear isend message handles */
  while (*ids_out >= 0)
    {msgwait(*ids_out++);}

}
/*
new w/daxpy


new w/ivec_Add
test.log.8.8:      64       0     225b  2.9443E-03  ainv MAX time
test.log.8.8:      64       0     225b  2.1902E-03  XXt MAX time
test.log.8.8:      64       0     225b  2.1836E-03  XXt2 MAX time

original
test.log.8.8.1:      64       0     225b  2.9565E-03  ainv MAX time
test.log.8.8.1:      64       0     225b  2.2564E-03  XXt MAX time
test.log.8.8.1:      64       0     225b  2.2511E-03  XXt2 MAX time
*/
