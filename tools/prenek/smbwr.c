
/***************************************************************************
SPARSE MATRIX BANDWIDTH REDUCTION PACKAGE: 
driver.c, smbwr.c, and Makefile.

Author:
Henry M. Tufo III

e-mail: 
hmt@cs.brown.edu

snail-mail:
Division of Applied Mathematics
Brown University
Providence, RI 02912

Last Modification: 
04.06.96
***************************************************************************/



/******************************************************************************
Notes On Usage: sm_bandwidth_reduction(input ptr,ouput ptr, processing type)

smbwr.c contains the sparse matrix bandwidth reduction code. The only 
entry point is through the function sm_bandwidth_reduction() which takes as 
input a pointer to the array containing the sparse representation of the 
matrix (LDA,# of (i,j |i!=j) pairs,i_1,j_1,...,i_#,j_#), a pointer to the 
premutation array (LDA integers), and the type of seed to begin 
Cuthill-Mckee. If an error occurs in processing a 0 is returned. else if
connected the bandwidth of the permuted matrix, else the -1*bandwidth of 
the permuted matrix. For seed types the following work:
 1..LDA -> use that vertice as seed.
 0      -> use first vertice encountered with minimum degree.
-1      -> use pseudo-peripherial vertex.
-2      -> pass all LDA vertices and choose one that minimizes bandwidth.     

Note: if you pass disconnected matrices into sm_... be forewarned ... there
is no guarantee that it will work or that it will produce outstanding results.
***************************************************************************/



/* included C libraries */
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "bss_malloc.h"
#include "error.h"



/* defines used in this modeule  - note stack depth for non-recursive  qs */
#define   MAX_STACK    100
#define   MAX_LEN      80
#define   OPT          38
#define   PASS	       1
#define   ERROR	       0
#define   DC          -1
#define   TRUE	       1
#define   FALSE	       0
#define   IN	       1
#define   OUT	       0
#define   INT          sizeof(int)



/* macros */
#define   MAX(x,y)     (x>y) ? x : y



/* data structures used to implement adjacency list */
typedef struct node *adj_node_ptr;

typedef struct node{
  int elmt;
  adj_node_ptr next;
} adj_node;




/* Global Variables */
static adj_node_ptr adj_list = NULL;
static int N;
static int min_degree = INT_MAX, max_degree = -1, who_min = -1;



/******************************************************************************
Function sm_bandwidth_reduction():

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
int
sm_bandwidth_reduction(int *in_list,int *out_list,int type)
{
  /* function prototypes */
  static int load_sm(int *);
  static int bw(int *,int *,int);
  static int dc_graph(int *);
  static int ppv(int *);
  static int cm_sparse(int *,int);
  static void del_adj_list();
  static void error_msg(char *);


  /* vars */
  /* adj_node_ptr ptr; */
  int i, rt = PASS;
  int in_bw = -1, out_bw = -1, ppv_bw = -1;
  int min = INT_MAX, who_m = -1;


  who_min=-1;

  /* start processing */
  if (rt=load_sm(in_list))
    {
      in_bw = bw(in_list,out_list,IN);

      /* don't use std cm process if graph is disconnected   */
      if ((rt = dc_graph(out_list)) != DC)
	{
	  if ((who_min<0)||(who_min>=N))
	    {error_msg_fatal("sm_bandwidth_reduction() :: who_min=%d",who_min);}

	  if (type == 0)
            {rt = cm_sparse(out_list,who_min);}
	  else if ((type > 0) && (type <= N))
            {rt = cm_sparse(out_list,type - 1);}
	  else if (type == -1)
            {rt = cm_sparse(out_list,ppv(out_list));}
	  /*
	    else if (type == -3)
            {rt = cm_sparse_sp(out_list,who_min);}
	  */
	  else if (type == -2)
            {
	      type = 0;
	      while ((rt) && (type < N))
		{
		  if (rt = cm_sparse(out_list,type))
		    {
		      out_bw = bw(in_list,out_list,OUT);

		      /* keep track of the best */
		      if (out_bw < min)
			{min = out_bw; who_m = type;}
    
		      /* reset adj */     
		      for (i=0;i<N;i++)
			{(adj_list + i) -> elmt = -(adj_list + i) -> elmt;}

		      type++;
		    }
		}
	      /* dump the best */
	      if (rt)	
		{
		  rt = cm_sparse(out_list,ppv(out_list));
		  ppv_bw = bw(in_list,out_list,OUT);

		  /* reset adj */        
		  for (i=0;i<N;i++)
		    {(adj_list + i) -> elmt = -(adj_list + i) -> elmt;}

		  rt = cm_sparse(out_list,who_m);
		}
	    }
          else
	    {rt=ERROR;}
	}

      else
	{
	  error_msg_warning("DC matrix ...\n");
	  /*
	    if (type > 0)
            {rt = cm_sparse(out_list,type-1);}
	    else
            {
	    error_msg_warning("Program Supports Only Pos Types For DC Case!\n");
	    del_adj_list();
	    return(ERROR);
	    }
	  */
	  del_adj_list();
	  return(ERROR);
	}
    }

  /* return space */
  del_adj_list();

  /* close out */
  if (rt == PASS)
    {
      /*
	if (min != INT_MAX)
	{
	printf("MIN BW OVER ALL SEEDS :: %d \n",min);
	printf("PPV BW ...            :: %d \n",ppv_bw);
	}
      */

      out_bw = bw(in_list,out_list,OUT);

      /* how good was our processing? */
      /*      printf("IN  BW                :: %d \n",in_bw);
	      printf("OUT BW                :: %d \n",out_bw);
	      */

      for (i=0;i<N;i++)
        {*(out_list+i)+=1;}

      return(out_bw);
    }
  else if (rt == DC)
    {
      out_bw = bw(in_list,out_list,OUT);
   
      /* how good was our processing? */
      printf("IN  BW                :: %d \n",in_bw);
      printf("OUT BW                :: %d \n",out_bw);

      for (i=0;i<N;i++)
	{*(out_list+i)+=1;}

 
      return(-1*out_bw);
    }
  else
    {return(ERROR);}


  /* check adj and ...*/
  /* EXTRA DEBUG CODE */
  /*
    printf("MAX: %d \n",max_degree);
    printf("MIN: %d \n",min_degree);



    for (i=0;i<N;i++)
    {
    printf("%d: %d ",i, (adj_list+i)->elmt);
    ptr = (adj_list+i)->next;
    while (ptr != NULL)
    {
    printf("%d ", (ptr)->elmt);
    ptr = ptr->next;
    }
    printf("\n");
    }

  */
}



/******************************************************************************
Function load_sm(): created the adj linked list structure.

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
int
load_sm(int *in_list)
{
  /* protos */
  static int init_adj_list();
  static int insert_adj_list(int ,int );
  static void error_msg(char *);


/* vars */
  int i,j,k;
  int MNNZ;
  int rt;
  int count = 0, lda_i = 0, lda_j = 0, lda = 0;



/* dim matrix and number of non-zero entries excluding diags and symm pairs */
  N = *(in_list);
  MNNZ = 2 + 2 * *(in_list + 1);


/* get space and initialize adj lists - dependent on global param N */
  if (rt=init_adj_list())
    {
      /* load (i,j) pairs into adj lists   */
      /* note 1...N converted to 0...(N-1) */
      k=2;
      while ((k<MNNZ) && rt)
	{
          i = *(in_list+k++) - 1;
          j = *(in_list+k++) - 1;

          /* check to see if in bounds */
          if ((i<N) && (i>=0) && (j<N) && (j>=0))
	    {
	      /* remove (i,i) pairs - they should not be there ... */
	      if (i != j)
		{
		  rt = (insert_adj_list(i,j) && insert_adj_list(j,i));
		  /* check that lda given is lda found */
		  lda_i = MAX(lda_i,i);
		  lda_j = MAX(lda_j,j);
 
		  /* check NZ and given are same */
		  count++;
		}
	    }
          else
	    {rt = ERROR; error_msg_warning("i OR j OUT OF BOUNDS\n"); }
	}
    }


  if (rt)
    { 
      lda = MAX(lda_i,lda_j);

      if (++lda != N)
	{rt = ERROR; error_msg_warning("LDA MISMATCH\n"); }
      /*     if (++lda_i != N)
	     {error_msg_warning("LDA_i MISMATCH\n"); } */
      if (++lda_j != N)
	{error_msg_warning("LDA_j MISMATCH\n"); }
      if (count != *(in_list + 1))
	{rt = ERROR; error_msg_warning("NZ MISMATCH\n"); }
    }


  /* ok the input looks good and no mem. alloc. pblms */
  if (rt)
    {
      min_degree = INT_MAX;
      max_degree = -1;

      /* seems ok so get min degree and first min and max */
      for (i=0;i<N;i++)
	{
          if (abs((adj_list+i)->elmt) < min_degree)
            {min_degree = abs((adj_list+i)->elmt); who_min = i;}
  
          if (abs((adj_list+i)->elmt) > max_degree)
            {max_degree = abs((adj_list+i)->elmt);}
        }
     
      /* ok return -> pass */
      return(PASS);
    }

  /* oops */
  return(ERROR);
}



/******************************************************************************
Function init_adj_list():

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
int
init_adj_list()
{
  static void error_msg(char *error);

  int i;


  /* space - N verticies */
  if ((adj_list = (adj_node_ptr) bss_malloc(N*sizeof(adj_node))) != NULL)
    {
      /* Init  - first node holds degree of vertex */
      for (i=0;i<N;i++)
	{(adj_list+i)->elmt = 0; (adj_list+i)->next = NULL;}
      return(PASS);
    }

  error_msg_warning("NOT ENOUGH SPACE TO INIT ADJ LIST!!!\n");

  return(ERROR);
}



/******************************************************************************
Function insert_adj_list():

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
int
insert_adj_list(int vertice, int edge)
{
  static void error_msg(char *error);


  adj_node_ptr new;


  /* non-sorted insertion                                   */
  /* note negative value -> not included in CM construction */
  --  (adj_list+vertice)->elmt;
  if ((new = (adj_node_ptr) bss_malloc(sizeof(adj_node))) != NULL)
    {
      /* insert at head of list */
      new->elmt = edge;
      new->next = (adj_list+vertice)->next;
      (adj_list+vertice)->next = new;
  
      return(PASS);
    }

  error_msg_warning("NOT ENOUGH SPACE TO BUILD ADJ LISTS !!!\n");
  return(ERROR);
}



/******************************************************************************
Function del_adj_list():

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
void
del_adj_list()
{
  adj_node_ptr ptr, fr;
  int i;


  /* return linked portion */
  if (adj_list != NULL)
    {
      for (i=0;i<N;i++)
	{
          fr = ptr = (adj_list+i)->next;
          while (ptr != NULL)
            {
	      ptr = ptr->next;
	      bss_free(fr);
	      fr = ptr;
	    }
	}
    }

  /* return list - N verticies */
  bss_free(adj_list);
}



/******************************************************************************
Function dc_graph():

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
int
dc_graph(int *out)
{
  int i,j;
  adj_node_ptr ptr;
  int where=1;
  int num_adj_list;
  int rt = PASS;



  /* initialize space -> hit -1 then graph is disconnected */
  for (i=1;i<N;i++)
    {*(out+i)=-1;}


  /* begin search with seed 0 */
  *(out)=0;
  (adj_list + *(out))->elmt = abs((adj_list + *(out))->elmt);

  i=0;
  while ((rt!=DC) && (i<N))
    {
      if (*(out +i) < 0)
	{return(DC);} 
	/*	{rt = DC;} */
      else
	{
          ptr = (adj_list + *(out +i));
          num_adj_list = ptr->elmt;
  
          for (j=0;j<num_adj_list;j++)
            {
	      ptr = ptr->next;
	      if ((adj_list + (ptr->elmt))->elmt < 0)
		{     
		  (adj_list + (ptr->elmt))->elmt =
		    -(adj_list + (ptr->elmt))->elmt;
		  *(out+where++) = ptr->elmt;
		}
	    }
          i++;
	}
    }

  /* fix adj list */
  for (i=0;i<N;i++)
    {
      if ((adj_list + i)->elmt > 0)
	{(adj_list + i)->elmt = -(adj_list + i)->elmt;}
    }

  return(rt);
}



/******************************************************************************
Function ppv():

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
int 
ppv(int *out)
{
  int i,j;
  adj_node_ptr ptr;
  int *sr;
  int head, tail;
  /* int where=1; */
  int num_adj_list;
  /* int rt = PASS; */
  int l,u;
  int level = -1, size = 0;
  int max_l = -1, max_who = -1;


  /* search space */
  sr = (int *) bss_malloc(N*INT);

  /* initialize space */
  for (i=0;i<N;i++)
    {*(sr+i) = *(out+i) = -1;}

  /* set search space to min degree node */
  head = tail = 0;
  *(sr) = who_min;

  /* build level structure for target */
  while (head <= tail)
    {

      /* printf("ppv:: head %d tail %d max_who %d \n",head,tail,max_who); */

      /* init list */
      for (j=0;j<N;j++)
	{
          *(out+j) = -1; 
          if ((adj_list +j)->elmt > 0)
            {(adj_list +j)->elmt = -((adj_list +j)->elmt);}
	}

      *(out)= *(sr+head);
      /*     printf("NEXT :: %d \n",*(out)= *(sr+head)); */
      (adj_list + *(out))->elmt = abs((adj_list + *(out))->elmt);

      l=u=1;
      level=0;
      i=0;
      head++;

      /* grab first level */   
      while(u<N)
	{
          l=u;
          size=0;
    
          while (i<l)
	    {
	      ptr = (adj_list + *(out +i));
	      num_adj_list = ptr->elmt;

	      for (j=0;j<num_adj_list;j++)
		{
		  ptr = ptr->next;
		  if ((adj_list + (ptr->elmt))->elmt < 0)
		    {     
		      (adj_list + (ptr->elmt))->elmt =
			-(adj_list + (ptr->elmt))->elmt;
		      *(out+u) = ptr->elmt;
		      size++; u++;
		    }
		}
	      i++;
	    }
  
	  level++;
	}


      if (level > max_l)
	{
	  max_l = level;
	  max_who = *(sr + head -1);
	  for (j=0;j<size;j++)
	    {
	      if ((adj_list+*(out))->elmt< (min_degree+2))
		{
		  tail++;
		  *(sr+tail) = *(out+N-1-j);
		}
	    }
	}

      if ((level == max_l) && ((adj_list+max_who)->elmt > (adj_list + *(out))->elmt))
	{max_who = *(out);}
    }

  /* fix adj list */
  for (i=0;i<N;i++)
    {
      if ((adj_list + i)->elmt > 0)
	{(adj_list + i)->elmt = -(adj_list + i)->elmt;}
    }


  bss_free(sr);

  /* printf("\n PPV :: %d \n",max_who);  */
  return(max_who);
}



/******************************************************************************
Function cm_sparse():

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
static
int
cm_sparse(int *out,int seed)
{
  static int sort(int *,int *, int);


  int i,j;
  adj_node_ptr ptr;
  int *degree, *nodes;
  int num_sort=0;
  int where=1;
  int num_adj_list;
  int tmp;
  int rt = PASS;



  /* get space to sort by degree */
  /*
    degree = (int *) bss_malloc(max_degree*INT);
    nodes  = (int *) bss_malloc(max_degree*INT);
  */
  degree = (int *) bss_malloc(N*INT);
  nodes  = (int *) bss_malloc(N*INT);

  /* initialize space -> hit -1 then graph is disconnected */
  for (i=0;i<N;i++)
    {*(out+i)=-1;}


  /* begin CM with seed */
  *(out)=seed;
  if ((adj_list + seed)->elmt != 0)
    {(adj_list + seed)->elmt = abs((adj_list + seed)->elmt);}
  else
    {(adj_list + seed)->elmt = INT_MAX;}


  /* start modified CM */
  i=0;
  while ((rt) && (i < N))
    {
      if (*(out +i) < 0)
	{
          rt = -1;
          num_sort = INT_MAX;
          num_adj_list = -1;
          for (j=N-1;j>=0;j--)
	    {
	      if ((tmp = (adj_list + j)->elmt) <= 0)
		{
		  if (abs(tmp) < num_sort)
		    {
		      num_sort = abs(tmp);
		      num_adj_list = j;
		    }
		}
	    }
          if (num_adj_list >= 0)
	    {
	      *(out +i) = num_adj_list;
	      if ((adj_list+num_adj_list)->elmt != 0)
		{(adj_list+num_adj_list)->elmt=
		   abs((adj_list + num_adj_list) -> elmt);}
	      else
		{(adj_list+num_adj_list)->elmt = INT_MAX;}
	      /*              printf("\nDISCONNECTED BUT PUSHED ON!!!\n"); */
	    }
          else
	    {
	      printf("\nDISCONNECTED BUT NO ONE TO GO TO!?!\n");
	      return(ERROR);
	    }
          where ++;
	}

      num_sort = 0;
      ptr = (adj_list + *(out +i));
      if (ptr->elmt != INT_MAX)
	{num_adj_list = ptr->elmt;}
      else
	{num_adj_list = 0;}

      for (j=0;j<num_adj_list;j++)
	{
          ptr = ptr->next;
          if ((adj_list + (ptr->elmt))->elmt < 0)
	    {
	      *(degree+num_sort) =  (adj_list + (ptr->elmt))->elmt
		= -(adj_list + (ptr->elmt))->elmt;
	      *(nodes+num_sort)  = ptr->elmt;
	      num_sort++;
	    }
	}

      /* sort in order of increasing degree */
      if (rt=sort(degree,nodes,num_sort))
	{
          /* insert into permutation array and index */
          for (j=0;j<num_sort;j++)
	    {*(out+where+j) = *(nodes + j);}
          where += num_sort;
	}
      i++;
    }

  bss_free(degree);
  bss_free(nodes);
  return(rt);
}


/******************************************************************************
Function bw()

Input : addresses of input space and ouput space.
Output: none.
Return: bandwidth
Description: computes bandwidth of either sparse matrix before or after 
permutation depending on input flag.
******************************************************************************/
static
int
bw(int *in_list,int *out_list,int flag)
{
  int i,j,k,tmp = 0;
  int *temp;
  int MNNZ;



  /* number to scan */
  MNNZ = 2 + 2 * *(in_list +1);


  /* unpermuted */
  if (flag)
    {
      for (k=2;k<MNNZ;)
	{
          i = *(in_list+k++) - 1;
          j = *(in_list+k++) - 1;
          if (abs(i-j) > tmp)
            {tmp = abs(i-j);}
	}
    }

  /* permuted -> actually construct inverse map */
  else
    {
      temp = bss_malloc(N*INT);

      for(i=0;i<N;i++)
	{*(temp + *(out_list + i)) = i;}

      for (k=2;k<MNNZ;)
	{
          i = *(temp + *(in_list+k++) - 1);
          j = *(temp + *(in_list+k++) - 1);
          if (abs(i-j) > tmp)
            {tmp = abs(i-j);}
	}
      bss_free(temp);
    }


  return(tmp);
}



/******************************************************************************
Function error_msg():

Input : pointer to message.
Output: prints message to stdio.
Return: none.
Description: print error message ...
******************************************************************************/
static 
void 
error_msg(char *error)
{
  printf("\n%s\n",error);
}



/******************************************************************************
Function sort():

Input : offset of list to be sorted, offset of companion, size of list.
Output: sorted list and companion.
Return: error value.
Description: This function used a combination of non-recursive quicksort 
w/pivot on median of three and bubble sort to sort a list of numbers in 
ascending order. The companion list passed in through the second argument
is permuted in the same fashion as the sorted list. The only possible sorce 
of error is exhaustion of stack space.
******************************************************************************/
static
int
sort(register int *array,register int *vert, int size)
{
  register int up, bottom, down, temp;
  int mid, pivot_value, stack_size=2;
  int stack[MAX_STACK+4];



  /* init. stack */
  stack[0]=0;
  stack[1]=size;


  /* process until stack is empty or overflow seems iminent */
  while ((stack_size) && (stack_size < MAX_STACK))
    {
      /* pop next off stack */
      size=*(stack+ --stack_size);
      bottom=*(stack+ --stack_size);
  

      /* bottom out on bubble */
      if (size<OPT) 
	{
          mid = TRUE;
          pivot_value = size+bottom; 
          for (up = 0; mid && ++up < size;)
	    {
	      mid = FALSE;
	      for (down = bottom; down < pivot_value - up; down++) 
		{
		  if (*(array+down) > *(array+down+1)) 
		    {
		      temp = *(array+down);
		      *(array+down) = *(array+down+1);
		      *(array+down+1) = temp;
		      temp = *(vert+down);
		      *(vert+down) = *(vert+down+1);
		      *(vert+down+1) = temp;
		      mid = TRUE;
		    }
		}
	    }
	}

      else
	{

          /* else another time through qs             */
          /* start at bottom and top of array segment */
          down=bottom;
          up=size-1+bottom;
          mid=size/2+bottom;


          /* order first last */
          if (*(array+bottom) > *(array+up))
	    {
	      temp = *(array+bottom);
	      *(array+bottom) = *(array+up);
	      *(array+up) = temp;
	      temp = *(vert+bottom);
	      *(vert+bottom) = *(vert+up);
	      *(vert+up) = temp;
	    }

          /* up is middle */
          if (*(array+mid) > *(array+up))
	    {
	      temp = *(array+bottom);
	      *(array+bottom) = *(array+up);
	      *(array+up) = temp;
	      temp = *(vert+bottom);
	      *(vert+bottom) = *(vert+up);
	      *(vert+up) = temp;
	    }

          /* mid is middle */
          else if (*(array+mid) > *(array+bottom))
	    {
	      temp = *(array+mid);
	      *(array+mid) = *(array+bottom);
	      *(array+bottom) = temp;
	      temp = *(vert+mid);
	      *(vert+mid) = *(vert+bottom);
	      *(vert+bottom) = temp;
	    }

	  /* else we're already done */


	  /* pivot value is now located in first position */
	  pivot_value=*(array+bottom);

	  /* partition: find place and divide */
	  while (down<up)
	    {
	      /* skip over values less than or equal to pivot */
	      while ((down<up) && (*(array+down)<=pivot_value)) 
		{down++;}          

	      /* skip over values greater than pivot */
	      while (*(array+up)>pivot_value)
		{up--;}

	      /* swap -> divide */
	      if (down<up)
		{
		  temp=*(array+down);
		  *(array+down)=*(array+up);
		  *(array+up)=temp;
		  temp=*(vert+down);
		  *(vert+down)=*(vert+up);
		  *(vert+up)=temp;
		}
	    }



	  /* place pivot value */
	  *(array+bottom)=*(array+up);
	  *(array+up)=pivot_value;
          temp = *(vert+bottom);
          *(vert+bottom) = *(vert+up);
          *(vert+up) = temp;



	  /* nonrecursive - put on stack */
          temp=up-bottom; 
          *(stack+ stack_size++)=bottom;
	  *(stack+ stack_size++)=temp++; 
          *(stack+ stack_size++)=down;
          *(stack+ stack_size++)=size - temp;

	}
    }

  /* check reason for termination */
  if (stack_size)
    {return(ERROR);}

  return(PASS);
}












