/********************************adj_list.c************************************
SPARSE MATRIX BANDWIDTH REDUCTION PACKAGE: driver.c, smbwr.c,
                                           adj_list.c, and Makefile.

Author: Henry M. Tufo III

e-mail: hmt@@cs.brown.edu (no I didn't stutter on the <at> but lgrind ...)

sn-mail: Division of Applied Mathematics, Brown University,Providence, RI 02912

Last Modification: 5.2.95
*********************************adj_list.c***********************************/



/********************************adj_list.c************************************
NOTES ON USAGE: 

*********************************adj_list.c***********************************/



/* C modules for I/O etc. */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


/* mine : const before types! */
#include "const.h"
#include "types.h"
#include "error.h"
#include "bss_malloc.h"
#include "adj_list.h"


/*	
extern int N;
extern adj_node_ptr adj_list;


#define   PASS	       1
#define   FAIL	       0
*/


/********************************smbwr.c**************************************
Function reset_adj_list():

Input : none.
Output: neg. vals for adj list lengths.
Return: none
Description: not much
*********************************smbwr.c*************************************/
void
reset_adj_list(adj_node_ptr adj_list, int N)
{
int i;


for (i=0;i<N;i++)
  {
     if ((adj_list + i) -> elmt > 0)
       {(adj_list + i) -> elmt = -(adj_list + i) -> elmt;}
   }

}



/********************************smbwr.c**************************************
Function insert_adj_list():

Input : head number and edge to insert in adj list.
Output: inserted edge in adj list.
Return: error flag.
Description: 
*********************************smbwr.c*************************************/
int
insert_adj_list(adj_node_ptr adj_list, int vertice, int edge, double val)
{
adj_node_ptr new, tmp;



/* sorted insertion	                                   */
/* note negative value -> not included in CM construction  */
/*--  (adj_list+vertice)->elmt;*/
++  (adj_list+vertice)->elmt;
if ((new = (adj_node_ptr) bss_malloc(sizeof(adj_node))) != NULL)
  {
    /* search for position */
    tmp = (adj_list+vertice);
    while ((tmp->next!=NULL)&&(edge>(tmp->next)->elmt))
       {tmp=tmp->next;}

    /* insert  */
    /* set new block */
    new->elmt = edge;
    new->val  = val;
    new->next = tmp->next;
    tmp->next = new;

    /* ok we're done */
    return(PASS);
  }

error_msg_fatal("NOT ENOUGH SPACE TO BUILD ADJ LISTS !!!");
return(FAIL);
}



/********************************smbwr.c**************************************
Function insert_adj_list():

Input : head number and edge to insert in adj list.
Output: inserted edge in adj list.
Return: error flag.
Description: 
*********************************smbwr.c*************************************/
int
nc_insert_adj_list(adj_node_ptr adj_list, int vertice, int edge, double val)
{
adj_node_ptr new, tmp, empty;



/* sorted insertion	                                   */
/* note negative value -> not included in CM construction  */
/*--  (adj_list+vertice)->elmt;*/
++  (adj_list+vertice)->elmt;
if ((new = (adj_node_ptr) bss_malloc(sizeof(adj_node))) != NULL)
  {
    /* search for position */
    empty = tmp = (adj_list+vertice);
    

    while ((tmp->next!=NULL)&&(edge>=(tmp->next)->elmt))
       {tmp=tmp->next;}

    if (tmp==empty)
      {
	/* set new block */
	new->elmt = edge;
	new->val  = val;
	new->next = tmp->next;
	tmp->next = new;
	
      }
    else if (edge!=(tmp->elmt))
	{
	  /* set new block */
	  new->elmt = edge;
	  new->val  = val;
	  new->next = tmp->next;
	  tmp->next = new;
	}
    else
      {
	--  (adj_list+vertice)->elmt;
	bss_free(new);
      }
    
    /* ok we're done */
    return(PASS);
  }

error_msg_fatal("NOT ENOUGH SPACE TO BUILD ADJ LISTS !!!");
return(FAIL);
}



/********************************smbwr.c**************************************
Function free_adj_list():

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
void
free_adj_list(adj_node_ptr adj_list, int N)
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



/********************************smbwr.c**************************************
Function print_adj_list_hb():

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
void
print_adj_list_hb(adj_node_ptr adj_list, int N, int sda, int type)
{
adj_node_ptr ptr;
int i;
int k;
double val;
int nnz=0;
int tmp;
int totl, npl, nil, nvl;

/* return linked portion */
if (adj_list != NULL)
  {
    for (i=0;i<N;i++)
      {nnz+=(adj_list+i)->elmt;}

    npl = (N+1)/10;
    if ((N+1)%10)
      {npl++;}


    nil = nnz/10;
    if (nnz%10)
      {nil++;}

    nvl = nnz/3;
    if (nnz%3)
      {nvl++;}
    totl = npl+nil+nvl;

    
    printf("%14d",totl);
    printf("%14d",npl);
    printf("%14d",nil);
    printf("%14d",nvl);
    printf("%14d          \n",0);


    if (type == UT)
      {
	printf("%14s","UT            ");
	printf("%14d",N);
	printf("%14d",sda);
	printf("%14d",nnz);
	printf("%14d          \n",0);
      }
    else if (type == LT)
      {
	printf("%14s","LT            ");
	printf("%14d",N);
	printf("%14d",sda);
	printf("%14d",nnz);
	printf("%14d          \n",0);
      }
    else if (type == FULL)
      {
	printf("%14s","F             ");
	printf("%14d",N);
	printf("%14d",sda);
	printf("%14d",nnz);
	printf("%14d          \n",0);
      }
    else if (type == SYMM)
      {
	printf("%14s","S             ");
	printf("%14d",N);
	printf("%14d",sda);
	printf("%14d",nnz);
	printf("%14d          \n",0);
      }
    else
      {error_msg_fatal("unrecognized type");}


    printf("%s","(10I8)          (10I8)          (3E25.17)");
    printf("%s\n","                                       ");


    /* pointers */
    k=0;
    tmp=1;
    for (i=0;i<N;i++)
      {
	printf("%8d",tmp);
	tmp+=(adj_list+i)->elmt;
	k++;
	if (!(k%10))
	  {printf("\n");}
      } 
    printf("%8d",tmp);
    k++;
    if (k%10)
      {
	for(i=0;i<(10-(k%10));i++)
	  {printf("        ");}
      }
    printf("\n");

   
    /* indicies */
    k=0;
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	
	while (ptr != NULL)
	  {
	    printf("%8d",ptr->elmt+1);
	    ptr = ptr->next;
	    k++;
	    if (!(k%10))
	      {printf("\n");}
	  }
      }
    if (k%10)
      {
	for(i=0;i<(10-(k%10));i++)
	  {printf("        ");}
	printf("\n");
      }


    /* values */
    k=0;
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	
	while (ptr != NULL)
	  {
	    printf("%26.17E",ptr->val);
	    ptr = ptr->next;
	    k++;
	    if (!(k%3))
	      {printf("  \n");}
	  }
      }
    if (k%3)
      {
	for(i=0;i<(3-(k%3));i++)
	  {printf("                          ");}
	printf("  \n");
      }
  }
}




/********************************smbwr.c**************************************
Function init_adj_list():

Input : none.
Output: space.
Return: error flag.
Description: allocates space for array of linked list heads.
*********************************smbwr.c*************************************/
adj_node_ptr
init_adj_list(int N)
{
int i;
adj_node_ptr adj_list;


/* space - N verticies */
if ((adj_list = (adj_node_ptr) bss_malloc(N*sizeof(adj_node))) != NULL)
  {
    /* Init  - first node holds degree of vertex */
    for (i=0;i<N;i++)
      {(adj_list+i)->elmt=0;(adj_list+i)->val=0;(adj_list+i)->next = NULL;}
    return(adj_list);
  }

error_msg_fatal("NOT ENOUGH SPACE TO INIT ADJ LIST!!!");
return(NULL);
}



/********************************smbwr.c**************************************
Function print_adj_list():

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
void
print_adj_list_matlab(adj_node_ptr adj_list, int N)
{
adj_node_ptr ptr;
int i;
int k;
double val;


/* return linked portion */
if (adj_list != NULL)
  {
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	while (ptr != NULL)
	  {
	    /*	    printf("%8d%8d%25.16E\n",ptr->elmt+1,i+1,ptr->val); */
	    printf("%8d%8d%25.16E\n",i+1,ptr->elmt+1,ptr->val); 
	    ptr = ptr->next;
	  }
      }
  }
else
  {
    printf("%d not yet donut time!!!\n",N);
    fflush(stdout);
  }
}




/********************************smbwr.c**************************************
Function is_adj_list_symm()

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
int
is_adj_list_symm(adj_node_ptr adj_list, int N)
{
adj_node_ptr ptr;
int i;
int k;
double val;


/* return linked portion */
if (adj_list != NULL)
  {
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	while (ptr != NULL)
	  {
	    /*	    printf("%8d%8d%25.16E\n",ptr->elmt+1,i+1,ptr->val); */
	    if (!adj_find(adj_list,  N, ptr->elmt,i,ptr->val))
	      {return(FALSE);}
	    ptr = ptr->next;
	  }
      }
  }
else
  {
    error_msg_warning("Empty Matrix?!?");
    return(ERROR);
  }

return(TRUE);
}



/********************************smbwr.c**************************************
Function make_adj_list_symm()

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
int
make_adj_list_symm(adj_node_ptr adj_list, int N)
{
adj_node_ptr ptr;
int i;
int k;
double val;


/* return linked portion */
if (adj_list != NULL)
  {
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	while (ptr != NULL)
	  {
	    if (ptr->elmt > i)
	      {
		if (!adj_remove(adj_list,  N, ptr->elmt,i))
		  {return(FAIL);}
	      }
	    ptr = ptr->next;
	  }
      }
  }
else
  {
    error_msg_warning("Empty Matrix?!?");
    return(FAIL);
  }

return(PASS);
}



/********************************smbwr.c**************************************
Function adj_remove()

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
int
adj_remove(adj_node_ptr adj_list, int N, int i, int j)
{
  adj_node_ptr ptr, tmp;
  int found = FALSE;

  tmp = (adj_list+i);
  ptr = (adj_list+i)->next;
  while (ptr != NULL)
    {
      if (ptr->elmt==j)
	{
	  tmp->next = ptr->next;
	  ((adj_list+i)->elmt)--;
	  bss_free(ptr);
	  return(PASS);
	}
      tmp = tmp->next;
      ptr = ptr->next;
    }
  return(FAIL);
}



/********************************smbwr.c**************************************
Function adj_find()

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
int
adj_find(adj_node_ptr adj_list, int N, int i, int j, double val)
{
  adj_node_ptr ptr;
  int found = FALSE;

  
  ptr = (adj_list+i)->next;
  while (ptr != NULL)
    {
      if ((ptr->elmt==j)&&((ptr->val - val) < EPS))
	{return(TRUE);}
      ptr = ptr->next;
    }
  return(FALSE);
}



/********************************smbwr.c**************************************
Function position()

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
int
position(adj_node_ptr adj_list, int i, int j)
{
  adj_node_ptr ptr;
  int ct = 0;
  

  ptr = (adj_list+i)->next;
  while (ptr != NULL)
    {

      if (ptr->elmt==j)
	{return(ct);}
      ct++;       
      ptr = ptr->next;
    }
  return(-1);
}




/********************************smbwr.c**************************************
Function fprint_adj_list_hb_no_vals():

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
void
fprint_adj_list_hb_no_vals(FILE *ofp,adj_node_ptr adj_list, int N, int sda,
			   int type) 
{
adj_node_ptr ptr;
int i;
int k;
double val;
int nnz=0;
int tmp;
int totl, npl=0, nil=0, nvl=0;

/* return linked portion */
if (adj_list != NULL)
  {
    for (i=0;i<N;i++)
      {nnz+=(adj_list+i)->elmt;}

    npl = (N+1)/10;
    if ((N+1)%10)
      {npl++;}


    nil = nnz/10;
    if (nnz%10)
      {nil++;}

/*  #ifdef 0
    nvl = nnz/3;
    if (nnz%3)
      {nvl++;}
#endif    */
    totl = npl+nil+nvl;

    
    fprintf(ofp,"%14d",totl);
    fprintf(ofp,"%14d",npl);
    fprintf(ofp,"%14d",nil);
    fprintf(ofp,"%14d",nvl);
    fprintf(ofp,"%14d          \n",0);


    if (type == UT)
      {
	fprintf(ofp,"%14s","UT            ");
	fprintf(ofp,"%14d",N);
	fprintf(ofp,"%14d",sda);
	fprintf(ofp,"%14d",nnz);
	fprintf(ofp,"%14d          \n",0);
      }
    else if (type == LT)
      {
	fprintf(ofp,"%14s","LT            ");
	fprintf(ofp,"%14d",N);
	fprintf(ofp,"%14d",sda);
	fprintf(ofp,"%14d",nnz);
	fprintf(ofp,"%14d          \n",0);
      }
    else if (type == FULL)
      {
	fprintf(ofp,"%14s","F             ");
	fprintf(ofp,"%14d",N);
	fprintf(ofp,"%14d",sda);
	fprintf(ofp,"%14d",nnz);
	fprintf(ofp,"%14d          \n",0);
      }
    else if (type == SYMM)
      {
	fprintf(ofp,"%14s","S             ");
	fprintf(ofp,"%14d",N);
	fprintf(ofp,"%14d",sda);
	fprintf(ofp,"%14d",nnz);
	fprintf(ofp,"%14d          \n",0);
      }
    else
      {error_msg_fatal("unrecognized type");}


    fprintf(ofp,"%s","(10I8)          (10I8)          (3E25.17)");
    fprintf(ofp,"%s\n","                                       ");


    /* pointers */
    k=0;
    tmp=1;
    for (i=0;i<N;i++)
      {
	fprintf(ofp,"%8d",tmp);
	tmp+=(adj_list+i)->elmt;
	k++;
	if (!(k%10))
	  {fprintf(ofp,"\n");}
      } 
    fprintf(ofp,"%8d",tmp);
    k++;
    if (k%10)
      {
	for(i=0;i<(10-(k%10));i++)
	  {fprintf(ofp,"        ");}
      }
    fprintf(ofp,"\n");

   
    /* indicies */
    k=0;
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	
	while (ptr != NULL)
	  {
	    fprintf(ofp,"%8d",ptr->elmt+1);
	    ptr = ptr->next;
	    k++;
	    if (!(k%10))
	      {fprintf(ofp,"\n");}
	  }
      }
    if (k%10)
      {
	for(i=0;i<(10-(k%10));i++)
	  {fprintf(ofp,"        ");}
	fprintf(ofp,"\n");
      }


    /* NO VALUES!!! */
/*  #ifdef 0
    k=0;
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	
	while (ptr != NULL)
	  {
	    fprintf(ofp,"%26.17E",ptr->val);
	    ptr = ptr->next;
	    k++;
	    if (!(k%3))
	      {fprintf(ofp,"  \n");}
	  }
      }
    if (k%3)
      {
	for(i=0;i<(3-(k%3));i++)
	  {fprintf(ofp,"                          ");}
	fprintf(ofp,"  \n");
      }
#endif    */
  }
}




/********************************smbwr.c**************************************
Function fprint_adj_list_hb():

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
void
fprint_adj_list_hb(FILE *ofp,adj_node_ptr adj_list, int N, int sda, int type)
{
adj_node_ptr ptr;
int i;
int k;
double val;
int nnz=0;
int tmp;
int totl, npl, nil, nvl;

/* return linked portion */
if (adj_list != NULL)
  {
    for (i=0;i<N;i++)
      {nnz+=(adj_list+i)->elmt;}

    npl = (N+1)/10;
    if ((N+1)%10)
      {npl++;}


    nil = nnz/10;
    if (nnz%10)
      {nil++;}

    nvl = nnz/3;
    if (nnz%3)
      {nvl++;}
    totl = npl+nil+nvl;

    
    fprintf(ofp,"%14d",totl);
    fprintf(ofp,"%14d",npl);
    fprintf(ofp,"%14d",nil);
    fprintf(ofp,"%14d",nvl);
    fprintf(ofp,"%14d          \n",0);


    if (type == UT)
      {
	fprintf(ofp,"%14s","UT            ");
	fprintf(ofp,"%14d",N);
	fprintf(ofp,"%14d",sda);
	fprintf(ofp,"%14d",nnz);
	fprintf(ofp,"%14d          \n",0);
      }
    else if (type == LT)
      {
	fprintf(ofp,"%14s","LT            ");
	fprintf(ofp,"%14d",N);
	fprintf(ofp,"%14d",sda);
	fprintf(ofp,"%14d",nnz);
	fprintf(ofp,"%14d          \n",0);
      }
    else if (type == FULL)
      {
	fprintf(ofp,"%14s","F             ");
	fprintf(ofp,"%14d",N);
	fprintf(ofp,"%14d",sda);
	fprintf(ofp,"%14d",nnz);
	fprintf(ofp,"%14d          \n",0);
      }
    else if (type == SYMM)
      {
	fprintf(ofp,"%14s","S             ");
	fprintf(ofp,"%14d",N);
	fprintf(ofp,"%14d",sda);
	fprintf(ofp,"%14d",nnz);
	fprintf(ofp,"%14d          \n",0);
      }
    else
      {error_msg_fatal("unrecognized type");}


    fprintf(ofp,"%s","(10I8)          (10I8)          (3E25.17)");
    fprintf(ofp,"%s\n","                                       ");


    /* pointers */
    k=0;
    tmp=1;
    for (i=0;i<N;i++)
      {
	fprintf(ofp,"%8d",tmp);
	tmp+=(adj_list+i)->elmt;
	k++;
	if (!(k%10))
	  {fprintf(ofp,"\n");}
      } 
    fprintf(ofp,"%8d",tmp);
    k++;
    if (k%10)
      {
	for(i=0;i<(10-(k%10));i++)
	  {fprintf(ofp,"        ");}
      }
    fprintf(ofp,"\n");

   
    /* indicies */
    k=0;
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	
	while (ptr != NULL)
	  {
	    fprintf(ofp,"%8d",ptr->elmt+1);
	    ptr = ptr->next;
	    k++;
	    if (!(k%10))
	      {fprintf(ofp,"\n");}
	  }
      }
    if (k%10)
      {
	for(i=0;i<(10-(k%10));i++)
	  {fprintf(ofp,"        ");}
	fprintf(ofp,"\n");
      }


    /* values */
    k=0;
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	
	while (ptr != NULL)
	  {
	    fprintf(ofp,"%26.17E",ptr->val);
	    ptr = ptr->next;
	    k++;
	    if (!(k%3))
	      {fprintf(ofp,"  \n");}
	  }
      }
    if (k%3)
      {
	for(i=0;i<(3-(k%3));i++)
	  {fprintf(ofp,"                          ");}
	fprintf(ofp,"  \n");
      }
  }
}



/********************************smbwr.c**************************************
Function fprint_adj_list():

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
void
fprint_adj_list_matlab(FILE *ofp,adj_node_ptr adj_list, int N)
{
adj_node_ptr ptr;
int i;
int k;
double val;


/* return linked portion */
if (adj_list != NULL)
  {
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	while (ptr != NULL)
	  {
	    /*	    printf("%8d%8d%25.16E\n",ptr->elmt+1,i+1,ptr->val); */
	    if (ptr->val != 0.0)
	      {fprintf(ofp,"%8d%8d%25.16E\n",i+1,ptr->elmt+1,ptr->val);}
	    else
	      {fprintf(ofp,"%8d%8d\n",i+1,ptr->elmt+1);}
	    ptr = ptr->next;
	  }
      }
  }
else
  {
    printf("%d not yet donut time!!!\n",N);
    fflush(stdout);
  }
}




/********************************smbwr.c**************************************
Function nnz_adj_list():

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
int
nnz_adj_list(adj_node_ptr adj_list, int N)
{
adj_node_ptr ptr;
int i, sum=0;


/* return linked portion */
if (adj_list != NULL)
  {
    for (i=0;i<N;i++)
      {sum+= (adj_list+i)->elmt;}
    return(sum);
  }

/* return list - N verticies */
return(-1);
}




/********************************smbwr.c**************************************
Function sfprint_adj_list():

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
void
sfprint_adj_list_matlab(FILE *ofp,adj_node_ptr adj_list, int N)
{
adj_node_ptr ptr;
int i;
int k;
double val;


/* return linked portion */
if (adj_list != NULL)
  {
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	while (ptr != NULL)
	  {
	    fprintf(ofp,"%8d\n",ptr->elmt+1);
	    ptr = ptr->next;
	  }
      }
  }
else
  {
    printf("%d not yet donut time!!!\n",N);
    fflush(stdout);
  }
}



/********************************smbwr.c**************************************
Function sfprint_adj_list():

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
void
ssfprint_adj_list_matlab(FILE *ofp,adj_node_ptr adj_list, adj_node_ptr rhs,
			 adj_node_ptr sol, int N)
{
adj_node_ptr ptr;
int i;
int k, row;
double r, s;
double val;


/* return linked portion */
if (adj_list != NULL)
  {
    for (i=0;i<N;i++)
      {
	ptr = (adj_list+i)->next;
	while (ptr != NULL)
	  {
	    row = ptr->elmt;
	    /*
	    fprintf(ofp,"%8d\n",row+1);
	    */
	    r = ((rhs+row)->next)->val;
	    s = ((sol+row)->next)->val;
	    fprintf(ofp,"%8d%26.17E%26.17E\n",row+1,r,s);
	    ptr = ptr->next;
	  }
      }
  }
else
  {
    printf("%d not yet donut time!!!\n",N);
    fflush(stdout);
  }
}



/********************************smbwr.c**************************************
Function is_adj_list_symm()

Input : pointer to adj head list.
Output: none
Return: none. 
Description: returns space.
*********************************smbwr.c*************************************/
int *
extract_od(adj_node_ptr adj_list, int N)
{
  adj_node_ptr ptr;
  int i;
  int k;
  double val;
  int *pairs, *iptr, nnz;


  nnz = nnz_adj_list(adj_list, N) - N;
  pairs = iptr = (int *) bss_malloc((nnz+2)*INT_LEN);
  nnz>>=1;
  *iptr++ = N;
  *iptr++ = nnz;

  /* return linked portion */
  if (adj_list != NULL)
    {
      for (i=0;i<N;i++)
	{
	  ptr = (adj_list+i)->next;
	  while (ptr != NULL)
	    {
	      if (ptr->elmt>i)
		{*iptr++ = i+1; *iptr++ = ptr->elmt+1; nnz--;}
	      ptr = ptr->next;
	    }
	}
    }
  else
    {error_msg_fatal("extract_od() :: adj_list == NULL");}

  if (nnz)
    {error_msg_fatal("extract_od() :: nnz not zero :: %d\n",nnz);}

  return(pairs);
}
