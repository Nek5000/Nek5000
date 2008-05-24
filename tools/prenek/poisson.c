/***************************************************************************
Henry M. Tufo III
AM258
Professor Fischer
HW#2
03.16.93
***************************************************************************/



/***************************************************************************
Develop conjugate gradient solver to use on two dim poisson problem:
***************************************************************************/



/* std c modules */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



/* set definitions for this module */
#define SIZE           sizeof(double)
#define LOG2	       log10(2)
#define CG_LIMIT       1000000 /* pow(n,6.0) */
#define RES_LIMIT      pow(2.0,-50.0)
#define ERROR_LIMIT    pow(10.0,-3.0)
#define PI             acos(-1.0)



/* c defined function proto. Definition in respective .h file */
FILE *fopen(char *name,char *mode);
int fclose(FILE *fp);
double pow(double x,double y);
int atoi(const char *s);
int fprintf(FILE * stream,const char *format, ...);
int printf(const char *format, ...);
double fabs(double x);
double sin(double x);
double acos(double x);
double sqrt(double x);
double atof(const char *s);
double log10(double x);

/* for some reason CFM gcc having problems with these 
void *   calloc(size_t nmemb, size_t size );
void free(void *p);
void  exit(int status);
*/



/* defined function proto from this module */
void init_files(int ac,char **av);
void load_data();
void clean_up();
void compute_b();
void compute_A();
void conjugate_grad(double *A,double *b);
double nodal_error(double * x);
double v_mult(double *x1,double *x2);
double * v_add(double *x1,double *x2);
double * v_const_mult(double x1,double *x2);
double * matrix_mult(double *A,double *b);
double * Init_guess();


/* global vars */
int  n;
FILE *ofp;
double a;
double *A, *soln, *b;



/* main function call */
main(int argc,char *argv[])
{
init_files(argc, argv);
load_data();
compute_A();
compute_b();
conjugate_grad(A,b); 
clean_up();
exit(0);
}

     
/*************************************************************************** 
Function init_files: sets ifp and ofp to uses input files. error checking
done on format and file creation.
***************************************************************************/
void init_files(int ac,char **av)
{

if (ac != 4) 
   {
      printf("\n Please provide the file name. 
              Format: PE <N> <a> <outputfilename> \n");
      exit(1);
   }

else
   {
      n = atoi(*++av);
      a  = atof(*++av);
      if ((ofp = fopen(*++av,"w")) == NULL)
         {printf("\n Can't open output file!!! \n"); exit(1);}
   }


/* file id */
/*
fprintf(ofp,"Henry M. Tufo III \n");
fprintf(ofp,"AM258 \n");
fprintf(ofp,"Professor Fischer \n");
fprintf(ofp,"HW#2 \n");
fprintf(ofp,"03.16.93 \n\n\n");
fprintf(ofp,"Coef: %.5f\n",a);
fprintf(ofp,"N_x: %d, N_y: %d \n\n",n,n);
*/
}



/*************************************************************************** 
Function load_data: allocates space for the program:  A,b 
***************************************************************************/
void load_data()
{

/* allocate space for Stiffness matrix A */
if ((A = (double *) calloc(5*n*n,SIZE)) == NULL)
  {fprintf(ofp,"Could not allocate space for A!?!\n"); exit(1);}

/* allocate space for vector resultant from mass matrixb */
if ((b = (double *) calloc(n*n,SIZE)) == NULL)
  {fprintf(ofp,"Could not allocate space for b!?!\n"); exit(1);}
}



/*************************************************************************** 
Function compute_A: given nodal points and h_i's the stiffness matrix A is
loaded. Recall A_ij = a(phi_i',phi_j') + (phi_i,phi_j).
***************************************************************************/
void compute_A()
{
int i,j,k;
double diag, rt, kitty;


/* first through n-1th row with look forward sym */
for (i=0;i<(n*(n-1));i++)
   {
     *(A+5*i)   = diag  = (a*8.0/3.0 + 4.0/(9.0*n*n));	
     *(A+5*i+1) = rt    = (a* -1.0/3.0 + 1.0/(9.0*n*n));
     *(A+5*i+2) = kitty = (a* -1.0/3.0 + 1.0/(36.0*n*n));
     *(A+5*i+3) = rt;
     *(A+5*i+4) = kitty;
   }

/* last row except last nodal point */
for (i=(n*(n-1));i<(n*n-1);i++)
   {
     *(A+5*i)   = diag/2.0;
     *(A+5*i+1) = rt/2.0;
   }

/* last nodal point */
*(A+5*n*n-5) = diag/4.0;

/* first nodal point in row has no left top neighbor */
for (i=0;i<n-1;i++)
   {*(A+5*i*n+2) = 0;}

/* nodal points in last column have no right neighbors and 1/2 area */
for (i=0;i<n-1;i++)
   {
      *(A+5*n*(i+1)-1) = 0;
      *(A+5*n*(i+1)-2) /= 2.0;
      *(A+5*n*(i+1)-4) = 0;
      *(A+5*n*(i+1)-5) /= 2.0;
    }
}



/*************************************************************************** 
Function compute_b(): This function is dependent upon f. b_i = (f,phi_i) 
where f = sin(PI*x/2)sin(PI*y/2) and phi(x,y)=phi(x)phi(y) , where phi(x) 
and phi(y) are one dim hat fcts.
***************************************************************************/
void compute_b()
{
int i,j;
double * temp;

/* Load all but last row and column elements */
for (i=0;i<(n-1);i++)
   for (j=0;j<(n-1);j++)
     {
       *(b+i*n+j) =16*n*n/pow(PI,4.0);
       *(b+i*n+j)*=(sin(PI*(j+2)/(2*n))-2*sin(PI*(j+1)/(2*n))+sin(PI*j/(2*n)));
       *(b+i*n+j)*=(sin(PI*(i+2)/(2*n))-2*sin(PI*(i+1)/(2*n))+sin(PI*i/(2*n)));
      }

/* Load the last column minus elm which is in last row and col */
j = (n-1);
for (i=0;i<(n-1);i++)
  {
     *(b+i*n+j)  = -16*n*n/pow(PI,4.0);
     *(b+i*n+j) *= (1 - sin(PI*j/(2*n)));
     *(b+i*n+j) *= sin(PI*(i+2)/(2*n))-2*sin(PI*(i+1)/(2*n))+sin(PI*i/(2*n));
   }

/* Load the last row minus elm which is in last row and col */
i = (n-1);
for (j=0;j<(n-1);j++)
  {
     *(b+i*n+j)  = -16*n*n/pow(PI,4.0); 
     *(b+i*n+j) *= (1 - sin(PI*i/(2*n)));
     *(b+i*n+j) *= (sin(PI*(j+2)/(2*n))-2*sin(PI*(j+1)/(2*n))+sin(PI*j/(2*n)));
   }

/* Last row last column element */
i = j = (n-1);
*(b+i*n+j)  = 16*n*n/pow(PI,4.0); 
*(b+i*n+j) *= (1 - sin(PI*i/(2*n)));
*(b+i*n+j) *= (1 - sin(PI*j/(2*n)));

}



/*************************************************************************** 
Function conjugate_grad(): solves Ax=b for A,b given and sparse representation 
for A. initial x_0 is the zero vector.
***************************************************************************/
void conjugate_grad(double *A,double *z)
{
int k;
double *r_k, *r_k1, *r_k2, *p_k, *p_k1, *x_k, *x_k1, *tpmm, *tpcm, *tp;
double alpha, beta, temp;
double *Init_guess();


/* initialize */
k = 0;
x_k = x_k1 = Init_guess();
r_k = z;

/*repeat until either desired error limit attained or we have solved sys s.t
the residual is essentially zero or have surpassed iteration limit */
while (/*(nodal_error(x_k)>ERROR_LIMIT) &&*/ (k<CG_LIMIT) 
                                      && (v_mult(r_k,r_k) > RES_LIMIT))
  {
     k++;
     if (k==1)
       {r_k1 = r_k2 = p_k = p_k1 = r_k;}
     else
       {
          tp = r_k2;
          r_k2 = r_k1;
          r_k1 = r_k;
          if (k>3) {free(tp);}
          tp = p_k1;
          p_k1=p_k;
          if (k>3) {free(tp);}
          tp=x_k1;
          x_k1=x_k;
          if (k>2) {free(tp);}
          beta = v_mult(r_k1,r_k1)/v_mult(r_k2,r_k2);
          tpcm = v_const_mult(beta,p_k1);
          p_k  = v_add(r_k1,tpcm);
          free(tpcm);
	}
     tpmm=matrix_mult(A,p_k);
     alpha = v_mult(r_k1,r_k1)/v_mult(p_k,tpmm);
     tpcm=v_const_mult(alpha,p_k);
     x_k = v_add(x_k1,tpcm);
     free(tpcm);
     tpcm=v_const_mult(-1.0*alpha,tpmm);
     r_k=v_add(r_k1,tpcm);
     free(tpcm);
     free(tpmm);

/* error, iteration count and problem size */
fprintf(ofp,"%.15f %.15f \n",log10(k),log10(nodal_error(x_k)));
     }

free(r_k);
free(r_k1);
free(r_k2);
free(p_k);
free(p_k1);
free(x_k1);

/* last iteration yields solution */
soln = x_k;     
}



/*************************************************************************** 
Function nodal_error(): this computes ||u_exact--u_soln||_infinity
***************************************************************************/
double nodal_error(double *y)
{
int i,j;
double error=0, err2=0;
double temp;

/*
for (i=0;i<n;i++)
  {
     for (j=0;j<n;j++)
        {fprintf(ofp,"%.3f, ",(2.0/(PI*PI*a+2)*
                    (sin(PI*(j+1)/(2*n))*sin(PI*(i+1)/(2*n)))));}
     fprintf(ofp,"\n");
   }


     fprintf(ofp,"\n");

for (i=0;i<n;i++)
  {
     for (j=0;j<n;j++)
        {fprintf(ofp,"%.3f, ",*(y+i*n+j));}
     fprintf(ofp,"\n");
   }


     fprintf(ofp,"\n");     
     fprintf(ofp,"\n");
*/


for (i=0;i<n;i++)
  for (j=0;j<n;j++)
    {
       temp=fabs(*(y+i*n+j)
                   -(2.0/(a*PI*PI+2)*sin(PI*(j+1)/(2*n))*sin(PI*(i+1)/(2*n))));
       error = (error > temp) ? error : temp;
     }

for (i=0;i<n;i++)
  for (j=0;j<n;j++)
    {
       temp=fabs(2.0/(a*PI*PI+2)*sin(PI*(j+1)/(2*n))*sin(PI*(i+1)/(2*n)));
       err2 = (err2 > temp) ? err2 : temp;
     }

return(error/err2);
}



/*************************************************************************** 
Function clean_up(): released space and closes open files.
***************************************************************************/
void clean_up()
{

free(A);
free(b);
free(soln);

if (fclose(ofp)!= 0)
  {printf("\n Problem closing output file!!! \n");}
}



/*************************************************************************** 
Function Init_guess(): seed x_0 for CG. Take it to be zero.
***************************************************************************/
double *Init_guess()
{
int i,j;
double *temp;


/* allocate space for vector resultant from mass matrixb */
if ((temp = (double *) calloc(n*n,SIZE)) == NULL)
  {fprintf(ofp,"Could not allocate space for Init_guess temp!?!\n"); exit(1);}

for (i=0;i<n;i++)
   for (j=0;j<n;j++)
      {*(temp+i*n+j) = 0.0;} 

return(temp);
}



/*************************************************************************** 
Function v_mult(): vector vector dot product: <x1,x2> = x1^t.x2
***************************************************************************/
double  v_mult(double *x1,double *x2)
{
int i;
double temp;

temp=0;
for (i=0;i<(n*n);i++)
  {temp+= *(x1+i)* *(x2+i);}

return(temp);
                                                                 
}




/*************************************************************************** 
Function v_add(): adds two vectors: x1+x2
***************************************************************************/
double *v_add(double *x1,double *x2)
{
int i;
double *temp;


/* allocate space for vector resultant from mass matrixb */
if ((temp = (double *) calloc(n*n,SIZE)) == NULL)
  {fprintf(ofp,"Could not allocate space for v_add temp!?!\n"); exit(1);}


for (i=0;i<(n*n);i++)
  {*(temp + i) = *(x1+i) + *(x2+i);}

return(temp);             
}



/*************************************************************************** 
Function v_const_mult(): multiplies a vector by a real constant: alph*ax2_i.
***************************************************************************/
double *v_const_mult(double alpha,double *x2)
{
int i;
double *temp;


/* allocate space for vector resultant from mass matrix b */
if ((temp = (double *) calloc(n*n,SIZE)) == NULL)
  {fprintf(ofp,"Could not allocate space for v_add temp!?!\n"); exit(1);}


for (i=0;i<(n*n);i++)
  {*(temp + i) = alpha * *(x2+i);}

return(temp);             
}




/*************************************************************************** 
Function matrix_mult(): matrix vector product. Sparse matrix representation.
***************************************************************************/
double *matrix_mult(double *A,double *x)
{
int i,j;
double *temp;

/* allocate space for vector resultant from mass matrixb */
if ((temp = (double *) calloc(n*n,SIZE)) == NULL)
  {fprintf(ofp,"Could not allocate space for matrix_mult temp!?!\n"); exit(1);}


/* do all but first and last rows and first last elmt in each row */
for (i=1;i<(n-1);i++)
  for (j=1;j<(n-1);j++)
    {
       *(temp+i*n+j) = *(A+i*n*5+(j-n-1)*5+4) * *(x+i*n+j-n-1);	
       *(temp+i*n+j) += *(A+i*n*5+(j-n)*5+3) * *(x+i*n+j-n);
       *(temp+i*n+j) += *(A+i*n*5+(j-n+1)*5+2) * *(x+i*n+j-n+1);
       *(temp+i*n+j) += *(A+i*n*5+(j-1)*5+1) * *(x+i*n+j-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5) * *(x+i*n+j);
       *(temp+i*n+j) += *(A+i*n*5+j*5+1) * *(x+i*n+j+1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+2) * *(x+i*n+j+n-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+3) * *(x+i*n+j+n);
       *(temp+i*n+j) += *(A+i*n*5+j*5+4) * *(x+i*n+j+n+1);
     }


/* first elm in all but first and last two rows */
j=0;
for (i=2;i<(n-2);i++)
    {
       *(temp+i*n+j) = *(A+i*n*5+(j-n-1)*5+4) * *(x+i*n+j-n-1);	
       *(temp+i*n+j) += *(A+i*n*5+(j-n)*5+3) * *(x+i*n+j-n);
       *(temp+i*n+j) += *(A+i*n*5+(j-n+1)*5+2) * *(x+i*n+j-n+1);
       *(temp+i*n+j) += *(A+i*n*5+(j-1)*5+1) * *(x+i*n+j-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5) * *(x+i*n+j);
       *(temp+i*n+j) += *(A+i*n*5+j*5+1) * *(x+i*n+j+1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+2) * *(x+i*n+j+n-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+3) * *(x+i*n+j+n);
       *(temp+i*n+j) += *(A+i*n*5+j*5+4) * *(x+i*n+j+n+1);
     }

/* last elm in all but first and last two rows */
j=n-1;
for (i=2;i<(n-2);i++)
    {
       *(temp+i*n+j) = *(A+i*n*5+(j-n-1)*5+4) * *(x+i*n+j-n-1);	
       *(temp+i*n+j) += *(A+i*n*5+(j-n)*5+3) * *(x+i*n+j-n);
       *(temp+i*n+j) += *(A+i*n*5+(j-n+1)*5+2) * *(x+i*n+j-n+1);
       *(temp+i*n+j) += *(A+i*n*5+(j-1)*5+1) * *(x+i*n+j-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5) * *(x+i*n+j);
       *(temp+i*n+j) += *(A+i*n*5+j*5+1) * *(x+i*n+j+1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+2) * *(x+i*n+j+n-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+3) * *(x+i*n+j+n);
       *(temp+i*n+j) += *(A+i*n*5+j*5+4) * *(x+i*n+j+n+1);
     }

/* first elm in next to last row */
i=n-2;
j=0;
    {
       *(temp+i*n+j) = *(A+i*n*5+(j-n-1)*5+4) * *(x+i*n+j-n-1);	
       *(temp+i*n+j) += *(A+i*n*5+(j-n)*5+3) * *(x+i*n+j-n);
       *(temp+i*n+j) += *(A+i*n*5+(j-n+1)*5+2) * *(x+i*n+j-n+1);
       *(temp+i*n+j) += *(A+i*n*5+(j-1)*5+1) * *(x+i*n+j-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5) * *(x+i*n+j);
       *(temp+i*n+j) += *(A+i*n*5+j*5+1) * *(x+i*n+j+1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+2) * *(x+i*n+j+n-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+3) * *(x+i*n+j+n);
       *(temp+i*n+j) += *(A+i*n*5+j*5+4) * *(x+i*n+j+n+1);
     }

/* last elm in second row */
i=1;
j=n-1;
    {
       *(temp+i*n+j) = *(A+i*n*5+(j-n-1)*5+4) * *(x+i*n+j-n-1);	
       *(temp+i*n+j) += *(A+i*n*5+(j-n)*5+3) * *(x+i*n+j-n);
       *(temp+i*n+j) += *(A+i*n*5+(j-n+1)*5+2) * *(x+i*n+j-n+1);
       *(temp+i*n+j) += *(A+i*n*5+(j-1)*5+1) * *(x+i*n+j-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5) * *(x+i*n+j);
       *(temp+i*n+j) += *(A+i*n*5+j*5+1) * *(x+i*n+j+1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+2) * *(x+i*n+j+n-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+3) * *(x+i*n+j+n);
       *(temp+i*n+j) += *(A+i*n*5+j*5+4) * *(x+i*n+j+n+1);
     }

/* first row plus first elm in second */
i=0;
for (j=0;j<(n+1);j++)
    {
       *(temp+i*n+j) = *(A+i*n*5+j*5) * *(x+i*n+j);
       *(temp+i*n+j) += *(A+i*n*5+j*5+1) * *(x+i*n+j+1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+2) * *(x+i*n+j+n-1);
       *(temp+i*n+j) += *(A+i*n*5+j*5+3) * *(x+i*n+j+n);
       *(temp+i*n+j) += *(A+i*n*5+j*5+4) * *(x+i*n+j+n+1);
       if ((j-n-1) >= 0)
          {*(temp+i*n+j) += *(A+i*n*5+(j-n-1)*5+4) * *(x+i*n+j-n-1);}
       if ((j-n) >= 0)	
          {*(temp+i*n+j) += *(A+i*n*5+(j-n)*5+3) * *(x+i*n+j-n);}
       if ((j-n+1) >= 0)
          {*(temp+i*n+j) += *(A+i*n*5+(j-n+1)*5+2) * *(x+i*n+j-n+1);}
       if ((j-1) >= 0)
          {*(temp+i*n+j) += *(A+i*n*5+(j-1)*5+1) * *(x+i*n+j-1);}
     }

/* last row plus last elm in prev row */
i=n-1;
for (j=-1;j<n;j++)
    {
       *(temp+i*n+j) = *(A+i*n*5+j*5) * *(x+i*n+j);
       *(temp+i*n+j) += *(A+i*n*5+(j-n-1)*5+4) * *(x+i*n+j-n-1);
       *(temp+i*n+j) += *(A+i*n*5+(j-n)*5+3) * *(x+i*n+j-n);
       *(temp+i*n+j) += *(A+i*n*5+(j-n+1)*5+2) * *(x+i*n+j-n+1);
       *(temp+i*n+j) += *(A+i*n*5+(j-1)*5+1) * *(x+i*n+j-1);
       if ((j+1) < n)
          {*(temp+i*n+j) += *(A+i*n*5+j*5+1) * *(x+i*n+j+1);}
       if ((j+n-1) >= 0)	
          {*(temp+i*n+j) += *(A+i*n*5+j*5+2) * *(x+i*n+j+n-1);}
       if ((j+n) >= 0)
          {*(temp+i*n+j) += *(A+i*n*5+j*5+3) * *(x+i*n+j+n);}
       if ((j+n+1) >= 0)
          {*(temp+i*n+j) += *(A+i*n*5+j*5+4) * *(x+i*n+j+n+1);}

     }

return(temp);             
}





