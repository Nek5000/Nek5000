#ifndef FNAME_H
#define FNAME_H

#ifdef UPCASE
#  define FORTRAN_NAME(low,up) up
#else
#ifdef UNDERSCORE
#  define FORTRAN_NAME(low,up) low##_
#else
#  define FORTRAN_NAME(low,up) low
#endif
#endif

#endif

#define k10_mxm FORTRAN_NAME(k10_mxm, K10_MXM)

void tune_mxm888    (double*, double*, double*);
void tune_mxm8864   (double*, double*, double*);
void tune_mxm6488   (double*, double*, double*);
void tune_mxm101010 (double*, double*, double*);
void tune_mxm1010100(double*, double*, double*);
void tune_mxm1001010(double*, double*, double*);


int k10_mxm(double* a, int* sz1, double* b, int* sz2, double* c, int* sz3)
{
    int m, k, n;

    m = *sz1;
    k = *sz2;
    n = *sz3;

    if ((unsigned)a%16 != 0 || (unsigned)b%16 != 0 || (unsigned)c%16 != 0){
       return 1;}

    // 8,8,8  8,8,64  64,8,8
    if (k == 8){
       if (m == 8){
          if     (n ==  8){tune_mxm888 (a,b,c); return 0;}
          else if(n == 64){tune_mxm8864(a,b,c); return 0;}
       }
       else if (m == 64 && n == 8){tune_mxm6488(a,b,c); return 0;}
    }

    // 10,10,10  10,10,100  100,10,10
    if (k == 10){
       if (m == 10){
          if     (n ==  10){tune_mxm101010 (a,b,c); return 0;}
          else if(n == 100){tune_mxm1010100(a,b,c); return 0;}
       }
       else if (m == 100 && n == 10){tune_mxm1001010(a,b,c); return 0;}
    }

    return 1;
}
