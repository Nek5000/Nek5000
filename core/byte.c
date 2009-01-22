#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>


#define READ     1
#define WRITE    2
#define MAX_NAME 80


#define SWAP(a,b)       temp=(a); (a)=(b); (b)=temp;


static FILE *fp=NULL;
static int  flag=0;
static char name[MAX_NAME+1];

int bytesw_write=0;
int bytesw_read=0;

/*************************************byte.c***********************************/

void cexitt()
{
#ifdef UNDERSCORE
  exitt_();
#else
  exitt();
#endif

}

void
#ifdef UPCASE
BYTE_REVERSE (float *buf, int *nn)
#elif  UNDERSCORE
byte_reverse_(float *buf, int *nn)
#else
byte_reverse(float *buf, int *nn)
#endif
{
  int n;
  char temp, *ptr;

  if (*nn<0)
  {
    printf("byte_reverse() :: n must be positive\n"); 
    cexitt();
  }
  
  for (ptr=(char *)buf,n=*nn; n--; ptr+=4)
  {
     SWAP(ptr[0],ptr[3])
     SWAP(ptr[1],ptr[2])
  }
}


void
#ifdef UPCASE
BYTE_OPEN (char *n)
#elif  UNDERSCORE
byte_open_(char *n)
#else
byte_open(char *n)
#endif
{
  int  i,len,istat;
  char slash;
  char dirname[MAX_NAME+1];

  len = strlen(n);
  
  if (len<0)
  {
    printf("byte_open() :: file name has negative length!\n"); 
    cexitt();
  }

  if (len>MAX_NAME)
  {
    printf("byte_open() :: file name too long!\n"); 
    cexitt();
  }

  strcpy(name,n);
  strcpy(dirname,n);

  for (i=1;dirname[i]!='\0';i++)
  {
     if (i>0 && dirname[i]=='/')
     {
       slash = name[i];
       dirname[i] = '\0';
       istat = mkdir(dirname,0755);
     }
  }
}


void
#ifdef UPCASE
BYTE_CLOSE ()
#elif  UNDERSCORE
byte_close_()
#else
byte_close()
#endif
{
  if (!fp) return;

  if (fclose(fp))
  {
    printf("byte_close() :: couldn't fclose file!\n");
    cexitt();
  }

  fp=NULL;
}


void
#ifdef UPCASE
BYTE_REWIND ()
#elif  UNDERSCORE
byte_rewind_()
#else
byte_rewind()
#endif
{
  if (!fp) return;

  rewind(fp);
}


void
#ifdef UPCASE
BYTE_WRITE (float *buf, int *n)
#elif  UNDERSCORE
byte_write_(float *buf, int *n)
#else
byte_write(float *buf, int *n)
#endif
{
  int flags;
  mode_t mode;

  if (*n<0)
  {
    printf("byte_write() :: n must be positive\n"); 
    cexitt();
  }

  if (!fp)
  {
    if (!(fp=fopen(name,"wb")))
    {
      printf("byte_write() :: fopen failure!\n"); 
      cexitt();
    }
    flag=WRITE;
  }

  if (flag==WRITE)
    {
      if (bytesw_write == 1)
#ifdef UPCASE
        BYTE_REVERSE (buf,n);
#elif  UNDERSCORE
        byte_reverse_(buf,n);
#else
        byte_reverse (buf,n);
#endif
      fwrite(buf,sizeof(float),*n,fp);
    }
  else
  {
      printf("byte_write() :: can't fwrite after freading!\n"); 
      cexitt();
  }
}


void
#ifdef UPCASE
BYTE_READ (float *buf, int *n)
#elif  UNDERSCORE
byte_read_(float *buf, int *n)
#else
byte_read(float *buf, int *n)
#endif
{
  int flags;
  mode_t mode;

  if (*n<0)
    {printf("byte_read() :: n must be positive\n"); cexitt();}

  if (!fp)
  {
     if (!(fp=fopen(name,"rb")))
     {
        printf("%s\n",name);
        printf("byte_read() :: fopen failure2!\n"); cexitt();
     }
     flag=READ;
  }

  if (flag==READ)
  {
     if (bytesw_read == 1)
#ifdef UPCASE
        BYTE_REVERSE (buf,n);
#elif  UNDERSCORE
        byte_reverse_(buf,n);
#else
        byte_reverse (buf,n);
#endif
     fread(buf,sizeof(float),*n,fp);
     if (ferror(fp))
     {
       printf("ABORT: Error reading %s\n",name);
       cexitt();
     }
     else if (feof(fp))
     {
       printf("ABORT: EOF found while reading %s\n",name);
       cexitt();
     }

  }
  else
  {
     printf("byte_read() :: can't fread after fwriting!\n"); 
     cexitt();
  }
}


void
#ifdef UPCASE
SET_BYTESW_WRITE (int *pa)
#elif  UNDERSCORE
set_bytesw_write_(int *pa)
#else
set_bytesw_write (int *pa)
#endif
{
    if (*pa != 0)
       bytesw_write = 1;
    else
       bytesw_write = 0;
}


void
#ifdef UPCASE
SET_BYTESW_READ (int *pa)
#elif  UNDERSCORE
set_bytesw_read_(int *pa)
#else
set_bytesw_read (int *pa)
#endif
{
    if (*pa != 0)
       bytesw_read = 1;
    else
       bytesw_read = 0;
}


void
#ifdef UPCASE
GET_BYTESW_WRITE (int *pa)
#elif  UNDERSCORE
get_bytesw_write_(int *pa)
#else
get_bytesw_write (int *pa)
#endif
{
    *pa = bytesw_write;
}


void
#ifdef UPCASE
GET_BYTESW_READ (int *pa)
#elif  UNDERSCORE
get_bytesw_read_(int *pa)
#else
get_bytesw_read (int *pa)
#endif
{
    *pa = bytesw_read;
}
