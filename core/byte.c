#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef CWRITE
#include <fcntl.h>
#endif


#define READ     1
#define WRITE    2
#define MAX_NAME 80


#define SWAP(a,b)       temp=(a); (a)=(b); (b)=temp;


#ifdef CWRITE
static int  fp=0;
#else
static FILE *fp=NULL;
#endif
static int  flag=0;
static char name[MAX_NAME+1];

int bytesw_write=0;
int bytesw_read=0;

/*************************************byte.c***********************************/

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


#ifdef SAFE
  if (*nn<0)
    {printf("byte_reverse() :: n must be positive\n"); abort();}
#endif
  
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
    {printf("byte_open() :: file name has negative length!\n"); abort();}

  if (len>MAX_NAME)
    {printf("byte_open() :: file name too long!\n"); abort();}

  strcpy(name,n);
  strcpy(dirname,n);

  for(i=1;dirname[i]!='\0';i++)
  {
     if(i>0 && dirname[i]=='/')
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
  if (!fp)
    {return;}



#ifdef CWRITE  
  if (close(fp))
    {printf("byte_close() :: couldn't close file!\n"); abort();}

  fp=0;
#else
  if (fclose(fp))
    {printf("byte_close() :: couldn't fclose file!\n"); abort();}

  fp=NULL;
#endif
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
  if (!fp)
    {return;}

#ifdef CWRITE
  lseek(fp, 0, SEEK_SET);
#else
  rewind(fp);
#endif
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

#ifdef SAFE
  if (*n<0)
    {printf("byte_write() :: n must be positive\n"); abort();}
#endif

#ifdef CWRITE
  if (!fp)
    {
      flags = O_WRONLY|O_CREAT;
      mode = S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH;
      if ((fp=open(name, flags, mode))==-1)
	{printf("%s\n",name);
         error_msg_fatal("byte_write() :: open failure!\n");}
      flag=WRITE;
    }

  if (flag==WRITE)
    {write(fp,(char *)buf,sizeof(float)* *n);}
  else
    {printf("byte_write() :: can't write after reading!\n"); abort();}
#else
  if (!fp)
    {
      if (!(fp=fopen(name,"wb")))
        {printf("byte_write() :: fopen failure!\n"); abort();}
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
    {printf("byte_write() :: can't fwrite after freading!\n"); abort();}
#endif
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

#ifdef SAFE
  if (*n<0)
    {printf("byte_reverse() :: n must be positive\n"); abort();}
#endif

#ifdef CWRITE
  if (!fp)
    {
      flags = O_RDONLY;
      mode = S_IRUSR|S_IRGRP|S_IROTH;
      if ((fp=open(name, flags, mode))==-1)
	{printf("%s\n",name);
	 error_msg_fatal("byte_read() :: open failure1!\n");}
      flag=READ;
    }

  if (flag==READ)
    {read(fp,(char *) buf,sizeof(float)* *n);}
  else
    {printf("byte_read() :: can't read after writing!\n"); abort();}
#else
  if (!fp)
    {
      if (!(fp=fopen(name,"rb")))
	{printf("%s\n",name);
         printf("byte_read() :: fopen failure2!\n"); abort();}
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
     }
  else
    {printf("byte_read() :: can't fread after fwriting!\n"); abort();}
#endif
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


#ifdef TEST_DRIVER

main(int argc, char **argv)
{
  int i,j,k;
  int inplace;
  int input_ch;
  char ch, *name;
  FILE *ifp,*ofp;
  float rf[]={1.0,2.0,3.0,4.0};

  if (argc==2)
    {
      inplace=TRUE;
      printf("warning: original file will be destroyed!!!\n"); abort();
      printf("continue? (y/n) : ");
      input_ch = getchar();
      ch=toupper(input_ch);
      if (ch!='Y')
	{exit(0);}
    }
  else if (argc==3)
    {
      inplace=FALSE;
      printf("warning: new file will be created!!!\n");
      printf("continue? (y/n) : ");
      input_ch = getchar();
      ch=toupper(input_ch);
      if (ch!='Y')
	{exit(0);}
    }
  else
    {
      printf("usage: rbyte i/o_file\n");
      printf("or\n");
      printf("usage: rbyte i_file o_file\n");
      exit(0);
    }
  printf("continuing ...\n");

  if (inplace)
    {
      if ((name=tmpnam(NULL))==NULL)
	{printf("tmpname() failure!\n"); exit(0);}
    }
  else
    {
      name = *(argv+2);
    }

  printf("%s,%s\n",*(argv+1),*(argv+2));

  ifp = fopen(*(argv+1),"rb");
  ofp = fopen(name,"wb");

  if (!ifp||!ofp)
    {printf("fopen() failure!\n"); exit(0);}

  fwrite(&rf,sizeof(float),4,ofp);
  fclose(ifp);
  fclose(ofp);

  if (inplace)
    {
      if (remove(*(argv+1))||rename(name,*(argv+1)))
	{printf("remove() or rename() failure!\n"); exit(0);}
    }

  exit(0);
}
#endif
