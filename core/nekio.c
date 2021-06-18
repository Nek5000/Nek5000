#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mpi.h>
#include "gslib.h"
#include "name.h"

#define fNEK_File_open    FORTRAN_UNPREFIXED(nek_file_open, NEK_FILE_OPEN )
#define NEK_File_read     FORTRAN_UNPREFIXED(nek_file_read, NEK_FILE_READ )
#define NEK_File_write    FORTRAN_UNPREFIXED(nek_file_write,NEK_FILE_WRITE)
#define NEK_File_close    FORTRAN_UNPREFIXED(nek_file_close,NEK_FILE_CLOSE)

#define READ       0
#define READWRITE  1
#define WRITE      2
#define MAX_NAME 132

#define SWAP(a,b)       temp=(a); (a)=(b); (b)=temp;

typedef struct NEK_File_handle {
    // General
    char      name[MAX_NAME+1];  // name of the file
    int       mpiio;             // 1 if use mpiio
    int       cbnodes;           // TODO: to implement
    MPI_Comm  comm;             // MPI comm
    // MPIIO specific
    MPI_File *mpifh;             // mpi file pointer
    int       mpimode;           // MPI_MODE_RDONLY, MPI_MODE_RDWR, or MPI_MODE_WRONLY
    // Byte file specific
    FILE     *file;              // byte file pointer
    int       bmode;             // READ, READWRITE, or WRITE
    char      bytemode[4];       // "rb", "rwb", "wb" correspond to bmode = 0,1,2 
} nekfh;


#ifdef UNDERSCORE
  void exitt_();
#else
  void exitt();
#endif

int NEK_File_open(const MPI_Comm fcomm, void *handle, char *filename, int *amode, int *ifmpiio, int *cb_nodes, int nlen)
{
    int i,istat,ierr;
    char dirname[MAX_NAME+1];
    MPI_File *mpi_fh = malloc(sizeof(MPI_File));
    nekfh *nek_fh; 

    *((nekfh*)handle) = *((nekfh*) malloc(sizeof(nekfh)));
    nek_fh = (nekfh*) handle;

    strncpy(nek_fh->name,filename,MAX_NAME);
    for (i=nlen-1; i>0; i--) if ((nek_fh->name)[i] != ' ') break;
    (nek_fh->name)[i+1] = '\0';
    nek_fh->cbnodes = *cb_nodes;
    nek_fh->comm = fcomm;
    
    if (*ifmpiio) {
        // Use MPIIO
        nek_fh->mpiio = 1;
        switch (*amode) {
            case READ:
                nek_fh->mpimode = MPI_MODE_RDONLY;
                break;
            case READWRITE:
                nek_fh->mpimode = MPI_MODE_RDWR;
                break;
            case WRITE:
                nek_fh->mpimode = MPI_MODE_WRONLY;
                break;
            default:
                nek_fh->mpimode = MPI_MODE_RDONLY;
        }

        MPI_File_open(fcomm,nek_fh->name,nek_fh->mpimode,MPI_INFO_NULL,mpi_fh);
        nek_fh->mpifh = mpi_fh;
    } else {
        // Use byte.c 
        // TODO: clean up bytemode
        nek_fh->mpiio = 0;
        nek_fh->bmode = *amode;
        switch (*amode) {
            case READ:
                strncpy(nek_fh->bytemode,"rb",4);
                break;
            case READWRITE:
                strncpy(nek_fh->bytemode,"rwb",4);
                break;
            case WRITE:
                strncpy(nek_fh->bytemode,"wb",4);
                break;
        }
        if (!((nek_fh->file)=fopen(nek_fh->name,nek_fh->bytemode))) {
            printf("%s\n",nek_fh->name);
            printf("Nek_File_read() :: fopen failure2!\n");
            return 1;
        }
        for (i=nlen-1; i>0; i--) if ((nek_fh->name)[i] == '/') break;
        if (i > 0) {
            strncpy(dirname,nek_fh->name,i);
            dirname[i] = '\0';
            istat = mkdir(dirname,0755);
        }
    }

    return 0;
}

// Wrapper for Fortran
void fNEK_File_open(const int *fcomm, void *handle, char *filename, int *amode, int *ifmpiio, int *cb_nodes, int *ierr, int nlen) {
    *ierr = 1;
    comm_ext c = MPI_Comm_f2c(*fcomm);
    *ierr = NEK_File_open(c,handle,filename,amode,ifmpiio,cb_nodes,nlen);
}



void NEK_File_read(void *handle, void *buf, long long int *count, long long int *offset, int *ierr)
{
    nekfh *nek_fh = (nekfh*) handle;

    if (*count < 0) {
        printf("Nek_File_read() :: count must be positive\n");
        *ierr = 1;
        return;
    }
    
    if (nek_fh->mpiio) {
        // MPIIO read
        if (nek_fh->mpimode == MPI_MODE_WRONLY) {
            printf("Nek_File_read :: file is write only\n");
            *ierr = 1;
            return;
        }
        if (*(nek_fh->mpifh) == NULL) {
            printf("Nek_File_read :: no file opened");
        }
        MPI_File_set_view(*(nek_fh->mpifh),*offset,MPI_BYTE,MPI_BYTE,"native",MPI_INFO_NULL);
        MPI_File_read_all(*(nek_fh->mpifh),buf,*count,MPI_REAL,MPI_STATUS_IGNORE);
    
    } else {
        // byte read
        if ((nek_fh->bmode)==READ || (nek_fh->bmode)==READWRITE) {
            // TODO: deal with reverse byte
            fseek(nek_fh->file,(*offset)*sizeof(float),SEEK_CUR);
            fread(buf,sizeof(float),*count,nek_fh->file);
            if (ferror(nek_fh->file)) {
                printf("ABORT: Error reading %s\n",nek_fh->name);
                *ierr=1;
                return;
            }
            else if (feof(nek_fh->file)) {
                printf("ABORT: EOF found while reading %s\n",nek_fh->name);
                *ierr=1;
                return;
            }
        } else {
            printf("Nek_File_read :: invalid access mode\n");
            *ierr=1;
            return;
        }

    }
    *ierr = 0;
}

void NEK_file_write(void *handle, void *buf, long long int *count, long long int *offset, int *ierr) 
{
    nekfh *nek_fh = (nekfh*) handle;

    if (*count < 0) {
        printf("Nek_File_write() :: count must be positive\n");
        *ierr = 1;
        return;
    }
    
    if (nek_fh->mpiio) {
        // Use MPIIO
         if (nek_fh->mpimode == MPI_MODE_RDONLY) {
            printf("Nek_File_write :: file is read only\n");
            *ierr = 1;
            return;
        }
        if (*(nek_fh->mpifh) == NULL) {
            printf("Nek_File_write :: no file opened\n");
        }
        MPI_File_write_all(*(nek_fh->mpifh),buf,*count,MPI_REAL,MPI_STATUS_IGNORE);
        
    } else {
        // byte write
        if ((nek_fh->bmode)==WRITE || (nek_fh->bmode)==READWRITE) {
            fseek(nek_fh->file,(*offset)*sizeof(float),SEEK_SET);
            fwrite(buf,sizeof(float),*count,(nek_fh->file));
            if (ferror(nek_fh->file))
            {
                printf("ABORT: Error writing %s\n",nek_fh->name);
                *ierr=1;
                return;
            }
        } else {
            printf("Nek_File_write() :: invalid access mode.\n");
            *ierr=1;
            return;
        }
    }

    *ierr = 0;
}

void NEK_File_close(void *handle, int *ierr)
{
    nekfh *nek_fh = (nekfh*) handle;
    int rank;
    MPI_Comm_rank(nek_fh->comm,&rank);

    if (nek_fh->mpiio) {
        MPI_File_close(nek_fh->mpifh);
    } else {
        if (!nek_fh->file) {
            return;
        }

        if (fclose(nek_fh->file)) {
            printf("Nek_File_close() :: couldn't fclose file\n");
            *ierr=1;
            return;
        }
    }
    // TODO: free handle!
    // free(nek_fh->mpifh)
    // free(nek_fh);
    *ierr=0;
}
