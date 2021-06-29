#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mpi.h>
#include "gslib.h"
#include "name.h"

#define fNEK_File_open    FORTRAN_UNPREFIXED(nek_file_open, NEK_FILE_OPEN )
#define fNEK_File_read    FORTRAN_UNPREFIXED(nek_file_read, NEK_FILE_READ )
#define fNEK_File_write   FORTRAN_UNPREFIXED(nek_file_write,NEK_FILE_WRITE)
#define fNEK_File_close   FORTRAN_UNPREFIXED(nek_file_close,NEK_FILE_CLOSE)
#define fNEK_File_EOF     FORTRAN_UNPREFIXED(nek_file_eof  ,NEK_FILE_EOF  )

#define READ           0
#define READWRITE      1
#define WRITE          2
#define MAX_NAME       132
#define CB_BUFFER_SIZE 16777216

// Error code
#define NEKIO_RDERR_COUNT  1
#define NEKIO_RDERR_WRONL  2
#define NEKIO_RDERR_NONFL  3
#define NEKIO_RDERR_MPISV  4
#define NEKIO_RDERR_MPIRA  5
#define NEKIO_RDERR_ACCMD  6
#define NEKIO_RDERR_OVRLP  7
#define NEKIO_RDERR_CREAD  8
#define NEKIO_RDERR_GTPOS  8
#define NEKIO_RDERR_FSEEK  9

#define NEKIO_WRERR_COUNT  1
#define NEKIO_WRERR_RDONL  2
#define NEKIO_WRERR_NONFL  3
#define NEKIO_WRERR_MPISV  4
#define NEKIO_WRERR_MPIWA  5
#define NEKIO_WRERR_NSUPP  6

#define SWAP(a,b)       temp=(a); (a)=(b); (b)=temp;

typedef struct NEK_File_handle {
    // General
    char           name[MAX_NAME+1];  // name of the file
    int            mpiio;             // 1 if use mpiio
    int            cbnodes;           // Number of aggregators
    int            num_node;          // Number of compute nodes
    struct comm    comm;              // MPI comm
    struct crystal cr;                // crystal router
    // MPIIO specific
    MPI_File *mpifh;             // mpi file pointer
    int       mpimode;           // MPI_MODE_RDONLY, MPI_MODE_RDWR, or MPI_MODE_WRONLY
    // Byte file specific
    FILE        *file;              // byte file pointer
    int          bmode;             // READ, READWRITE, or WRITE
    struct array tarr;              // Array holding tuples (nektp)
    struct array io2parr;           // Array holding file blocks (nekfb)
    struct array tmp_buf;           // Array holding temporary buffer in ioranks
} nekfh;

// Struct for sarray_transfer
typedef struct NEK_File_block {
    uint proc;
    long long int global_offset;   
    long long int bytes;
    char buf[CB_BUFFER_SIZE];   // data start at the position    
} nekfb;

typedef struct NEK_File_tuple {
    uint iorank;              // the iorank in charge of this process
    long long int start;      // start of the file
    long long int end;        // end of the file
} nektp;

static int handle_n        = 0;
static int handle_max      = 0;
static nekfh **fhandle_arr = 0;

#ifdef UNDERSCORE
  void exitt_();
#else
  void exitt();
#endif

int NEK_File_open(const MPI_Comm fcomm, void *handle, char *filename, int *amode, int *ifmpiio, int *cb_nodes, int nlen)
{
    int i,istat,ierr,ferr;
    int rank,shmrank;
    char dirname[MAX_NAME+1];
    MPI_Comm shmcomm, nodecomm;
    MPI_File *mpi_fh;
    nekfh *nek_fh; 
    struct crystal crs;
    struct comm c;
    struct array tarr    = null_array;
    struct array io2parr = null_array;
    struct array tmp_buf = null_array;
    char cbnodes_str[12];
    char bytemode[4];       // "rb", "rwb", "wb" correspond to bmode = 0,1,2 
    int num_node;
    
    nek_fh = (nekfh*) handle;

    strncpy(nek_fh->name,filename,MAX_NAME);
    for (i=nlen-1; i>0; i--) if ((nek_fh->name)[i] != ' ') break;
    (nek_fh->name)[i+1] = '\0';

    comm_init(&c, fcomm);
    nek_fh->comm = c;
    crystal_init(&crs, &nek_fh->comm);
    nek_fh->cr = crs;
    rank = c.id;

    MPI_Comm_split_type(c.c, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &shmrank);
    MPI_Comm_split(c.c, shmrank, 0, &nodecomm);
    MPI_Comm_size(nodecomm, &num_node);
    MPI_Comm_free(&shmcomm);
    MPI_Comm_free(&nodecomm); 
    nek_fh->num_node = num_node;
    
    nek_fh->cbnodes = *cb_nodes;
    sprintf(cbnodes_str, "%d", ((*cb_nodes == 0) ? num_node : *cb_nodes));

    nek_fh->tarr    = tarr;
    nek_fh->io2parr = io2parr;
    nek_fh->tmp_buf = tmp_buf;

    if (*ifmpiio) {
        // Use MPIIO
        mpi_fh = malloc(sizeof(MPI_File));
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

        ferr = MPI_File_open(c.c,nek_fh->name,nek_fh->mpimode,MPI_INFO_NULL,mpi_fh);
        if (ferr != MPI_SUCCESS) {
            // collective call 
            if (rank == 0) {
                printf("Nek_File_open() :: MPI_File_open failure!\n");
            }
            return 1;
        }
        nek_fh->mpifh = mpi_fh;
    } else {
        // Use byte.c 
        nek_fh->mpiio = 0;
        nek_fh->bmode = *amode;
        switch (*amode) {
            case READ:
                strncpy(bytemode,"rb",4);
                break;
            case READWRITE:
                strncpy(bytemode,"rwb",4);
                break;
            case WRITE:
                strncpy(bytemode,"wb",4);
                break;
        }
        if (!((nek_fh->file)=fopen(nek_fh->name,bytemode))) {
            // collective call 
            if (rank == 0) {
                printf("Nek_File_open() :: fopen failure! Filename: %s\n",nek_fh->name);
            }
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

int is_iorank(int rank, int iorank_interval, int num_iorank) {
    int rank_id = rank/iorank_interval;
    // For every iorank_interval, the first rank will be iorank. (Uniform sampling)
    return (rank % iorank_interval == 0 && rank_id < num_iorank); 
}

long long int get_nbyte(long long int nbyte_g, int num_iorank, int rank, int iorank_interval) {
    if (!is_iorank(rank,iorank_interval,num_iorank)) { return 0; }
    long long int nbyte = nbyte_g/num_iorank;
    int mid = num_iorank-(nbyte_g%num_iorank); // rank/iorank_interval < mid then has nbyte elem, otherwise nbyte+1 elems
    return (rank/iorank_interval < mid) ? nbyte : nbyte+1;
}

long long int get_start_io(long long int start_g, long long int nbyte_g, int num_iorank, int rank, int iorank_interval) {
    if (!is_iorank(rank,iorank_interval,num_iorank)) { return 0; }
    int mid = num_iorank-(nbyte_g%num_iorank);
    int rank_id = rank/iorank_interval;
    return (rank_id < mid) ? start_g+rank_id*(nbyte_g/num_iorank) : start_g+rank_id*(nbyte_g/num_iorank)+(rank_id-mid);
}

/*
   offset: in number of bytes
   count: number of bytes to read
*/
int NEK_File_read(void *handle, void *buf, long long int *count, long long int *offset)
{
    nekfh *nek_fh = (nekfh*) handle;
    int rank;
    // MPI rank on the compute node, compute node id, number of compute nodes, number of compute nodes doing io
    int num_node, num_iorank, iorank_interval;
    int nproc;
    int ferr,ierr_p,ierr_g;
    struct comm c = nek_fh->comm;
    MPI_Comm comm = c.c;
    
    rank = c.id;
    MPI_Comm_size(comm,&nproc);
    num_node   = nek_fh->num_node;
    num_iorank = ((nek_fh->cbnodes >= num_node) || (nek_fh->cbnodes <= 0)) ? num_node : nek_fh->cbnodes;
    iorank_interval = nproc/num_iorank;

    ierr_p = 0;
    ierr_g = 0;

    ierr_p = (*count < 0);
    MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
    if (ierr_g) { return NEKIO_RDERR_COUNT; }

    if (nek_fh->mpiio) {
        // MPIIO read
        ierr_p = (nek_fh->mpimode == MPI_MODE_WRONLY);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_WRONL; }

        ierr_p = (*(nek_fh->mpifh) == NULL);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_NONFL; }

        ferr = MPI_File_set_view(*(nek_fh->mpifh),*offset,MPI_BYTE,MPI_BYTE,"native",MPI_INFO_NULL);
        ierr_p = (ferr != MPI_SUCCESS);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_MPISV; }

        ferr = MPI_File_read_all(*(nek_fh->mpifh),buf,*count,MPI_BYTE,MPI_STATUS_IGNORE);
        ierr_p = (ferr != MPI_SUCCESS);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_MPIRA; }
        
    } else {
        // byte read
        // starting byte and ending byte of global, each iorank, each process, and each batch
        // total number of bytes to read in global, each iorank, and each batch
        long long int start_g, end_g, start_io, end_io, start_p, end_p, start_b, end_b, sid, eid, nbyte_g, nbyte_t, nbyte_b, nbyte;
        int num_pass;
        long long int nr, idx;
        char val;
        char  *tmp_buffer;
        nektp *p, *e;
        nekfb *d, *s; 
        start_g = LLONG_MAX;
        end_g   = LLONG_MIN;
        nbyte_g = 0;
        
        // Check access mode
        ierr_p = (!((nek_fh->bmode)==READ || (nek_fh->bmode)==READWRITE));
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_ACCMD; }
        
        // Determine total number of bytes to read and global start/end index 
        start_p = *offset;
        end_p   = *offset+*count;
        MPI_Allreduce(&start_p,&start_g,1,MPI_LONG_LONG,MPI_MIN,comm);
        MPI_Allreduce(&end_p  ,&end_g  ,1,MPI_LONG_LONG,MPI_MAX,comm);
        MPI_Allreduce(count   ,&nbyte_t,1,MPI_LONG_LONG,MPI_SUM,comm);

        // Determine byte to read for each iorank
        nbyte_g  = end_g-start_g;
        nbyte    = get_nbyte(nbyte_g, num_iorank, rank, iorank_interval);
        start_io = get_start_io(start_g, nbyte_g, num_iorank, rank, iorank_interval);
        end_io   = start_io + nbyte;

        // Check overlapping
        ierr_p = (nbyte_t != nbyte_g);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_OVRLP; }
        
        // Tuple list on each process, add the iorank correspond to current process
        p = array_reserve(nektp, &(nek_fh->tarr), num_iorank);

        p  = (nek_fh->tarr).ptr;
        nr = 0;
        for (int io = 0; io < num_iorank; ++io) {
            sid = get_start_io(start_g, nbyte_g, num_iorank, io*iorank_interval, iorank_interval);
            eid = sid+get_nbyte(nbyte_g, num_iorank, io*iorank_interval, iorank_interval);
            // If overlap
            if (sid <= end_p && start_p <= eid) {
                p->start  = start_p;
                p->end    = end_p;
                p->iorank = io*iorank_interval;
                nr++;
                p = p+1;
            }
        }
        (nek_fh->tarr).n = nr;
        
        // Pass tuple list to io ranks
        sarray_transfer(nektp,&(nek_fh->tarr),iorank,1,&(nek_fh->cr));

        // Send read data in CB_BUFFER_SIZE chunks
        num_pass = (nbyte_g % num_iorank == 0) ? (nbyte_g/num_iorank-1)/CB_BUFFER_SIZE+1 : (nbyte_g/num_iorank)/CB_BUFFER_SIZE+1;
        for (int n_pass = 0; n_pass < num_pass; ++n_pass) {
            nbyte_b = (n_pass < num_pass-1) ? CB_BUFFER_SIZE : nbyte-n_pass*CB_BUFFER_SIZE; 
            start_b = start_io+n_pass*CB_BUFFER_SIZE;
            end_b   = start_io+n_pass*CB_BUFFER_SIZE+nbyte_b;

            // If we are in a io node, read blocks with size <= CB_BUFFER_SIZE to tmp_buf
            if (is_iorank(rank,iorank_interval,num_iorank)) {
                tmp_buffer = array_reserve(char,&(nek_fh->tmp_buf),nbyte_b);
                ierr_p = fseek(nek_fh->file,start_io+n_pass*CB_BUFFER_SIZE,SEEK_SET);
                fread(tmp_buffer,1,nbyte_b,nek_fh->file);
            }

            // Check for file error
            ierr_p = (ierr_p || ferror(nek_fh->file) || feof(nek_fh->file));
            MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
            if (ierr_g) { return NEKIO_RDERR_CREAD; }
            
            // In io ranks, traverse through each byte read before, and put into sarray_transform if
            // appears in the tuple list
            d = array_reserve(nekfb, &(nek_fh->io2parr), nproc);
            d = (nek_fh->io2parr).ptr;
            s = (nek_fh->io2parr).ptr;
            p = (nek_fh->tarr).ptr;
            e = (nek_fh->tarr).ptr;
            nr = 0;
            
            // For each child of iorank
            for (int k = 0; k < (nek_fh->tarr).n; ++k) {
                e = p+k;

                // If the current block overlaps with the processor block, send the blocks back to processes
                if (start_b <= e->end && e->start <= end_b) {
                    sid = (e->start < start_b) ? start_b : e->start;
                    eid = (e->end   > end_b  ) ? end_b   : e->end;
                    s->proc          = e->iorank;
                    s->global_offset = sid;
                    s->bytes         = eid-sid;
                    // Fill buffer
                    idx = s->global_offset-start_b;

                    for (int i = 0; i < s->bytes; ++i) {
                        s->buf[i] = tmp_buffer[idx+i];
                    }
                    nr++;
                    s = s+1;
                }

            }
            (nek_fh->io2parr).n = nr;

            // Pass data in iorank to individual processes
            sarray_transfer(nekfb,&(nek_fh->io2parr),proc,1,&(nek_fh->cr));
     
            // Traverse through io_to_proc_array, and put byte into buffer
            d = (nek_fh->io2parr).ptr;
            s = (nek_fh->io2parr).ptr;

            for (int k = 0; k < (nek_fh->io2parr).n; k++) {
                s = d+k;
                memcpy(buf+s->global_offset-start_p,s->buf,s->bytes);
            }
        }
    }
    return 0;
}

int NEK_File_write(void *handle, void *buf, long long int *count, long long int *offset) 
{   
    int ierr_p, ierr_g, ferr;
    nekfh *nek_fh = (nekfh*) handle;
    struct comm c = nek_fh->comm;
    MPI_Comm comm = c.c;
    ierr_p = 0;
    ierr_g = 0;

    ierr_p = (*count < 0);
    MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
    if (ierr_g) { return NEKIO_WRERR_COUNT; }

    if (nek_fh->mpiio) {
        // Use MPIIO
        ierr_p = (nek_fh->mpimode == MPI_MODE_RDONLY);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_WRERR_RDONL; }

        ierr_p = (*(nek_fh->mpifh) == NULL);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_WRERR_NONFL; }

        ferr = MPI_File_set_view(*(nek_fh->mpifh),*offset,MPI_BYTE,MPI_BYTE,"native",MPI_INFO_NULL);
        ierr_p = (ferr != MPI_SUCCESS);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_WRERR_MPISV; }

        ferr = MPI_File_write_all(*(nek_fh->mpifh),buf,*count,MPI_BYTE,MPI_STATUS_IGNORE);
        ierr_p = (ferr != MPI_SUCCESS);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_WRERR_MPIWA; }
        
    } else {
        // byte write
        // TODO: parallel write
        return NEKIO_WRERR_NSUPP;
        /*
        // Serial write
        if ((nek_fh->bmode)==WRITE || (nek_fh->bmode)==READWRITE) {
            fseek(nek_fh->file,*offset,SEEK_CUR);
            fwrite(buf,1,*count,(nek_fh->file));
            if (ferror(nek_fh->file))
            {
                printf("ABORT: Error writing %s\n",nek_fh->name);
                return;
            }
        } else {
            printf("Nek_File_write() :: invalid access mode.\n");
            return;
        }
        */
    }
    return 0;
}

void NEK_File_close(void *handle, int *ierr)
{
    nekfh *nek_fh = (nekfh*) handle;
    int ferr, rank;
    struct comm c = nek_fh->comm;
    MPI_Comm comm = c.c;
    rank = c.id;

    if (nek_fh->mpiio) {
        ferr = MPI_File_close(nek_fh->mpifh);
        free(nek_fh->mpifh);
        if (ferr != MPI_SUCCESS) {
            // collective call 
            if (rank == 0) {
                printf("Nek_File_close() :: MPI_File_close failure!\n");
            }
            *ierr=1;
        }
    } else {
        if (!nek_fh->file) {
            return;
        }

        if (fclose(nek_fh->file)) {
            // collective call
            if (rank == 0) {
                printf("Nek_File_close() :: couldn't fclose file\n");
            }
            *ierr=1;
            return;
        }
    }

    crystal_free(&(nek_fh->cr));
    array_free(&(nek_fh->tarr));
    array_free(&(nek_fh->io2parr));
    array_free(&(nek_fh->tmp_buf));
    comm_free(&(nek_fh->comm));
    *ierr=0;
}

int NEK_File_EOF(void *handle, int* ifeof) {
    nekfh *nek_fh = (nekfh*) handle;
    void *tmp_buf;
    int eof_g, eof_p, ierr_p, ierr_g;
    MPI_Offset curr, end;
    struct comm c = nek_fh->comm;
    MPI_Comm comm = c.c;
    eof_p = 0;
    eof_g = 0;
    ierr_p = 0;
    ierr_g = 0;

    if (nek_fh->mpiio) {
        ierr_p = (*(nek_fh->mpifh) == NULL);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_NONFL; }

        ierr_p = MPI_File_get_position(*(nek_fh->mpifh),&curr);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_GTPOS; }
        
        ierr_p = MPI_File_seek(*(nek_fh->mpifh),0,MPI_SEEK_END);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_FSEEK; }

        ierr_p = MPI_File_get_position(*(nek_fh->mpifh),&end);
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_GTPOS; }
 
        ierr_p = MPI_File_seek(*(nek_fh->mpifh),curr,MPI_SEEK_SET); // Put individual file pointer back
        MPI_Allreduce(&ierr_p,&ierr_g,1,MPI_INT,MPI_LOR,comm);
        if (ierr_g) { return NEKIO_RDERR_FSEEK; }
        
        eof_p = (curr == end);
        MPI_Allreduce(&eof_p,&eof_g,1,MPI_INT,MPI_LOR,comm);
        *ifeof = eof_g;
    } else {
        // Test read to active eof indicator, safe operation since every
        // fread in nekio is preceded by fseek
        fread(tmp_buf,1,1,nek_fh->file); 
        eof_p = feof(nek_fh->file);
        MPI_Allreduce(&eof_p,&eof_g,1,MPI_INT,MPI_LOR,comm);
        *ifeof = eof_g;
    }

    return 0;
}

// Wrapper for Fortran
void fNEK_File_open(const int *fcomm, char *filename, int *amode, int *ifmpiio, int *cb_nodes, int *handle, int *ierr, int nlen) {
    *ierr = 1;
    comm_ext c = MPI_Comm_f2c(*fcomm);
    if (handle_n == handle_max) {
        handle_max += handle_max/2+1;
        fhandle_arr = trealloc(nekfh*,fhandle_arr,handle_max);
    }
    fhandle_arr[handle_n]= (nekfh*) tmalloc(nekfh,1);
    *ierr = NEK_File_open(c,fhandle_arr[handle_n],filename,amode,ifmpiio,cb_nodes,nlen);
    *handle = handle_n++;
}


void fNEK_File_read(int *handle, long long int *count, long long int *offset, void *buf, int *ierr)
{
    int err_code,rank;
    nekfh *fh = fhandle_arr[*handle];
    rank = (fh->comm).id;
    err_code = NEK_File_read(fh,buf,count,offset);

    switch (err_code) {
        case 0:
            *ierr = 0;
            break;
        case NEKIO_RDERR_COUNT:
            if (rank == 0) { printf("Nek_File_read() :: count must be positive\n"); }
            *ierr = 1;
            break;
        case NEKIO_RDERR_WRONL:
            if (rank == 0) { printf("Nek_File_read :: file is write only\n"); }
            *ierr = 1;
            break;
        case NEKIO_RDERR_NONFL:
            if (rank == 0) { printf("Nek_File_read :: no file opened"); }
            *ierr = 1;
            break;
        case NEKIO_RDERR_MPISV:
            if (rank == 0) { printf("Nek_File_read :: MPI_File_set_view failure!\n"); }
            *ierr = 1;
            break;
        case NEKIO_RDERR_MPIRA:
            if (rank == 0) { printf("Nek_File_read :: MPI_File_read_all failure!\n"); }
            *ierr = 1;
            break;
        case NEKIO_RDERR_ACCMD:
            if (rank == 0) { printf("Nek_File_read :: invalid access mode\n"); }
            *ierr = 1;
            break;
        case NEKIO_RDERR_OVRLP:
            if (rank == 0) { printf("ABORT: nekio doesn't support overlapping partitions! \n"); }
            *ierr = 1;
            break;
        case NEKIO_RDERR_CREAD:
            if (rank == 0) { printf("ABORT: Error reading %s\n",fh->name); }
            *ierr = 1;
            break;
        default:
            *ierr = 1;
    }
}

void fNEK_file_write(int *handle, long long int *count, long long int *offset, void *buf, int *ierr) 
{
    int err_code,rank;
    nekfh *fh = fhandle_arr[*handle];
    rank = (fh->comm).id;
    err_code = NEK_File_write(fh,buf,count,offset);

    switch (err_code) {
        case 0:
            *ierr = 0;
            break;
        case NEKIO_WRERR_COUNT:
            if (rank == 0) { printf("Nek_File_write :: count must be positive\n"); }
            *ierr = 1;
            break;
        case NEKIO_WRERR_RDONL:
            if (rank == 0) { printf("Nek_File_write :: file is read only\n"); }
            *ierr = 1;
            break;
        case NEKIO_WRERR_NONFL:
            if (rank == 0) { printf("Nek_File_write :: no file opened"); }
            *ierr = 1;
            break;
        case NEKIO_WRERR_MPISV:
            if (rank == 0) { printf("Nek_File_write :: MPI_File_set_view failure!\n"); }
            *ierr = 1;
            break;
        case NEKIO_WRERR_MPIWA:
            if (rank == 0) { printf("Nek_File_write :: MPI_File_write_all failure!\n"); }
            *ierr = 1;
            break;
        case NEKIO_WRERR_NSUPP:
            if (rank == 0) { printf("Parallel write currently not supported!"); }
            *ierr = 1;
            break;
        default:
            *ierr = 1;
    } 
}

void fNEK_File_close(int *handle, int *ierr)
{
    NEK_File_close(fhandle_arr[*handle],ierr);
    free(fhandle_arr[*handle]);
    fhandle_arr[*handle] = 0;
}

void fNEK_File_EOF(int *handle, int *ifeof, int *ierr)
{
   int err_code,rank;
   nekfh *fh = fhandle_arr[*handle];
   rank = (fh->comm).id;
   err_code = NEK_File_EOF(fh,ifeof);
    
    switch (err_code) {
        case 0:
            *ierr = 0;
            break;
        case NEKIO_RDERR_NONFL:
            if (rank == 0) { printf("Nek_File_EOF :: no file opened"); }
            *ierr = 1;
            break;
        case NEKIO_RDERR_GTPOS:
            if (rank == 0) { printf("Nek_File_EOF :: NEK_File_get_position error"); }
            *ierr = 1;
            break;
        case NEKIO_RDERR_FSEEK:
            if (rank == 0) { printf("Nek_File_EOF :: NEK_File_seek error"); }
            *ierr = 1;
            break;
        default:
            *ierr = 1;
    }
}
