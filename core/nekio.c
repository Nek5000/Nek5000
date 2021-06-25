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

#define READ           0
#define READWRITE      1
#define WRITE          2
#define MAX_NAME       132
#define MAX_FHANDLE    100
#define CB_BUFFER_SIZE 16777216     // default 16 MB. TODO: pass in as user param?

#define SWAP(a,b)       temp=(a); (a)=(b); (b)=temp;

typedef struct NEK_File_handle {
    // General
    char           name[MAX_NAME+1];  // name of the file
    int            mpiio;             // 1 if use mpiio
    int            cbnodes;           // Number of aggregators
    MPI_Comm       comm;              // MPI comm
    MPI_Info       info;              // MPI info
    MPI_Comm       shmcomm;           // MPI comm within compute nodes
    MPI_Comm       nodecomm;          // MPI comm across compute nodes
    struct crystal cr;                // crystal router
    // MPIIO specific
    MPI_File *mpifh;             // mpi file pointer
    int       mpimode;           // MPI_MODE_RDONLY, MPI_MODE_RDWR, or MPI_MODE_WRONLY
    // Byte file specific
    FILE        *file;              // byte file pointer
    int          bmode;             // READ, READWRITE, or WRITE
    struct array tarr;              // Array holding tuples (nektp)
    struct array io2parr;           // Array holding file blocks (nekfb)
} nekfh;

// Struct for sarray_transfer
typedef struct NEK_File_block {
    uint proc;
    long long int global_offset;   
    long long int bytes;
    uint8_t buf[CB_BUFFER_SIZE];   // data start at the position    
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
    int i,istat,ierr;
    char dirname[MAX_NAME+1];
    MPI_Info info;
    MPI_Comm new_comm;
    MPI_File *mpi_fh = malloc(sizeof(MPI_File));
    nekfh *nek_fh; 
    struct crystal crs;
    struct comm c;
    char cbnodes_str[12];
    char bytemode[4];       // "rb", "rwb", "wb" correspond to bmode = 0,1,2 
    int num_node;
    
    nek_fh = (nekfh*) handle;

    strncpy(nek_fh->name,filename,MAX_NAME);
    for (i=nlen-1; i>0; i--) if ((nek_fh->name)[i] != ' ') break;
    (nek_fh->name)[i+1] = '\0';
    MPI_Comm_dup(fcomm, &new_comm);
    nek_fh->comm = new_comm;
    comm_init(&c, new_comm);
    crystal_init(&crs, &c);
    nek_fh->cr = crs;

    int shmrank;
    MPI_Comm shmcomm, nodecomm; 
    MPI_Comm_split_type(nek_fh->comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &shmrank);
    MPI_Comm_split(nek_fh->comm, shmrank, 0, &nodecomm);
    nek_fh->shmcomm  = shmcomm;
    nek_fh->nodecomm = nodecomm;
    MPI_Comm_size(nek_fh->nodecomm, &num_node);
    
    nek_fh->cbnodes = *cb_nodes;
    sprintf(cbnodes_str, "%d", ((*cb_nodes == 0) ? num_node : *cb_nodes));
    MPI_Info_create(&info);
    MPI_Info_set(info, "cb_nodes", cbnodes_str);
    nek_fh->info = info;

    struct array tarr    = null_array;
    struct array io2parr = null_array;
    nek_fh->tarr    = tarr;
    nek_fh->io2parr = io2parr;

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

        MPI_File_open(nek_fh->comm,nek_fh->name,nek_fh->mpimode,nek_fh->info,mpi_fh);
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
            printf("%s\n",nek_fh->name);
            printf("Nek_File_open() :: fopen failure2!\n");
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
void NEK_File_read(void *handle, void *buf, long long int *count, long long int *offset, int *ierr)
{
    nekfh *nek_fh = (nekfh*) handle;

    int rank;
    // MPI rank on the compute node, compute node id, number of compute nodes, number of compute nodes doing io
    int shmrank, nodeid, num_node, num_iorank, iorank_interval;
    int nproc;

    MPI_Comm_rank(nek_fh->comm,&rank);
    MPI_Comm_size(nek_fh->comm,&nproc);
    MPI_Comm_rank(nek_fh->shmcomm, &shmrank);
    MPI_Comm_rank(nek_fh->nodecomm, &nodeid);
    MPI_Comm_size(nek_fh->nodecomm, &num_node);
    num_iorank = ((nek_fh->cbnodes >= num_node) || (nek_fh->cbnodes <= 0)) ? num_node : nek_fh->cbnodes;
    iorank_interval = nproc/num_iorank;

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
            *ierr = 1;
            return;
        }
        MPI_File_set_view(*(nek_fh->mpifh),*offset,MPI_BYTE,MPI_BYTE,"native",nek_fh->info);
        MPI_File_read_all(*(nek_fh->mpifh),buf,*count,MPI_BYTE,MPI_STATUS_IGNORE);
    
    } else {
        // byte read
        uint8_t *tmp_buf;
        // starting byte and ending byte of global, each iorank and each process
        // total number of bytes to read in global, and each iorank
        long long int start_g, end_g, start_io, end_io, start_p, end_p, sid, eid, nbyte_g, nbyte_t, nbyte;
        int num_buffer;
        long long int nr, idx;
        uint8_t val;
        nektp *p, *e;
        nekfb *d, *s; 
        start_g = LLONG_MAX;
        end_g   = LLONG_MIN;
        nbyte_g = 0;
        
        // Check access mode
        if (!((nek_fh->bmode)==READ || (nek_fh->bmode)==READWRITE)) {
            printf("Nek_File_read :: invalid access mode\n");
            *ierr=1;
            return;
        }
        
        // Determine total number of bytes to read and global start/end index 
        start_p = *offset;
        end_p   = *offset+*count;
        MPI_Allreduce(&start_p,&start_g,1,MPI_LONG_LONG,MPI_MIN,nek_fh->comm);
        MPI_Allreduce(&end_p  ,&end_g  ,1,MPI_LONG_LONG,MPI_MAX,nek_fh->comm);
        MPI_Allreduce(count   ,&nbyte_t,1,MPI_LONG_LONG,MPI_SUM,nek_fh->comm);

        // Determine byte to read for each iorank
        nbyte_g = end_g-start_g;
        nbyte   = get_nbyte(nbyte_g, num_iorank, rank, iorank_interval);

        // Check overlapping
        if (nbyte_t == nproc*nbyte_g) {
            num_iorank      = 1; // If all processes read the same chunk, only one io rank do the read
            iorank_interval = nproc;
            nbyte           = get_nbyte(nbyte_g, num_iorank, rank, iorank_interval);
        } 

        // If we are in a io node, read file to tmp_buf
        start_io = get_start_io(start_g, nbyte_g, num_iorank, rank, iorank_interval);
        end_io   = start_io + nbyte;
        if (is_iorank(rank,iorank_interval,num_iorank)) {
            tmp_buf = (uint8_t*) malloc(nbyte);
            fseek(nek_fh->file,start_io,SEEK_SET);
            fread(tmp_buf,1,nbyte,nek_fh->file);
            fseek(nek_fh->file,0,SEEK_SET);        // Move file pointer back
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
        }
        
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
        
        // In io ranks, traverse through each byte read before, and put into sarray_transform if
        // appears in the tuple list
        d = array_reserve(nekfb, &(nek_fh->io2parr), ((end_io-start_io-1)/CB_BUFFER_SIZE+1)*nproc);
        d = (nek_fh->io2parr).ptr;
        s = (nek_fh->io2parr).ptr;
        p = (nek_fh->tarr).ptr;
        e = (nek_fh->tarr).ptr;
        nr = 0;

        // For each child of iorank
        for (int k = 0; k < (nek_fh->tarr).n; ++k) {
            e = p+k;
            sid = (e->start < start_io) ? start_io : e->start;
            eid = (e->end   > end_io  ) ? end_io   : e->end;
            num_buffer = (eid-sid-1)/CB_BUFFER_SIZE+1;
            // Send num_buffer blocks back to processes
            for (int bi = 0; bi < num_buffer; ++bi) {
                s->proc          = e->iorank;
                s->global_offset = sid+bi*CB_BUFFER_SIZE;
                s->bytes         = (bi < num_buffer-1) ? CB_BUFFER_SIZE : (eid-sid)-(num_buffer-1)*CB_BUFFER_SIZE;
                // Fill buffer
                idx = s->global_offset-start_io;   // start idx in tmp_buf
                for (int i = 0; i < s->bytes; ++i) {
                    s->buf[i] = tmp_buf[idx+i];
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
        

        // Free memory
        if (is_iorank(rank,iorank_interval,num_iorank)) {
            free(tmp_buf);
        }
    }
    *ierr = 0;
}

// TODO: parallel write
void NEK_File_write(void *handle, void *buf, long long int *count, long long int *offset, int *ierr) 
{
    nekfh *nek_fh = (nekfh*) handle;

    if (*count < 0) {
        printf("Nek_File_write() :: count must be positive\n");
        *ierr = 1;
        return;
    }

    // TODO: didn't test
    if (nek_fh->mpiio) {
        // Use MPIIO
         if (nek_fh->mpimode == MPI_MODE_RDONLY) {
            printf("Nek_File_write :: file is read only\n");
            *ierr = 1;
            return;
        }
        if (*(nek_fh->mpifh) == NULL) {
            printf("Nek_File_write :: no file opened\n");
            *ierr = 1;
            return;
        }
        MPI_File_set_view(*(nek_fh->mpifh),*offset,MPI_BYTE,MPI_BYTE,"native",nek_fh->info);
        MPI_File_write_all(*(nek_fh->mpifh),buf,*count,MPI_BYTE,MPI_STATUS_IGNORE);
        
    } else {
        // byte write
        if ((nek_fh->bmode)==WRITE || (nek_fh->bmode)==READWRITE) {
            fseek(nek_fh->file,*offset,SEEK_CUR);
            fwrite(buf,1,*count,(nek_fh->file));
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

    *ierr=0;
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
    nekfh *fh = fhandle_arr[*handle];
    NEK_File_read(fh,buf,count,offset,ierr);
}

void fNEK_file_write(int *handle, long long int *count, long long int *offset, void *buf, int *ierr) 
{
    nekfh *fh = fhandle_arr[*handle];
    NEK_File_write(fh,buf,count,offset,ierr);
}

void fNEK_File_close(int *handle, int *ierr)
{
    nekfh *fh = fhandle_arr[*handle];
    NEK_File_close(fh,ierr);
    MPI_Comm_free(&(fh->shmcomm));
    MPI_Comm_free(&(fh->nodecomm)); 
    crystal_free(&(fh->cr));
    array_free(&(fh->tarr));
    array_free(&(fh->io2parr));
    if (fhandle_arr[*handle]->mpiio) {
        free(fhandle_arr[*handle]->mpifh);
    }
    free(fhandle_arr[*handle]);
    fhandle_arr[*handle] = 0;
}

