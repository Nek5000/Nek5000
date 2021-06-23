#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mpi.h>
#include "gslib.h"
#include "name.h"

#define fNEK_File_open    FORTRAN_UNPREFIXED(nek_file_open, NEK_FILE_OPEN )
#define fNEK_File_read    FORTRAN_UNPREFIXED(nek_file_read, NEK_FILE_READ )
#define fNEK_File_write   FORTRAN_UNPREFIXED(nek_file_write,NEK_FILE_WRITE)
#define fNEK_File_close   FORTRAN_UNPREFIXED(nek_file_close,NEK_FILE_CLOSE)

#define READ        0
#define READWRITE   1
#define WRITE       2
#define MAX_NAME    132
#define MAX_FHANDLE 100

#define SWAP(a,b)       temp=(a); (a)=(b); (b)=temp;

typedef struct NEK_File_handle {
    // General
    char      name[MAX_NAME+1];  // name of the file
    int       mpiio;             // 1 if use mpiio
    int       cbnodes;           // Number of aggregators
    MPI_Comm  comm;              // MPI comm
    MPI_Info  info;              // MPI info
    MPI_Comm  shmcomm;           // MPI comm within compute nodes
    MPI_Comm  nodecomm;          // MPI comm across compute nodes
    sint     *cr;                // crystal
    // MPIIO specific
    MPI_File *mpifh;             // mpi file pointer
    int       mpimode;           // MPI_MODE_RDONLY, MPI_MODE_RDWR, or MPI_MODE_WRONLY
    // Byte file specific
    FILE     *file;              // byte file pointer
    int       bmode;             // READ, READWRITE, or WRITE
    char      bytemode[4];       // "rb", "rwb", "wb" correspond to bmode = 0,1,2 
} nekfh;

// Struct for sarray_transfer
typedef struct NEK_File_map {
    uint proc;              // processor id
    // TODO: long long int?
    int count;    // number of elements
    int offset;   // offset from file
} nekfp;

typedef struct NEK_File_data {
    uint proc;
    long long int pos;   // position in the buffer
    uint8_t data;        // data at the position    
} nekfd;

typedef struct NEK_File_tuple {
    uint iorank;              // the iorank in charge of this process
    long long int start;      // start of the file
    long long int end;        // end of the file
} nektp;

static int handle_n                    = 0;
static nekfh *fhandle_arr[MAX_FHANDLE] = {NULL};

#ifdef UNDERSCORE
  void exitt_();
#else
  void exitt();
#endif

int NEK_File_open(const MPI_Comm fcomm, sint *cr, void *handle, char *filename, int *amode, int *ifmpiio, int *cb_nodes, int nlen)
{
    int i,istat,ierr;
    char dirname[MAX_NAME+1];
    MPI_Info info;
    MPI_File *mpi_fh = malloc(sizeof(MPI_File));
    nekfh *nek_fh; 
    char cbnodes_str[12];
    sprintf(cbnodes_str, "%d", *cb_nodes);
    
    nek_fh = (nekfh*) handle;

    strncpy(nek_fh->name,filename,MAX_NAME);
    for (i=nlen-1; i>0; i--) if ((nek_fh->name)[i] != ' ') break;
    (nek_fh->name)[i+1] = '\0';
    nek_fh->cbnodes = *cb_nodes;
    nek_fh->comm = fcomm;
    nek_fh->cr = cr;

    MPI_Info_create(&info);
    MPI_Info_set(info, "cb_nodes", cbnodes_str);
    nek_fh->info = info;

    int shmrank;
    MPI_Comm shmcomm, nodecomm; 
    MPI_Comm_split_type(fcomm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &shmrank);
    MPI_Comm_split(fcomm, shmrank, 0, &nodecomm);
    nek_fh->shmcomm  = shmcomm;
    nek_fh->nodecomm = nodecomm;

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

        MPI_File_open(fcomm,nek_fh->name,nek_fh->mpimode,nek_fh->info,mpi_fh);
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

// TODO: num_ionode, rank long long int
long long int get_nbyte(long long int nbyte_g, int num_ionode, int rank) {
    if (rank >= num_ionode) { return 0; }
    long long int nbyte = nbyte_g/num_ionode;
    int mid = num_ionode-(nbyte_g%num_ionode); // rank < mid then has nbyte elem, otherwise nbyte+1 elems
    return (rank < mid) ? nbyte : nbyte+1;
}

long long int get_start_io(long long int start_g, long long int nbyte_g, int num_ionode, int rank) {
    if (rank >= num_ionode) { return 0; }
    int mid = num_ionode-(nbyte_g%num_ionode);
    return (rank < mid) ? start_g+rank*(nbyte_g/num_ionode) : start_g+rank*(nbyte_g/num_ionode)+(rank-mid);
}

/*
   offset: in number of bytes
   count: number of element (real/floats)
*/
void NEK_File_read(void *handle, void *buf, long long int *count, long long int *offset, int *ierr)
{
    nekfh *nek_fh = (nekfh*) handle;
    MPI_Comm fcomm = nek_fh->comm;

    int rank;
    // MPI rank on the compute node, compute node id, number of compute nodes, number of compute nodes doing io
    int shmrank, nodeid, num_node, num_ionode;
    long long int count_io, count_last; // number of elements to read for node 0,...,num_ionode-2. And num_ionode-1
    int size;

    MPI_Comm_rank(fcomm,&rank);
    MPI_Comm_size(fcomm,&size);
    MPI_Comm_rank(nek_fh->shmcomm, &shmrank);
    MPI_Comm_rank(nek_fh->nodecomm, &nodeid);
    MPI_Comm_size(nek_fh->nodecomm, &num_node);
    num_ionode = num_node/(nek_fh->cbnodes)+1;
    count_io   = *count/num_ionode;
    count_last = *count-count_io*(num_ionode-1);
    
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
        MPI_File_read_all(*(nek_fh->mpifh),buf,*count,MPI_REAL,MPI_STATUS_IGNORE);
    
    } else {
        // byte read
        if (!((nek_fh->bmode)==READ || (nek_fh->bmode)==READWRITE)) {
            printf("Nek_File_read :: invalid access mode\n");
            *ierr=1;
            return;
        }

        // Send info to processor 1 to determine total number of bytes to read
        uint8_t *tmp_buf;
        long long int start_g, end_g, nbyte_g, nbyte;
        struct array tuple_array = null_array;
        nektp *p, *e;
        
        p = array_reserve(nektp, &tuple_array, 1), tuple_array.n = 1;
        p->start  = *offset;
        p->end    = *offset+*count*sizeof(float);
        p->iorank = 0;
        printf("getting info from childs... rank %d, start %lld, end %lld\n", rank,p->start, p->end);

        sarray_transfer(nektp,&tuple_array,iorank,1,get_crystal(nek_fh->cr));
//        printf("!!!!!! After transfer, size of tuple_array is: %zu\n", tuple_array.n);
        
        p = tuple_array.ptr;
        e = tuple_array.ptr;
        start_g = LLONG_MAX;
        end_g   = LLONG_MIN;
        nbyte_g = 0;
        if (rank == 0) {
            for (int pid = 0; pid < size; ++pid) {
                e = p+pid;
//                printf("After sarray_transfer, at rank 0..., start %lld, end %lld\n",e->start, e->end);
                start_g = (start_g < e->start) ? start_g : e->start;
                end_g   = (end_g > e->end) ? end_g : e->end;
            }
//            printf("rank0, updated global start, global end, total_byte: %lld, %lld,%lld\n", start_g, end_g, nbyte_g);
        }
        
        // Pass back global start and global end to processes
        if (rank == 0) {
            for (int pid = 0; pid < size; ++pid) {
                e = p+pid;
                e->start = start_g;
                e->end   = end_g;
            }
        }

        sarray_transfer(nektp,&tuple_array,iorank,1,get_crystal(nek_fh->cr));
        start_g = p->start;
        end_g   = p->end;
//        printf("-------   rank %d, start_g is %lld, end_g is %lld\n", rank, start_g, end_g);

        // Determine byte to read for each iorank
        nbyte_g = end_g-start_g;
        nbyte = get_nbyte(nbyte_g, num_ionode, rank);


        // If we are in a io node, read file to tmp_buf
        long long int start_io, end_io;     // On iorank, the starting byte and ending byte
        start_io = 0;
        end_io = 0;
        start_io = get_start_io(start_g, nbyte_g, num_ionode, rank);
        end_io = start_io + nbyte;
//        printf("rank %d, before fread offset is %lld, read till %lld\n", rank, start_io, end_io);
        if (rank < num_ionode) {
            tmp_buf = (uint8_t*) malloc(nbyte);
            fseek(nek_fh->file,start_io,SEEK_SET);
            fread(tmp_buf,1,nbyte,nek_fh->file);
        }
        
        // Tuple list on each process, add the iorank correspond to current process        
        long long int start_p, end_p, sid, eid;    // On each process, the starting byte and ending byte
        start_p = *offset;
        end_p   = *offset+(*count)*sizeof(float);
//        printf("rank %d, before constrcut tuple list, start_p is %lld, end_p is %lld\n", rank, start_p, end_p);
        struct array tuple_array2 = null_array;  // TODO: duplicate variable
        p = array_reserve(nektp, &tuple_array2, num_ionode);

        int nr = 0;
        p = tuple_array2.ptr;
        for (int io = 0; io < num_ionode; ++io) {
            sid = get_start_io(start_g, nbyte_g, num_ionode, io);
            eid = sid+get_nbyte(nbyte_g, num_ionode, io);
            // If overlap
            if (sid <= end_p && start_p <= eid) {
                p->start  = start_p;
                p->end    = end_p;
                p->iorank = io;
                nr++;
                p = p+1;
            }
        }
        tuple_array2.n = nr;
        
        // Pass tuple list to io ranks
        sarray_transfer(nektp,&tuple_array2,iorank,1,get_crystal(nek_fh->cr));
//        printf(" ==== new tuple list at rank %d with length %zu\n",rank, tuple_array2.n);
        
        // In io ranks, traverse through each byte read before, and put into sarray_transform if
        // appears in the tuple list
        nekfd *d, *s;
        struct array io_to_proc_array = null_array;
        d = array_reserve(nekfd, &io_to_proc_array, nbyte*size); // TODO: over allocate?
        d = io_to_proc_array.ptr;
        s = io_to_proc_array.ptr;
        p = tuple_array2.ptr;
        e = tuple_array2.ptr;
        long long int t = 0;
        long long int idx;
        uint8_t val;
        for (long long int i = start_io; i < end_io; i++) {
            idx = i-start_io;
            val = tmp_buf[idx];
            for (int k = 0; k < tuple_array2.n; ++k) {
                e = p+k;
                // If the byte should be sent to process
                if ((i >= e->start) && (i < e->end)) {
                    s->data = val;
                    s->pos  = i;
                    s->proc = e->iorank;
                    t++;
                    s = s+1;
                }
            }
        }
        io_to_proc_array.n = t;

        // Pass data in iorank to individual processes
        sarray_transfer(nekfd,&io_to_proc_array,proc,1,get_crystal(nek_fh->cr));
        printf(" ==== io_to_proc_array at rank %d with length %zu\n",rank, io_to_proc_array.n);
               
        // Traverse through io_to_proc_array, and put byte into buffer
        d = io_to_proc_array.ptr;
        s = io_to_proc_array.ptr;
        for (int i = 0; i < io_to_proc_array.n; i++) {
            s = d+i;
            ((uint8_t*) buf)[(s->pos-start_p)] = s->data;
        }


        // Free memory
        array_free(&tuple_array);
        array_free(&tuple_array2);
        array_free(&io_to_proc_array);
        if (rank < num_ionode) {
            free(tmp_buf);
        }

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
/*
        // Serial read
        if ((nek_fh->bmode)==READ || (nek_fh->bmode)==READWRITE) {
            fseek(nek_fh->file,*offset,SEEK_SET);
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
*/

        /*
        // Here we assume we have 4 processes, rank 0 will read all the things, and distribute to processes 

        // Send offset and count info to processor 1
        struct array A = null_array;
        struct array B = null_array;
        nekfp *p, *e;
        nekfd *b_ptr, *b_ptr2;

        p = array_reserve(nekfp, &A, 1), A.n = 1;
        p->count  = (int) *count;
        p->offset = (int) *offset;
        p->proc   = 0;   
        printf("getting info from childs... rank %d, count %d, offset %d, from proc %d\n", rank,p->count, p->offset, p->proc);

        sarray_transfer(nekfp,&A,proc,1,get_crystal(nek_fh->cr));
        
        // Accumulate total count
        p = A.ptr;
        e = A.ptr;
        int total_count = 0;
        if (rank == 0) {
            for (int pid = 0; pid < 4; ++pid) {
                e = p + pid;
                printf("accumulating total count... rank 0, process %d, count %d, offset %d, from proc %d\n", pid, e->count, e->offset, e->proc);
                total_count += e->count;
            }
            printf("Total count is %d\n", total_count);
        }
        
        array_reserve(nekfd, &B, total_count), B.n = total_count;
        b_ptr = B.ptr;
        
        // io rank read file and put to B
        p = A.ptr;
        e = A.ptr;
        float* tmp_buf;
        if (rank == 0) {
            for (int pid = 0; pid < 4; ++pid) {
                e = p + pid;
                printf("freading... rank 0, process %d, count %d, offset %d, from proc %d\n", pid, e->count, e->offset, e->proc);
                fseek(nek_fh->file,e->offset,SEEK_SET);
                tmp_buf = malloc(e->count*sizeof(float));
                fread(tmp_buf,sizeof(float),e->count,nek_fh->file);
                
                for (int i = 0; i < e->count; ++i) {
                    b_ptr->data = tmp_buf[i];
                    b_ptr->proc = e->proc;
                    b_ptr->pos  = i;
                    //printf(" i %d, tmp_buf[i] %.2f, proc %d\n ", i, tmp_buf[i], e->proc);
                    b_ptr++;
                }
                free(tmp_buf);
            }
        }

       
        sarray_transfer(nekfp,&B,proc,1,get_crystal(nek_fh->cr));
        
        // Each process receive and put into buffer
        b_ptr = B.ptr;
        for (int i = 0; i < *count; i++) {
            b_ptr2 = b_ptr + i; 
            printf("receive ... rank 1, B's index %d, pos %d, data %.2f, from proc %d\n", i, b_ptr2->pos, b_ptr2->data, b_ptr2->proc);
            ((float*) buf)[b_ptr2->pos] = b_ptr2->data;
        } 

        // Free arrays
        array_free(&A);
        array_free(&B);
        */


/*

    // TODO: how to pass in long long int to MPI functions?
    float sbuf[count_io];
    float rbuf[*count];
    int   recvcounts[num_node];
    int   displs[num_node];
    for (int i = 0; i < num_node; i++) {
        if (i < num_ionode-1) {
            recvcounts[i] = count_io;
            displs[i]     = i*count_io;
        } else if (i == num_ionode-1) {
            recvcounts[i] = count_last;
            displs[i]     = i*count_io;
        } else {
            recvcounts[i] = 0;
            displs[i]     = *count;
        }
    }

        // For the first $cbnodes compute nodes, use a MPI rank to do IO
        if ((nek_fh->bmode)==READ || (nek_fh->bmode)==READWRITE) {
            if ((shmrank == 0) && (nodeid < num_node)) {
                fseek(nek_fh->file,*offset+nodeid*count_io*sizeof(float),SEEK_SET);
                if (nodeid < nek_fh->cbnodes-1) {
                    // If not in the last node, read count_io floats
                    fread(sbuf,sizeof(float),count_io,nek_fh->file);
                    //fseek(nek_fh->file,((num_ionode-nodeid-1)*count_io)*sizeof(float),SEEK_SET); // Move all file pointer to finish reading $count bytes
                } else {
                    // If in the last node, read count_last floats
                    fread(sbuf,sizeof(float),count_last,nek_fh->file);
                }
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
        } else {
            printf("Nek_File_read :: invalid access mode\n");
            *ierr=1;
            return;
        }

        MPI_Gatherv(sbuf,count_io,MPI_FLOAT,rbuf,recvcounts,displs,MPI_FLOAT,0,nek_fh->nodecomm);
        MPI_Bcast(rbuf,*count,MPI_FLOAT,0,nek_fh->comm);
        memcpy(buf,rbuf,*count*sizeof(float));

*/
    }
    *ierr = 0;
}

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
        MPI_File_write_all(*(nek_fh->mpifh),buf,*count,MPI_REAL,MPI_STATUS_IGNORE);
        
    } else {
        // byte write
        if ((nek_fh->bmode)==WRITE || (nek_fh->bmode)==READWRITE) {
            fseek(nek_fh->file,*offset,SEEK_CUR);
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

    *ierr=0;
}


// Wrapper for Fortran
// TODO: put output to the end of argument list!
void fNEK_File_open(const int *fcomm, sint *cr, int *handle, char *filename, int *amode, int *ifmpiio, int *cb_nodes, int *ierr, int nlen) {
    *ierr = 1;
    comm_ext c = MPI_Comm_f2c(*fcomm);
    // TODO: case when *handle > MAX_FHANDLE
    fhandle_arr[handle_n] = (nekfh*) malloc(sizeof(nekfh));
    *ierr = NEK_File_open(c,cr,fhandle_arr[handle_n],filename,amode,ifmpiio,cb_nodes,nlen);
    *handle = handle_n;
    handle_n++;
}


void fNEK_File_read(int *handle, void *buf, long long int *count, long long int *offset, int *ierr)
{
    nekfh *fh = fhandle_arr[*handle];
    NEK_File_read(fh,buf,count,offset,ierr);
}

void fNEK_file_write(int *handle, void *buf, long long int *count, long long int *offset, int *ierr) 
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
    if (fhandle_arr[*handle]->mpiio) {
        free(fhandle_arr[*handle]->mpifh);
    }
    free(fhandle_arr[*handle]);
}


