#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include "gslib.h"
#include "name.h"

#define fNEK_File_open FORTRAN_UNPREFIXED(nek_file_open, NEK_FILE_OPEN)
#define fNEK_File_read FORTRAN_UNPREFIXED(nek_file_read, NEK_FILE_READ)
#define fNEK_File_write FORTRAN_UNPREFIXED(nek_file_write, NEK_FILE_WRITE)
#define fNEK_File_close FORTRAN_UNPREFIXED(nek_file_close, NEK_FILE_CLOSE)
#define fNEK_File_EOF FORTRAN_UNPREFIXED(nek_file_eof, NEK_FILE_EOF)

#define READ 0
#define READWRITE 1
#define WRITE 2
#define MAX_NAME 132
#define CB_BUFFER_SIZE 16777216

// Error code
#define NEKIO_RDERR_COUNT 1
#define NEKIO_RDERR_WRONL 2
#define NEKIO_RDERR_NONFL 3
#define NEKIO_RDERR_MPISV 4
#define NEKIO_RDERR_MPIRA 5
#define NEKIO_RDERR_ACCMD 6
#define NEKIO_RDERR_OVRLP 7
#define NEKIO_RDERR_CREAD 8
#define NEKIO_RDERR_GTPOS 8
#define NEKIO_RDERR_FSEEK 9
#define NEKIO_RDERR_USRBUF 10

#define NEKIO_WRERR_COUNT 1
#define NEKIO_WRERR_RDONL 2
#define NEKIO_WRERR_NONFL 3
#define NEKIO_WRERR_MPISV 4
#define NEKIO_WRERR_MPIWA 5
#define NEKIO_WRERR_NSUPP 6

#define SWAP(a, b) \
  temp = (a);      \
  (a) = (b);       \
  (b) = temp;

static char tmp_buffer[CB_BUFFER_SIZE];

typedef struct NEK_File_handle {
  // General
  char *name;
  int mpiio;
  int cbnodes;  // Number of aggregators
  int num_node;
  struct comm comm;

  int agg_rank;
  struct comm aggcomm;

  struct crystal cr;

  // MPIIO specific
  MPI_File *mpifh;
  int mode;

  FILE *file;
  struct array tarr;     // Array holding tuples (nektp)
  struct array io2parr;  // Array holding file blocks (nekfb)
} nekfh;

typedef struct NEK_File_block {
  uint proc;
  long long int offset;
  long long int bytes;
  char buf[CB_BUFFER_SIZE];
} nekfb;

typedef struct NEK_File_tuple {
  uint proc;
  long long int start;  // byte start of the file
  long long int end;    // byte end of the file
} nektp;

static int handle_n = 0;
static int handle_max = 0;
static nekfh **fhandle_arr = 0;

int NEK_File_open(const MPI_Comm fcomm, void *handle, char *filename, int amode,
                  int ifmpiio, int cb_nodes) {
  int i, istat;
  int ferr;
  int ierr = 0;
  int rank, shmrank;
  MPI_Comm shmcomm, nodecomm;
  MPI_File *mpi_fh;
  nekfh *nek_fh;
  struct crystal crs;
  struct comm c;
  struct array tarr = null_array;
  struct array io2parr = null_array;
  char cbnodes_str[12];
  int num_node;

  nek_fh = (nekfh *)handle;
  nek_fh->name = filename;

  comm_init(&c, fcomm);
  nek_fh->comm = c;
  rank = c.id;

  if (!ifmpiio) {
    const int ratio = c.np / cb_nodes;
    const int aggid = rank / ratio;

    MPI_Comm aggcomm;
    MPI_Comm_split(c.c, aggid, 0, &aggcomm);
    comm_init(&nek_fh->aggcomm, aggcomm);
    if (nek_fh->aggcomm.id == 0) nek_fh->agg_rank = rank;
    MPI_Bcast(&nek_fh->agg_rank, 1, MPI_INT, 0, nek_fh->aggcomm.c);
  }

  MPI_Comm_split_type(c.c, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
  MPI_Comm_rank(shmcomm, &shmrank);
  MPI_Comm_split(c.c, shmrank, 0, &nodecomm);
  MPI_Comm_size(nodecomm, &num_node);
  MPI_Comm_free(&shmcomm);
  MPI_Comm_free(&nodecomm);
  nek_fh->num_node = num_node;

  crystal_init(&crs, &nek_fh->aggcomm);
  nek_fh->cr = crs;

  nek_fh->cbnodes = cb_nodes;
  sprintf(cbnodes_str, "%d", ((cb_nodes == 0) ? num_node : cb_nodes));

  nek_fh->tarr = tarr;
  nek_fh->io2parr = io2parr;

  if (ifmpiio) {
    // Use MPIIO
    mpi_fh = malloc(sizeof(MPI_File));
    nek_fh->mpiio = 1;
    switch (amode) {
      case READ:
        nek_fh->mode = MPI_MODE_RDONLY;
        break;
      case READWRITE:
        nek_fh->mode = MPI_MODE_RDWR;
        break;
      case WRITE:
        nek_fh->mode = MPI_MODE_WRONLY;
        break;
      default:
        nek_fh->mode = MPI_MODE_RDONLY;
    }

    ferr =
        MPI_File_open(c.c, nek_fh->name, nek_fh->mode, MPI_INFO_NULL, mpi_fh);
    if (ferr != MPI_SUCCESS) {
      // collective call
      if (rank == 0) {
        printf("Nek_File_open() :: MPI_File_open failure!\n");
      }
      return 1;
    }
    nek_fh->mpifh = mpi_fh;
  } else {
    nek_fh->mpiio = 0;
    nek_fh->mode = amode;
    char bytemode[4] = {0};
    switch (amode) {
      case READ:
        strncpy(bytemode, "rb", 2);
        break;
      case READWRITE:
        strncpy(bytemode, "rwb", 3);
        break;
      case WRITE:
        strncpy(bytemode, "wb", 2);
        break;
    }
    if (!((nek_fh->file) = fopen(nek_fh->name, bytemode))) {
      // collective call
      if (rank == 0) {
        printf("Nek_File_open() :: fopen failure! Filename: %s\n",
               nek_fh->name);
      }
      return 1;
    }
    {
      int i;
      const int nlen = strlen(nek_fh->name);
      char *dirname = (char *)calloc(nlen, sizeof(char));

      for (i = nlen - 1; i > 0; i--)
        if ((nek_fh->name)[i] == '/') break;
      if (i > 0) strncpy(dirname, nek_fh->name, i);
      dirname[i] = '\0';
      istat = mkdir(dirname, 0755);

      free(dirname);
    }
  }

  return 0;
}

int NEK_File_read(void *handle, void *buf, long long int count,
                  long long int offset) {
  nekfh *nek_fh = (nekfh *)handle;
  int rank;
  // MPI rank on the compute node, compute node id, number of compute nodes,
  // number of compute nodes doing io
  int num_node, num_iorank;
  int nproc;
  int ferr;
  int ierr = 0;
  struct comm c = nek_fh->comm;
  MPI_Comm comm = c.c;

  rank = c.id;
  MPI_Comm_size(comm, &nproc);
  num_node = nek_fh->num_node;
  num_iorank = nek_fh->cbnodes;

  ierr = (count < 0);
  MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
  if (ierr) {
    return NEKIO_RDERR_COUNT;
  }

  if (nek_fh->mpiio) {
    ierr = (nek_fh->mode == MPI_MODE_WRONLY);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_RDERR_WRONL;
    }

    ierr = (*(nek_fh->mpifh) == NULL);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_RDERR_NONFL;
    }

    ferr = MPI_File_set_view(*(nek_fh->mpifh), offset, MPI_BYTE, MPI_BYTE,
                             "native", MPI_INFO_NULL);
    ierr = (ferr != MPI_SUCCESS);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_RDERR_MPISV;
    }

    ferr = MPI_File_read_all(*(nek_fh->mpifh), buf, count, MPI_BYTE,
                             MPI_STATUS_IGNORE);
    ierr = (ferr != MPI_SUCCESS);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_RDERR_MPIRA;
    }
  } else {
    ierr = (!((nek_fh->mode) == READ || (nek_fh->mode) == READWRITE));
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_RDERR_ACCMD;
    }

    const long long int start_p = offset;
    const long long int end_p = start_p + count;

    // Add aggregator childs
    array_reserve(nektp, &(nek_fh->tarr), 1);
    (nek_fh->tarr).n = 0;
    if(count) {
      nektp *p = (nek_fh->tarr).ptr;
      p[0].start = start_p;
      p[0].end = end_p;
      p[0].proc = 0; // aggregator rank (always root of nek_fh->aggcomm.c)
      (nek_fh->tarr).n = 1;
    }
    sarray_transfer(nektp, &(nek_fh->tarr), proc, 1, &(nek_fh->cr));
    const int nchilds_agg = (nek_fh->tarr).n;

    // Determine byte to read for each aggregator
    long long int nbyte;
    MPI_Allreduce(&count, &nbyte, 1, MPI_LONG_LONG_INT, MPI_SUM,
                  nek_fh->aggcomm.c);

    // read + distribute data in chunks from aggreator to consumer
    int num_pass = (nbyte + CB_BUFFER_SIZE - 1) / CB_BUFFER_SIZE;
    MPI_Bcast(&num_pass, 1, MPI_INT, 0, nek_fh->aggcomm.c);

    for (int n_pass = 0; n_pass < num_pass; ++n_pass) {
      long long int nbyte_b = (n_pass < num_pass - 1)
                                  ? CB_BUFFER_SIZE
                                  : nbyte - n_pass * CB_BUFFER_SIZE;
      MPI_Bcast(&nbyte_b, 1, MPI_LONG_LONG_INT, 0, nek_fh->aggcomm.c);
      long long int start_b = start_p + n_pass * CB_BUFFER_SIZE;
      MPI_Bcast(&start_b, 1, MPI_LONG_LONG_INT, 0, nek_fh->aggcomm.c);
      long long int end_b = start_b + nbyte_b;

      ierr = 0;
      if (rank == nek_fh->agg_rank) {
        ierr = fseek(nek_fh->file, start_b, SEEK_SET);
        fread(tmp_buffer, 1, nbyte_b, nek_fh->file);
      }
      ierr = (ierr || ferror(nek_fh->file) || feof(nek_fh->file));
      MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX,
                    nek_fh->aggcomm.c);
      if (ierr) {
        return NEKIO_RDERR_CREAD;
      }

      // Distribute chunk to matching childs
      array_reserve(nekfb, &(nek_fh->io2parr), nchilds_agg);
      int nr = 0;
      for (int k = 0; k < nchilds_agg; ++k) {
        nektp *child = (nek_fh->tarr).ptr;

        // chunk overlaps (partly or completely) with child
        if (start_b <= child[k].end && child[k].start <= end_b) {
          const long long int sid =
              (child[k].start < start_b) ? start_b : child[k].start;
          const long long int eid =
              (child[k].end > end_b) ? end_b : child[k].end;

          nekfb *s = (nek_fh->io2parr).ptr;
          s[nr].proc = child[k].proc;
          s[nr].offset = sid;
          s[nr].bytes = eid - sid;

          const long long int idx = s[nr].offset - start_b;
          memcpy(s[nr].buf, tmp_buffer + idx, s[nr].bytes);

          nr++;
        }
      }
      (nek_fh->io2parr).n = nr;
      sarray_transfer(nekfb, &(nek_fh->io2parr), proc, 1, &(nek_fh->cr));

      // copy back into user buffer
      {
        long long int sum = 0;
        for (int k = 0; k < (nek_fh->io2parr).n; k++) {
          nekfb *s = (nek_fh->io2parr).ptr;
          sum += s[k].bytes;
          const long long int idx = s[k].offset - start_p;
          memcpy(buf + idx, s[k].buf, s[k].bytes);
        }
        ierr = 0;
        if (sum != count) ierr = 1;
        MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX,
                      nek_fh->aggcomm.c);
        if (ierr) {
          return NEKIO_RDERR_USRBUF;
        }
      }
    }
  }
  return 0;
}

int NEK_File_write(void *handle, void *buf, long long int count,
                   long long int offset) {
  int ierr;
  int ferr;
  nekfh *nek_fh = (nekfh *)handle;
  struct comm c = nek_fh->comm;
  MPI_Comm comm = c.c;
  ierr = 0;

  ierr = (count < 0);
  MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
  if (ierr) {
    return NEKIO_WRERR_COUNT;
  }

  if (nek_fh->mpiio) {
    // Use MPIIO
    ierr = (nek_fh->mode == MPI_MODE_RDONLY);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_WRERR_RDONL;
    }

    ierr = (*(nek_fh->mpifh) == NULL);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_WRERR_NONFL;
    }

    ferr = MPI_File_set_view(*(nek_fh->mpifh), offset, MPI_BYTE, MPI_BYTE,
                             "native", MPI_INFO_NULL);
    ierr = (ferr != MPI_SUCCESS);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_WRERR_MPISV;
    }

    ferr = MPI_File_write_all(*(nek_fh->mpifh), buf, count, MPI_BYTE,
                              MPI_STATUS_IGNORE);
    ierr = (ferr != MPI_SUCCESS);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_WRERR_MPIWA;
    }

  } else {
    // TODO: parallel write
    return NEKIO_WRERR_NSUPP;
  }
  return 0;
}

void NEK_File_close(void *handle, int *ierr) {
  nekfh *nek_fh = (nekfh *)handle;
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
      *ierr = 1;
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
      *ierr = 1;
      return;
    }
  }

  crystal_free(&(nek_fh->cr));
  array_free(&(nek_fh->tarr));
  array_free(&(nek_fh->io2parr));
  comm_free(&(nek_fh->comm));
  *ierr = 0;
}

int NEK_File_EOF(void *handle, int *ifeof) {
  nekfh *nek_fh = (nekfh *)handle;
  int eof_g, eof_p;
  int ierr = 0;
  MPI_Offset curr, end;
  struct comm c = nek_fh->comm;
  MPI_Comm comm = c.c;
  eof_p = 0;
  eof_g = 0;

  if (nek_fh->mpiio) {
    ierr = (*(nek_fh->mpifh) == NULL);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_RDERR_NONFL;
    }

    ierr = MPI_File_get_position(*(nek_fh->mpifh), &curr);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_RDERR_GTPOS;
    }

    ierr = MPI_File_seek(*(nek_fh->mpifh), 0, MPI_SEEK_END);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_RDERR_FSEEK;
    }

    ierr = MPI_File_get_position(*(nek_fh->mpifh), &end);
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_RDERR_GTPOS;
    }

    ierr = MPI_File_seek(*(nek_fh->mpifh), curr,
                         MPI_SEEK_SET);  // Put individual file pointer back
    MPI_Allreduce(MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm);
    if (ierr) {
      return NEKIO_RDERR_FSEEK;
    }

    eof_p = (curr == end);
    MPI_Allreduce(&eof_p, &eof_g, 1, MPI_INT, MPI_MAX, comm);
    *ifeof = eof_g;
  } else {
    // Test read to active eof indicator, safe operation since every
    // fread in nekio is preceded by fseek
    char *tmp_buf = malloc(sizeof(char));
    fread(tmp_buf, 1, 1, nek_fh->file);
    eof_p = feof(nek_fh->file);
    free(tmp_buf);
    MPI_Allreduce(&eof_p, &eof_g, 1, MPI_INT, MPI_MAX, comm);
    *ifeof = eof_g;
  }

  return 0;
}

// Wrapper for Fortran
void fNEK_File_open(const int *fcomm, char *filename, int *amode, int *ifmpiio,
                    int *cb_nodes, int *handle, int *ierr, int nlen) {
  *ierr = 1;
  comm_ext c = MPI_Comm_f2c(*fcomm);
  if (handle_n == handle_max) {
    handle_max += handle_max / 2 + 1;
    fhandle_arr = trealloc(nekfh *, fhandle_arr, handle_max);
  }
  fhandle_arr[handle_n] = (nekfh *)tmalloc(nekfh, 1);

  char fname[MAX_NAME + 1];
  {
    int i;
    strncpy(fname, filename, MAX_NAME);
    for (i = nlen - 1; i > 0; i--)
      if (fname[i] != ' ') break;
    fname[i + 1] = '\0';
  }
  *ierr = NEK_File_open(c, fhandle_arr[handle_n], fname, *amode, *ifmpiio,
                        *cb_nodes);
  *handle = handle_n++;
}

void fNEK_File_read(int *handle, long long int *count, long long int *offset,
                    void *buf, int *ierr) {
  nekfh *fh = fhandle_arr[*handle];
  *ierr = NEK_File_read(fh, buf, *count, *offset);
}

void fNEK_file_write(int *handle, long long int *count, long long int *offset,
                     void *buf, int *ierr) {
  nekfh *fh = fhandle_arr[*handle];
  *ierr = NEK_File_write(fh, buf, *count, *offset);
}

void fNEK_File_close(int *handle, int *ierr) {
  NEK_File_close(fhandle_arr[*handle], ierr);
  free(fhandle_arr[*handle]);
  fhandle_arr[*handle] = 0;
}

int fNEK_File_EOF(int *handle, int *ierr) {
  int ifeof;
  nekfh *fh = fhandle_arr[*handle];
  ifeof = 0;
  *ierr = NEK_File_EOF(fh, &ifeof);

  return ifeof;
}
