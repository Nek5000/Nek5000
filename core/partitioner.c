#include "gslib.h"

#ifdef PARMETIS

#include "parmetis.h"
#include "defs.h"

#define MAXNV 8 /* maximum number of vertices per element */
typedef struct {long long vtx[MAXNV]; long long eid; int proc;} vtx_data;


int parMETIS_partMesh(long long *elo, long long *vlo, int *nelo, 
                      long long *el , long long *vl , const int nel,
                      const int nv, int *opt, comm_ext ce)
{
  int i, j;
  int ierrm = METIS_OK;
  double time, time0;

  MPI_Comm comms;
  struct comm comm;
  int color = MPI_UNDEFINED;
  int ibuf;

  struct crystal cr;
  struct array A; 
  vtx_data *row;

  long long nell = nel;
  long long *nelarray;
  idx_t *elmdist;
  idx_t *evlptr;
  real_t *tpwgts;
  idx_t edgecut = 0;
  real_t ubvec = 1.02;
  idx_t *elmwgt = NULL; /* no weights */
  idx_t *part;
  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon = 1;
  idx_t ncommonnodes = 2;
  idx_t nparts;
  idx_t nelsm;
  idx_t options[10];

  part = (idx_t*) malloc(nel*sizeof(idx_t));

  if (sizeof(idx_t) != sizeof(long long)){
    printf("ERROR: invalid sizeof(idx_t)!\n");
    goto err;
  }
  if (nv != 4 && nv != 8){
    printf("ERROR: nv is %d but only 4 and 8 are supported!\n", nv);
    goto err;
  }

  if (nel > 0) color = 1;
  MPI_Comm_split(ce, color, 0, &comms);
  if (color == MPI_UNDEFINED)
    goto distribute;

  comm_init(&comm,comms);
  if (comm.id == 0) 
    printf("Running parMETIS ... "), fflush(stdout);

  nelarray = (long long*) malloc(comm.np*sizeof(long long));
  MPI_Allgather(&nell, 1, MPI_LONG_LONG_INT, nelarray, 1, MPI_LONG_LONG_INT, comm.c);
  elmdist = (idx_t*) malloc((comm.np+1)*sizeof(idx_t));
  elmdist[0] = 0;
  for (i=0; i<comm.np; ++i)
    elmdist[i+1] = elmdist[i] + (idx_t)nelarray[i];
  free(nelarray); 

  evlptr = (idx_t*) malloc((nel+1)*sizeof(idx_t));
  evlptr[0] = 0;
  for (i=0; i<nel; ++i) 
    evlptr[i+1] = evlptr[i] + nv;
  nelsm = elmdist[comm.id+1] - elmdist[comm.id];
  evlptr[nelsm]--;

  if (nv == 8) ncommonnodes = 4;
  nparts = comm.np;

  options[0] = 1;
  options[PMV3_OPTION_DBGLVL] = 0;
  options[PMV3_OPTION_SEED]   = 0;
  if (opt[0] != 0) {
    options[PMV3_OPTION_DBGLVL] = opt[1]; 
    if (opt[2] != 0) {
      options[3] = PARMETIS_PSR_UNCOUPLED;
      nparts = opt[2];
    }
  }  

  tpwgts = (real_t*) malloc(ncon*nparts*sizeof(real_t));
  for (i=0; i<ncon*nparts; ++i)
    tpwgts[i] = 1./(real_t)nparts;

  if (options[3] == PARMETIS_PSR_UNCOUPLED)
    for (i=0; i<nel; ++i) 
      part[i] = comm.id;

  comm_barrier(&comm); 
  time0 = comm_time();
  ierrm = ParMETIS_V3_PartMeshKway(elmdist,
                                   evlptr, 
                                   (idx_t*)vl,
                                   elmwgt,
                                   &wgtflag,
                                   &numflag,
                                   &ncon,
                                   &ncommonnodes,
                                   &nparts,
                                   tpwgts,
                                   &ubvec,
                                   options,
                                   &edgecut,
                                   part,
                                   &comm.c);

  time = comm_time() - time0;
  if (comm.id == 0) 
    printf("%lf sec\n", time), fflush(stdout);

  free(elmdist);
  free(evlptr);
  free(tpwgts);
  MPI_Comm_free(&comms);
  comm_free(&comm);

distribute: 
  comm_init(&comm,ce);
  comm_allreduce(&comm, gs_int, gs_min, &ierrm, 1, &ibuf); 
  if (ierrm != METIS_OK) goto err;

  array_init(vtx_data, &A, nel), A.n = nel;
  for(row = A.ptr, i = 0; i < A.n; ++i) {
    for(j = 0; j < nv; ++j) row[i].vtx[j] = vl[i*nv+j];  
    row[i].eid = el[i];
    row[i].proc = part[i];
  }
  free(part);

  crystal_init(&cr,&comm);
  sarray_transfer(vtx_data, &A, proc, 0, &cr);

  ierrm = 0;
  if (A.n > *nelo) ierrm = 1;
  comm_allreduce(&comm, gs_int, gs_add, &ierrm, 1, &ibuf);

  *nelo = A.n;

  if (ierrm > 0) {
    if (comm.id == 0)
      printf("ERROR: resulting parition requires lelt=%d!\n", *nelo);
    goto err;
  }

  for (row = A.ptr, i = 0; i < *nelo; ++i) {
    for (j = 0; j < nv; ++j) 
      vlo[i*nv+j] = row[i].vtx[j];  
    elo[i] = row[i].eid;
  }

  array_free(&A);
  crystal_free(&cr);
  comm_free(&comm);

  return 0;
                                 
err:
  fflush(stdout);
  return 1;
}

#define fparMETIS_partMesh FORTRAN_UNPREFIXED(fparmetis_partmesh,FPARMETIS_PARTMESH)
void fparMETIS_partMesh(long long *egl, long long *vl, int *negl,
                        long long *eglcon, long long *vlcon, int *neglcon,
                        int *nve, int *opt, int *comm, int *err)
{
  *err = 1;
  comm_ext c;

#if defined(MPI)
  c = MPI_Comm_f2c(*comm);
#else
  c = 0;
#endif

  *err = parMETIS_partMesh(egl, vl, negl,
                           eglcon, vlcon, *neglcon,
                           *nve, opt, c);
}
#endif

void printPartStat(long long *vtx, int nel, int nv, comm_ext ce)
{
  int i,j;
  int np, id;
  int numPoints = nel*nv;
  struct gs_data *gsh;
  buffer buf;
  long long *data;
  int neighborsCount = 0;
  int nelMin, nelMax;
  int ncMin, ncMax, ncSum;
  int b;
  struct comm comm;

  comm_init(&comm,ce);
  np = comm.np;
  id = comm.id;

  data= (long long*) malloc(numPoints*sizeof(long long));
  for(i = 0; i < numPoints; i++) data[i] = vtx[i]; 

  gsh = gs_setup(data, numPoints, &comm, 0, gs_pairwise, 0);

  buffer_init(&buf, 1024);

  for(i = 0; i < np; i++) {
    if(i != id) {
      for(j = 0; j < numPoints; j++) {
        data[j] = -1;
      }
    } else {
      for(j = 0; j < numPoints; j++) {
        data[j] = id + 1;
      }
    }

    gs(data, gs_long, gs_max, 0, gsh, &buf);

    for(j = 0; j < numPoints; j++) {
      if(data[j] > 0) {
        neighborsCount++;
        break;
      }
    }
  }

  gs_free(gsh);
  free(data);
  buffer_free(&buf);

  ncMax = neighborsCount;
  ncMin = neighborsCount;
  ncSum = neighborsCount;
  nelMax = nel;
  nelMin = nel;

  comm_allreduce(&comm, gs_int, gs_max, &ncMax , 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &ncMin , 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &ncSum , 1, &b);
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);

  if (id == 0) {
    printf(
      " Max neighbors: %d | Min neighbors: %d | Avg neighbors: %lf\n",
      ncMax, ncMin, (double)ncSum/np);
    printf(
      " Max elements: %d | Min elements: %d | Balance: %lf\n",
      nelMax, nelMin, (double)nelMax/nelMin);
    fflush(stdout);
  }

  comm_free(&comm);
}

#define fprintPartStat FORTRAN_UNPREFIXED(printpartstat,PRINTPARTSTAT)
void fprintPartStat(long long *vtx, int *nel, int *nv, int *comm)
{

#if defined(MPI)
  comm_ext c = MPI_Comm_f2c(*comm);
#else
  comm_ext c = 0;
#endif

  printPartStat(vtx, *nel, *nv, c);
}
