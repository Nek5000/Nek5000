#include "gslib.h"

#ifdef PARMETIS

#include "parmetis.h"
#include "defs.h"

#define MAXNV 8 /* maximum number of vertices per element */
typedef struct { long long vtx[MAXNV]; long long eid; int proc; } vtx_data;


int parMETIS_partMesh(long long *elo, long long *vlo, int *nelo, 
                      long long *el , long long *vl , const int nel,
                      const int nv, comm_ext ce)
{
  struct comm comm;
  struct crystal cr;

  comm_init(&comm,ce);
  int np = comm.np;
  int myid = comm.id;

  int i, j;

  comm_barrier(&comm);
  if(myid == 0) printf("Running parMETIS ... ");
  fflush(stdout);

  if (sizeof(idx_t) != sizeof(long long)){
    if (myid == 0) printf("ERROR: invalid sizeof(idx_t)!\n");
    return 1;
  }

  if (nv != 4 && nv != 8){
    if (myid == 0) printf("ERROR: nv is %d but only 4 and 8 are supported!\n", nv);
    return 1;
  }

  long long *nelarray = (long long*) malloc(np*sizeof(long long));
  long long buf = nel;
  MPI_Allgather(&buf, 1, MPI_LONG_LONG_INT, nelarray, 1, MPI_LONG_LONG_INT, comm.c);
  idx_t *elmdist = (idx_t*) malloc((np+1)*sizeof(idx_t));
  elmdist[0] = 0;
  for (i=0; i<np; ++i) elmdist[i+1] = elmdist[i] + (idx_t)nelarray[i];
  free(nelarray); 

  idx_t *evlptr = (idx_t*) malloc((nel+1)*sizeof(idx_t));
  evlptr[0] = 0;
  for (i=0; i<nel; ++i) evlptr[i+1] = evlptr[i] + nv;
  idx_t nelsm = elmdist[myid+1]- elmdist[myid];
  //printf("nelsm %d %d %d\n", elmdist[myid+1], elmdist[myid], nelsm);
  evlptr[nelsm]--;

  idx_t *elmwgt = NULL;
  idx_t wgtflag = 0; // no weights 
  idx_t numflag = 0; // all IDs start from one (Fortran-style)
  idx_t ncon = 1;

  idx_t ncommonnodes = 2; // face vertices
  if (nv == 8) ncommonnodes = 4;
  idx_t nparts = np; 
  
  real_t *tpwgts = (real_t*) malloc(ncon*nparts*sizeof(real_t));
  for (i=0; i<ncon*nparts; ++i){
    tpwgts[i] = 1./(real_t)nparts;
  }

  real_t ubvec = 1.02; //UNBALANCE_FRACTION;

  idx_t options[10];
  options[0] = 1;
  options[PMV3_OPTION_DBGLVL] = 0;
  options[PMV3_OPTION_SEED] = 0;
  idx_t edgecut = 0;

  idx_t *part = (idx_t*) malloc(nel*sizeof(idx_t));

  int ierr;

  comm_barrier(&comm); double t0 = comm_time();
  ierr = ParMETIS_V3_PartMeshKway(elmdist,
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

  if (ierr == METIS_OK) ierr = 0;

  free(elmdist);
  free(evlptr);
  free(tpwgts);

  // redistribute data for target proc 
  crystal_init(&cr,&comm);

  struct array A; 
  vtx_data *row;
  array_init(vtx_data, &A, nel), A.n = nel;
  for(row = A.ptr, i = 0; i < A.n; ++i) {
    for(j = 0; j < nv; ++j) row[i].vtx[j] = vl[i*nv+j];  
    row[i].eid = el[i];
    row[i].proc = part[i];
    //printf("send nid=%d, eid=%d proc=%d\n", myid, row[i].eid, row[i].proc);
  }

  free(part);
  sarray_transfer(vtx_data, &A, proc, 0, &cr);

  if (*nelo < A.n){
    printf("ERROR: nelo too small to hold resulting parition!\n");
    return 1;
  }
  *nelo = A.n;

  for(row = A.ptr, i = 0; i < *nelo; ++i) {
    for(j = 0; j < nv; ++j) vlo[i*nv+j] = row[i].vtx[j];  
    elo[i] = row[i].eid;
    //printf("recv nid=%d, elo=%d\n", myid, elo[i]);
  }

  comm_barrier(&comm); double time = comm_time() - t0;
  if(myid == 0) printf("%lf sec\n", time);

  array_free(&A);
  crystal_free(&cr);
  comm_free(&comm);
                                  
  return ierr;
}

#define fparMETIS_partMesh FORTRAN_UNPREFIXED(fparmetis_partmesh,FPARMETIS_PARTMESH)
void fparMETIS_partMesh(long long *egl, long long *vl, int *negl,
                        long long *eglcon, long long *vlcon, int *neglcon,
                        int *nve, int *comm, int *err)
{
  *err = 1;

#if defined(MPI)
  comm_ext c = MPI_Comm_f2c(*comm);
#else
  comm_ext c = 0;
#endif

  *err = parMETIS_partMesh(egl, vl, negl,
                           eglcon, vlcon, *neglcon,
                           *nve, c);
}
#endif

void printPartStat(long long *vtx, int nel, int nv, comm_ext ce)
{
  struct comm comm;
  comm_init(&comm,ce);

  int np = comm.np;
  int id = comm.id;

  int i,j;

  int numPoints = nel*nv;
  long long *data= (long long*) malloc(numPoints*sizeof(long long));
  for(i = 0; i < numPoints; i++) data[i] = vtx[i]; 

  struct gs_data *gsh;
  gsh = gs_setup(data, numPoints, &comm, 0, gs_pairwise, 0);

  buffer buf;
  buffer_init(&buf, 1024);

  int neighborsCount = 0;
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

  int ncMax = neighborsCount;
  int ncMin = neighborsCount;
  int ncSum = neighborsCount;

  int b;
  comm_allreduce(&comm, gs_int, gs_max, &ncMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &ncMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &ncSum, 1, &b);

  int nelMax = nel, nelMin = nel;
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);
  double imb = (double)nelMax/nelMin;

  if(id == 0) {
    printf(" Max neighbors: %d",ncMax);
    printf(" | Min neighbors: %d",ncMin);
    printf(" | Avg neighbors: %lf\n",(double)ncSum/np);
    printf(" Max elements: %d",nelMax);
    printf(" | Min elements: %d",nelMin);
    printf(" | Balance: %lf\n",imb);
  }

  fflush(stdout);
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
