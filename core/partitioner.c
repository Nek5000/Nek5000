#include "gslib.h"

#if defined(PARRSB)
#include "parRSB.h"
#endif

#define MAXNV 8 /* maximum number of vertices per element */
typedef struct {long long vtx[MAXNV]; long long eid; int proc;} edata;


#ifdef PARMETIS

#include "parmetis.h"
#include "defs.h"

int parMETIS_partMesh(int *part, long long *vl, int nel, int nv, int *opt, comm_ext ce)
{
  int i, j;
  int ierrm;
  double time, time0;

  MPI_Comm comms;
  struct comm comm;
  int color;
  int ibuf;

  struct crystal cr;
  struct array A; 
  edata *row;

  long long nell;
  long long *nelarray;
  idx_t *elmdist;
  idx_t *evlptr;
  idx_t *part_;
  real_t *tpwgts;
  idx_t edgecut;
  real_t ubvec;
  idx_t *elmwgt;
  idx_t wgtflag;
  idx_t numflag;
  idx_t ncon;
  idx_t ncommonnodes;
  idx_t nparts;
  idx_t nelsm;
  idx_t options[10];

  ierrm = METIS_OK;
  nell = nel;
  edgecut = 0;
  wgtflag = 0;
  numflag = 0;
  ncon = 1;
  ubvec = 1.02;
  elmwgt = NULL; /* no weights */
  ncommonnodes = 2;

  part_ = (idx_t*) malloc(nel*sizeof(idx_t));

  if (sizeof(idx_t) != sizeof(long long)){
    printf("ERROR: invalid sizeof(idx_t)!\n");
    goto err;
  }
  if (nv != 4 && nv != 8){
    printf("ERROR: nv is %d but only 4 and 8 are supported!\n", nv);
    goto err;
  }

  color = MPI_UNDEFINED;
  if (nel > 0) color = 1;
  MPI_Comm_split(ce, color, 0, &comms);
  if (color == MPI_UNDEFINED)
    goto end;

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
      part_[i] = comm.id;

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
                                   part_,
                                   &comm.c);

  time = comm_time() - time0;
  if (comm.id == 0) 
    printf("%lf sec\n", time), fflush(stdout);

  for (i=0; i<nel; ++i) 
    part[i] = part_[i];

  free(elmdist);
  free(evlptr);
  free(tpwgts);
  MPI_Comm_free(&comms);
  comm_free(&comm);

end: 
  comm_init(&comm,ce);
  comm_allreduce(&comm, gs_int, gs_min, &ierrm, 1, &ibuf); 
  if (ierrm != METIS_OK) goto err;
  return 0;
                                 
err:
  return 1;
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

#define fpartMesh FORTRAN_UNPREFIXED(fpartmesh,FPARTMESH)
void fpartMesh(long long *el, long long *vl, const int *lelt, int *nell, 
               const int *nve, comm_ext *fcomm, int *rtval) 
{
  struct comm comm;
  struct crystal cr;
  struct array eList;
  edata *data;

  int nel, nv;
  int e, n; 
  int count, ierr, ibuf;
  int *part;
  int opt[3];

  nel  = *nell;
  nv   = *nve;
  comm_init(&comm, MPI_Comm_f2c(*fcomm));

  part = (int*) malloc(nel * sizeof(int));

  ierr = 1;
#if defined(PARRSB)
  opt[0] = 1;
  opt[1] = 2; /* verbosity */
  opt[2] = 0;
  ierr = parRSB_partMesh(part, vl, nel, nv, opt, comm.c);
#elif defined(PARMETIS)
  opt[0] = 1;
  opt[1] = 0; /* verbosity */
  opt[2] = comm.np;
  ierr = parMETIS_partMesh(part, vl, nel, nv, opt, comm.c);
#endif
  if (ierr != 0) goto err; 

  /* redistribute data */
  array_init(edata, &eList, nel), eList.n = nel;
  for(data = eList.ptr, e = 0; e < nel; ++e) {
    data[e].proc = part[e];
    data[e].eid  = el[e];
    for(n = 0; n < nv; ++n) {
      data[e].vtx[n] = vl[e*nv + n];
    }
  }
  free(part);

  crystal_init(&cr, &comm);
  sarray_transfer(edata, &eList, proc, 0, &cr);
  crystal_free(&cr);

  nel = eList.n;

  count = 0;
  if (nel > *lelt) count = 1;
  comm_allreduce(&comm, gs_int, gs_add, &count, 1, &ibuf);
  if (count > 0) {
    if (comm.id == 0)
      printf("ERROR: resulting parition requires lelt=%d!\n", nel);
    goto err;
  }

  for(data = eList.ptr, e = 0; e < nel; ++e) {
    el[e] = data[e].eid;
    for(n = 0; n < nv; ++n) {
      vl[e*nv + n] = data[e].vtx[n];
    }
  }
  
  array_free(&eList);
  comm_free(&comm);

  *nell = nel;
  *rtval = 0;
  return;

err:
  fflush(stdout);
  *rtval = 1;
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
