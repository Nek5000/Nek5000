#ifdef MPI
  typedef MPI_Comm mpicomm_t;
#else
  typedef int      mpicomm_t;
#endif

#include "gs_hack.h"

/* communication matrix Q

   v [0:nloc-1] are local dofs
   ve[nloc:n-1] are neighbors' dofs
   
   the action of
     ve = Q v
   is to
    copy     ve[0:nloc-1] = v[0:nloc-1]
    retreive ve[nloc:n-1] from neighbors
 */

typedef struct {
  uint nloc, n;
  /*gs_data *gs;*/
  jl_gs_data *jgs;
} amg_Q;

#ifdef DIAGNOSTICS
static int check_Q(const amg_Q *Q)
{
  uint i;
  if(Q->nloc>Q->n) return 1;
  if(Q->gs) return 0;
  return 1;
}
#endif

/* ve = Q v */
static double apply_Q(real *ve, const amg_Q *Q, const real *v, mpicomm_t comm)
{
  double t0=0,t1=0;
  uint i, nloc=Q->nloc,n=Q->n;
  for(i=0;i<nloc;++i) ve[i]=v[i];
  for(;i<n;++i) ve[i]=0;
# ifdef MPI
# ifdef GS_BARRIER
  MPI_Barrier(comm);
# endif
# ifdef GS_TIMING
  t0 = MPI_Wtime();
# endif
  /*gs_op(ve,GS_OP_ADD,Q->gs);*/
#if !defined(USE_FLOAT)
  jl_gs_op(ve,gs_double,gs_add,0,Q->jgs,0);
#else
  jl_gs_op(ve,gs_float,gs_add,0,Q->jgs,0);
#endif

# ifdef GS_TIMING
  t1 = MPI_Wtime() - t0;
# endif
# endif
  return t1;
}

/* z := alpha y + beta Q^t x
   (x := Q Q^t x   as a side effect)
 */
static double apply_Qt(real *z, real alpha, const real *y,
                     real beta, const amg_Q *Q, real *x, mpicomm_t comm)
{
  double t0=0,t1=0;
  uint i, nloc=Q->nloc;
# ifdef MPI
# ifdef GS_BARRIER
  MPI_Barrier(comm);
# endif
# ifdef GS_TIMING
  t0 = MPI_Wtime();
# endif
  /*gs_op(x,GS_OP_ADD,Q->gs);*/
#if !defined(USE_FLOAT)
  jl_gs_op(x,gs_double,gs_add,1,Q->jgs,0);
#else
  jl_gs_op(x,gs_float,gs_add,1,Q->jgs,0);
#endif
# ifdef GS_TIMING
  t1 = MPI_Wtime() - t0;
# endif
# endif
  for(i=0;i<nloc;++i) *z++ = alpha*(*y++) + beta*(*x++);
  return t1;
}

/* sparse matrix */
typedef struct {
  uint rn, cn, *row_off, *col;
  real *pr;
} amg_mat;

#ifdef DIAGNOSTICS
static int check_M(const amg_mat *M)
{
  uint i,n;
  for(i=0;i<M->rn;++i) if(M->row_off[i+1]<M->row_off[i]) return 1;
  n = M->row_off[M->rn];
  for(i=0;i<n;++i) {
    if(M->col[i]>=M->cn) return 1;
    if(isnan(M->pr[i])) return 1;
  }
  return 0;
}
#endif

/* z = alpha y + beta M x */
static double apply_M(real *z, real alpha, const real *y,
                      real beta, const amg_mat *M, const real *x)
{
  uint i,rn=M->rn;
  const uint *row_off = M->row_off, *col = M->col;
  const real *pr = M->pr;
# ifdef GS_TIMING
    double time0 = MPI_Wtime();
# endif
  for(i=0;i<rn;++i) {
    uint j,je; real t = 0;
    for(j=row_off[i],je=row_off[i+1]; j<je; ++j)
      t += (*pr++) * x[*col++];
    *z++ = alpha*(*y++) + beta*t;
  }
# ifdef GS_TIMING
    return MPI_Wtime()-time0;
# else
    return 0;
# endif
}

/* z = M^t x */
static double apply_Mt(real *z, const amg_mat *M, const real *x)
{
  uint i,rn=M->rn,cn=M->cn;
  const uint *row_off = M->row_off, *col = M->col;
  const real *pr = M->pr;
# ifdef GS_TIMING
    double time0 = MPI_Wtime();
# endif

  for(i=0;i<cn;++i) z[i]=0;
  for(i=0;i<rn;++i) {
    uint j,je; real xi = *x++;
    for(j=row_off[i],je=row_off[i+1]; j<je; ++j)
      z[*col++] += (*pr++) * xi;
  }
# ifdef GS_TIMING
    return MPI_Wtime()-time0;
# else
    return 0;
# endif
}

typedef struct {
#ifdef MPI
  MPI_Comm comm;
#endif
  uint pid, np;

  real tn;
  uint un, cn;
  sint *perm;
  jl_gs_data *gs_top;
  
  unsigned levels;
  unsigned *cheb_m;
  real     *cheb_rho;
  unsigned null_space;

  uint *lvl_offset;
  real *Dff;
  amg_Q *Q_W, *Q_AfP, *Q_Aff;
  amg_mat *W, *AfP, *Aff;
  
  real *b, *x, *c, *c_old, *r, *buf;
  
  double *timing; uint timing_n;
} amg_data;

#ifdef DIAGNOSTICS
static void check_amg_data(const amg_data *data)
{
  uint i,n,lvl;
  if(isnan(data->tn)) printf("AMG error: tn not normal (%u)\n",data->pid);
  if(data->cn>=data->un)
    printf("AMG error: user -> condensed permutation wrongly sized (%u)\n",
      data->pid);
  for(i=0;i<data->un;++i)
    if(data->perm[i] < -1 || data->perm[i] >= (sint)data->cn)
      printf("AMG error: user -> condensed permutation not consistent (%u)\n",
        data->pid);
  for(lvl=0;lvl<data->levels-1;++lvl) {
    uint fn = data->lvl_offset[lvl+1] - data->lvl_offset[lvl];
    uint cn = data->lvl_offset[data->levels] - data->lvl_offset[lvl+1];
    if(data->W[lvl].rn != fn || data->Q_W[lvl].nloc != cn ||
       data->W[lvl].cn != data->Q_W[lvl].n)
      printf("AMG error: W %u matrix sized wrongly (%u)\n",lvl,data->pid);
    if(data->AfP[lvl].rn != fn || data->Q_AfP[lvl].nloc != cn ||
       data->AfP[lvl].cn != data->Q_AfP[lvl].n)
      printf("AMG error: AfP %u matrix sized wrongly (%u)\n",lvl,data->pid);
    if(data->Aff[lvl].rn != fn || data->Q_Aff[lvl].nloc != fn ||
       data->Aff[lvl].cn != data->Q_Aff[lvl].n)
      printf("AMG error: Aff %u matrix sized wrongly (%u)\n",lvl,data->pid);
  }
  n = 3*(data->levels-1);
  for(i=0;i<n;++i) {
    if(check_Q(&data->Q_W[i]))
      printf("AMG error: communication (%u) data not consistent (%u)\n",
        i,data->pid);
    if(check_M(&data->W[i]))
      printf("AMG error: matrix (%u) data not consistent (%u)\n",i,data->pid);
  }
  for(i=0;i<data->lvl_offset[data->levels];++i)
    if(isnan(data->Dff[i]))
      printf("AMG error: NAN in Dff matrix %u (%u)\n",i,data->pid);
  for(lvl=data->levels-1;lvl;) {
    unsigned ci,m; real alpha, gamma, beta;
    --lvl;
    m = data->cheb_m[lvl];
    if(m>1) {
      alpha = data->cheb_rho[lvl]/2, alpha *= alpha;
      gamma = 2*alpha/(1-2*alpha), beta = 1 + gamma;
      if(isnan(beta))
        printf("AMG error: Chebyshev coefficient %u NAN (%u)\n",lvl,data->pid);
    }
    for(ci=3;ci<=m;++ci) {
      gamma = alpha*beta/(1-alpha*beta), beta = 1 + gamma;
      if(isnan(beta) || isnan(gamma))
        printf("AMG error: Chebyshev coefficient %u NAN (%u)\n",lvl,data->pid);
    }
  }
}
#endif

void amg_stats(amg_data *data)
{
#ifdef GS_TIMING
  unsigned lvl, lm1 = data->levels - 1;
  double *global, *local;
  if(data->pid==0) global = tmalloc(double,9*lm1*data->np);
  MPI_Gather(data->timing,6*lm1,MPI_DOUBLE,
             global,      6*lm1,MPI_DOUBLE,
             0,data->comm);
  if(data->pid==0) {
#   define DOPRINT(name,i) do { \
      uint p; \
      double *timing = global; \
      for(p=0;p<data->np;++p,timing+=6*lm1) { \
        double *t = timing; \
        printf("%9d",(unsigned)p); \
        for(lvl=0;lvl<lm1;++lvl,t+=6) \
          printf(" %0.3e", t[i]/data->timing_n); \
        puts(" AMG " name); \
      } \
    } while(0)
    DOPRINT("Wt    work",1);
    DOPRINT("W,AfP work",3);
    DOPRINT("Aff   work",5);
    DOPRINT("Wt    comm",0);
    DOPRINT("W,AfP comm",2);
    DOPRINT("Aff   comm",4);
#   undef DOPRINT
  }
  local = tmalloc(double,9*lm1);
  for(lvl=0;lvl<9*lm1;++lvl) local[lvl]=0;
  for(lvl=0;lvl<lm1;++lvl) {
    if(data->Q_W[lvl].gs)
      gs_data_stats(local + 9*lvl + 0, data->Q_W[lvl].gs);
    if(data->Q_AfP[lvl].gs)
      gs_data_stats(local + 9*lvl + 3, data->Q_AfP[lvl].gs);
    if(data->Q_Aff[lvl].gs)
      gs_data_stats(local + 9*lvl + 6, data->Q_Aff[lvl].gs);
  }
  MPI_Gather(local ,9*lm1,MPI_DOUBLE,
             global,9*lm1,MPI_DOUBLE,
             0,data->comm);
  if(data->pid==0) {
#   define DOPRINT(name,i,fmt,type) do { \
      uint p; \
      double *gsdata = global; \
      for(p=0;p<data->np;++p,gsdata+=9*lm1) { \
        double *t = gsdata; \
        printf("%9d",(unsigned)p); \
        for(lvl=0;lvl<lm1;++lvl,t+=9) \
          printf(" " fmt, (type)t[i]); \
        puts(" AMG " name); \
      } \
    } while(0)
    DOPRINT("Wt    msgs",0,"%9d"  ,unsigned);
    DOPRINT("Wt    avwd",1,"%0.3e",double  );
    DOPRINT("Wt    mxwd",2,"%9d"  ,unsigned);
    DOPRINT("W,AfP msgs",3,"%9d"  ,unsigned);
    DOPRINT("W,AfP avwd",4,"%0.3e",double  );
    DOPRINT("W,AfP mxwd",5,"%9d"  ,unsigned);
    DOPRINT("Aff   msgs",6,"%9d"  ,unsigned);
    DOPRINT("Aff   avwd",7,"%0.3e",double  );
    DOPRINT("Aff   mxwd",8,"%9d"  ,unsigned);
#   undef DOPRINT
  }
  free(local);
  if(data->pid==0) free(global);
#endif
}

static void amg_exec(amg_data *data)
{
  unsigned lvl, levels = data->levels;
  const uint *off = data->lvl_offset;
  uint off_bot;
  real *b = data->b, *x = data->x;
  real *c = data->c, *c_old = data->c_old, *r = data->r;
  double *timing = data->timing;
#ifdef MPI
  mpicomm_t comm = data->comm;
#else
  mpicomm_t comm = 0;
#endif  
  for(lvl=0;lvl<levels-1;++lvl,timing+=6) {
    real *b_l = b+off[lvl], *b_lp1 = b+off[lvl+1];
    /* b_{l+1} += W^t b_l */
    timing[1] += apply_Mt(data->buf, &data->W[lvl],b_l);
    timing[0] += apply_Qt(b_lp1, 1,b_lp1, 1,&data->Q_W[lvl],data->buf,comm);
  }
  if(off[levels]-(off_bot=off[levels-1]))
    x[off_bot] = data->Dff[off_bot]*b[off_bot];
  for(lvl=levels-1;lvl;) {
    real *b_l, *x_l, *x_lp1, *d_l; uint i,n;
    unsigned ci,m; real alpha, gamma, beta;
    --lvl;
    timing = data->timing + lvl*6;
    b_l = b+off[lvl], x_l = x+off[lvl], x_lp1 = x+off[lvl+1];
    d_l = data->Dff+off[lvl];
    n = off[lvl+1]-off[lvl];
    m = data->cheb_m[lvl];
    /* buf = Q x_{l+1} */
    timing[2] += apply_Q(data->buf, &data->Q_AfP[lvl], x_lp1, comm);
    /* x_l = W x_{l+1} */
    /* careful: x_l not initialized; 0*x_l could contain nans */
    timing[3] += apply_M(x_l, 0,b_l, 1,&data->W[lvl],data->buf);
    /* b_l -= AfP x_{l+1} */
    timing[3] += apply_M(b_l, 1,b_l, -1,&data->AfP[lvl],data->buf);
    /* c_1 = D b_l */
    for(i=0;i<n;++i) c[i] = d_l[i] * b_l[i];
    if(m>1) {
      alpha = data->cheb_rho[lvl]/2, alpha *= alpha;
      gamma = 2*alpha/(1-2*alpha), beta = 1 + gamma;
      /* r_1 = b_l - Aff c_1 */
      timing[4] += apply_Q(data->buf, &data->Q_Aff[lvl],c, comm);
      timing[5] += apply_M(r, 1,b_l, -1,&data->Aff[lvl],data->buf);
      /* c_2 = (1+gamma) (c_1 + D r_1) */
      for(i=0;i<n;++i) c_old[i] = beta*(c[i]+d_l[i]*r[i]);
      { real *temp = c; c = c_old; c_old = temp; }
    }
    for(ci=3;ci<=m;++ci) {
      gamma = alpha*beta/(1-alpha*beta), beta = 1 + gamma;
      /* r_i = b_l - Aff c_i */
      timing[4] += apply_Q(data->buf, &data->Q_Aff[lvl],c, comm);
      timing[5] += apply_M(r, 1,b_l, -1,&data->Aff[lvl],data->buf);
      /* c_{i+1} = (1+gamma) (c_i + D r_i) - gamma c_{i-1} */
      for(i=0;i<n;++i) c_old[i] = beta*(c[i]+d_l[i]*r[i]) - gamma*c_old[i];
      { real *temp = c; c = c_old; c_old = temp; }
    }
    for(i=0;i<n;++i) x_l[i] += c[i];
  }
  data->timing_n++;
}

void amg_solve(real *x, amg_data *data, const real *b)
{
  uint i, un=data->un, cn=data->cn, ln=data->lvl_offset[data->levels];
  real *cb = data->b, *cx = data->x;
  for(i=0;i<cn;++i) cb[i]=0;
  for(i=0;i<un;++i) {
    sint p = data->perm[i];
    if(p!=-1) cb[p] += b[i];
  }
  jl_gs_op(cb,gs_double,gs_add,0,data->gs_top,0);
  amg_exec(data);
  if(data->null_space) {
    real avg = 0, sum;
    for(i=0;i<ln;++i) avg += cx[i];
#ifdef MPI
    MPI_Allreduce(&avg,&sum,1,REAL_MPI,MPI_SUM,data->comm);
#endif
    avg = sum/data->tn;
    for(i=0;i<ln;++i) cx[i] -= avg;
  }
  for(i=ln;i<cn;++i) cx[i]=0;
  jl_gs_op(cx,gs_double,gs_add,0,data->gs_top,0);
  for(i=0;i<un;++i) {
    sint p = data->perm[i];
    x[i] = (p == -1 ? 0 : cx[p]);
  }
}

/* reverse the byte order of the given double
   portable up to 129-byte doubles, standards conforming */
#define N sizeof(double)
static double byteswap(double x)
{
  char t, buf[N];
  memcpy(buf,&x,N);
#define SWAP1(i) if(N>2*(i)+1) t=buf[i],buf[i]=buf[N-1-(i)],buf[N-1-(i)]=t
#define SWAP2(i) SWAP1(i); SWAP1((i)+0x01); SWAP1((i)+0x02); SWAP1((i)+0x03)
#define SWAP3(i) SWAP2(i); SWAP2((i)+0x04); SWAP2((i)+0x08); SWAP2((i)+0x0c)
#define SWAP4(i) SWAP3(i); SWAP3((i)+0x10); SWAP3((i)+0x20); SWAP3((i)+0x30)
  SWAP4(0);
#undef SWAP1
#undef SWAP2
#undef SWAP3
#undef SWAP4
  memcpy(&x,buf,N);
  return x;
}
#undef N

typedef struct {
  FILE *fptr;
  int swap;  
} amg_file;

static void amg_fopen(amg_file *f, const char *filename)
{
  const double magic = 3.14159;
  double t;
  f->fptr = fopen(filename,"r");
  if(f->fptr==0) fail("AMG: could not open %s for reading\n",filename);
  fread(&t,sizeof(double),1,f->fptr);
  if(fabs(t-magic)>0.000001) {
    t=byteswap(t);
    if(fabs(t-magic)>0.000001)
      fail("AMG: magic number for endian test not found in %s\n",filename);
    f->swap = 1;
  } else {
    f->swap = 0;
  }
}

static void amg_fread(double *ptr, size_t n, amg_file *f)
{
  size_t nread = fread(ptr,sizeof(double),n,f->fptr);
  if(!nread && n) fail("AMG: failed reading %u doubles from disk\n",n);
  if(f->swap) {
    size_t i;
    for(i=0;i<n;++i) ptr[i]=byteswap(ptr[i]);
  }
}

static void amg_fclose(amg_file *f)
{
  fclose(f->fptr);
}

static void arrange_cdof_by_level(amg_data *data, tuple_list *cdof, buffer *buf)
{
  uint i; unsigned level;
  
  tuple_list_sort(cdof, cdof_level, buf);
  
  data->lvl_offset = tmalloc(uint, data->levels+1);
  data->lvl_offset[0] = 0; level=1;
  for(i=0;i<cdof->n;++i) {
    sint *t = &cdof->vi[cdof_mi*i];
    sint l=t[cdof_level], proc=t[cdof_proc];
    if(proc) break;
    t[cdof_index]=i;
    if(l>=level) for(;level<=l;++level) data->lvl_offset[level]=i;
  }
  for(;level<=data->levels;++level) data->lvl_offset[level]=i;
  for(;i<cdof->n;++i) cdof->vi[cdof_mi*i+cdof_index]=i;
}

static const int mat_mi=2, mat_ml=1, mat_mr=1;
static const int mat_ridx=0, mat_cidx=1;

static void add_ids(tuple_list *col_id, tuple_list *mat, buffer *buf)
{
  slong last_id = -1; uint last_cidx;
  uint i,n = mat->n;
  sint *mat_vi = mat->vi; slong *mat_vl = mat->vl;
  uint ci=0,cn = col_id->n;
  sint *col_vi = col_id->vi; slong *col_vl = col_id->vl;
  tuple_list_sort(col_id,1,buf);
  tuple_list_sort(mat,mat_mi,buf);
  for(i=0;i<n;++i,mat_vi+=mat_mi,mat_vl+=mat_ml) {
    ulong id = mat_vl[0];
    if(id==last_id) { mat_vi[mat_cidx] = last_cidx; continue; }
    last_id = id;
    while(ci<cn && col_vl[0]<id) ++ci,++col_vi,++col_vl;
    if(ci<cn && col_vl[0]==id)
      mat_vi[mat_cidx] = last_cidx = col_vi[0];
    else {
      uint cidn = col_id->n++;
      if(cidn==col_id->max) {
        tuple_list_grow(col_id);
        col_vi = col_id->vi+ci, col_vl = col_id->vl+ci;
      }
      mat_vi[mat_cidx] = col_id->vi[cidn] = last_cidx = cidn;
      col_id->vl[cidn] = id;
    }
  }
}

static void organize_matrix(amg_Q *Q, amg_mat *M,
                            uint rn, uint cnloc, real cnglob,
                            tuple_list *id, tuple_list *mat,
                            crystal_data *crystal, buffer *buf)
{
  uint ri,i,n, *row_off, *col; real *pr;
  sint *mat_vi; real *mat_vr;
  add_ids(id,mat,buf);
  tuple_list_sort(id,0,buf);
  
  Q->nloc = cnloc;
  Q->n = id->n;
  /*Q->gs = gs_data_setup(id->n,(ulong*)id->vl,1,crystal);*/
  {
#ifdef MPI
    jl_comm_t comm = {crystal->id,crystal->num,crystal->comm};
#else
    jl_comm_t comm = {0,1,0};
#endif
    uint i,ie;
    for(i=cnloc,ie=id->n;i!=ie;++i) id->vl[i] = -id->vl[i];
    Q->jgs = jl_gs_setup(id->vl,id->n,&comm);
    for(i=cnloc,ie=id->n;i!=ie;++i) id->vl[i] = -id->vl[i];
  }

  M->rn = rn;
  M->cn = id->n;
  row_off = M->row_off = tmalloc(uint,(rn+1)+mat->n);
  col     = M->col = M->row_off + (rn+1);
  pr      = M->pr = tmalloc(real,mat->n);
  tuple_list_sort(mat,mat_cidx,buf);
  tuple_list_sort(mat,mat_ridx,buf);
  mat_vi = mat->vi, mat_vr = mat->vr;
  ri = 0, n = mat->n;
  for(i=0;i<n;++i,mat_vi+=mat_mi,mat_vr+=mat_mr) {
    uint ridx = mat_vi[mat_ridx];
    *col++ = mat_vi[mat_cidx], *pr++ = mat_vr[0];
    for(;ri<=ridx;++ri) row_off[ri]=i;
  }
  for(;ri<=rn;++ri) row_off[ri]=n;
}

static void organize_data(amg_data *data, tuple_list *cdof,
                          uint *row_off[3], double *raw_mat[3],
                          crystal_data *crystal, buffer *buf)
{
  uint i; unsigned m,lvl,nl = data->levels;
  uint nloc;
  uint *perm_level;
  real *gln, *gln_in;
  tuple_list mat[3];
  tuple_list C_id, F_id; uint Fnloc, Cnloc; real Fnglob, Cnglob;
  for(i=0;i<cdof->n;++i)
    if(cdof->vi[cdof_mi*i+cdof_proc])
      cdof->vi[cdof_mi*i+cdof_level] = nl;

  arrange_cdof_by_level(data,cdof,buf);

  gln = tmalloc(real, nl);
  buffer_reserve(buf, nl*sizeof(real)); gln_in = buf->ptr;
  for(i=0;i<nl;++i) gln_in[i] = data->lvl_offset[i+1]-data->lvl_offset[i];
#ifdef MPI  
  MPI_Allreduce(gln_in,gln,nl,REAL_MPI,MPI_SUM,data->comm);
#endif
  if(data->pid==0)
    for(i=0;i<nl;++i)
      printf("AMG level %d F-vars: %g\n",i+1,gln[i]);

  nloc = data->lvl_offset[nl];
  data->Q_W = tmalloc(amg_Q,3*(nl-1));
  data->Q_AfP = data->Q_W   + (nl-1);
  data->Q_Aff = data->Q_AfP + (nl-1);
  data->W = tmalloc(amg_mat,3*(nl-1));
  data->AfP = data->W   + (nl-1);
  data->Aff = data->AfP + (nl-1);
  
  data->Dff = tmalloc(real, nloc);
  for(i=0;i<nloc;++i) data->Dff[i] = cdof->vr[i];
  
  tuple_list_init(&C_id,1,1,0);
  tuple_list_init(&F_id,1,1,0);
  for(m=0;m<3;++m) {
    uint size=0;
    for(i=0;i<nloc;++i) {
      uint len = row_off[m][i];
      row_off[m][i] = size;
      size += len;
    }
    row_off[m][i] = size;
    tuple_list_init(&mat[m],mat_mi,mat_ml,mat_mr);
  }

  tuple_list_sort(cdof, cdof_mi, buf);
  tuple_list_sort(cdof, cdof_proc, buf);
  perm_level = tmalloc(uint,nloc);
  buffer_reserve(buf,2*nloc*sizeof(sort_data));
  index_sort((uint*)cdof->vi+cdof_level,nloc,cdof_mi, perm_level, buf->ptr);
  tuple_list_sort(cdof, cdof_index, buf);

  for(lvl=0;lvl<nl-1;++lvl) {
    uint kb=data->lvl_offset[lvl],ke=data->lvl_offset[lvl+1],k;
    for(m=0;m<3;++m) mat[m].n=0;
    for(k=kb;k!=ke;++k) {
      uint i = perm_level[k];
      uint ridx = k-kb;
      for(m=0;m<3;++m) {
        uint j,jb=row_off[m][i],je=row_off[m][i+1];
        sint *vi; slong *vl; real *vr;
        tuple_list_reserve(&mat[m],mat[m].n+(je-jb));
        vi = mat[m].vi+mat_mi*mat[m].n;
        vl = mat[m].vl+mat_ml*mat[m].n;
        vr = mat[m].vr+mat_mr*mat[m].n;
        mat[m].n+=je-jb;
        for(j=jb;j!=je;++j,vi+=mat_mi,vl+=mat_ml,vr+=mat_mr)
          vi[mat_ridx]=ridx, vi[mat_cidx]=0,
          vl[0]=raw_mat[m][2*j], vr[0]=raw_mat[m][2*j+1];
      }
    }
    Fnloc = F_id.n = ke-kb; Fnglob = gln[lvl];
    tuple_list_reserve(&F_id,F_id.n);
    for(k=kb;k!=ke;++k) {
      F_id.vi[k-kb]=k-kb;
      F_id.vl[k-kb]=cdof->vl[k];
    }
    kb = ke; ke = data->lvl_offset[nl];
    Cnloc = C_id.n = ke-kb;
    Cnglob = 0;
    for(k=lvl+1;k<nl;++k) Cnglob += gln[k];
    tuple_list_reserve(&C_id,C_id.n);
    for(k=kb;k!=ke;++k) {
      C_id.vi[k-kb]=k-kb;
      C_id.vl[k-kb]=cdof->vl[k];
    }
    
    organize_matrix(&data->Q_W[lvl],&data->W[lvl],
                    Fnloc, Cnloc, Cnglob, &C_id, &mat[0],
                    crystal, buf);
    organize_matrix(&data->Q_AfP[lvl],&data->AfP[lvl],
                    Fnloc, Cnloc, Cnglob, &C_id, &mat[1],
                    crystal, buf);
    organize_matrix(&data->Q_Aff[lvl],&data->Aff[lvl],
                    Fnloc, Fnloc, Fnglob, &F_id, &mat[2],
                    crystal, buf);
  }
  
  free(gln);
  free(perm_level);
  for(m=0;m<3;++m) tuple_list_free(&mat[m]);
  tuple_list_free(&C_id);
  tuple_list_free(&F_id);
}

static void read_main_data(amg_data *data, amg_file *f, buffer *buf)
{
  unsigned n;
  if(data->pid==0) {
    double t;
    amg_fread(&t,1,f);
    data->levels = t;
    printf("AMG: %u levels\n", data->levels);
  }
#ifdef MPI
  MPI_Bcast(&data->levels,1,MPI_UNSIGNED,0,data->comm);
#endif
  n = data->levels-1;
  data->cheb_m   = tmalloc(unsigned,n);
  data->cheb_rho = tmalloc(real    ,n);
  if(data->pid==0) {
    double *t; unsigned i;
    buffer_reserve(buf,(2*n+1)*sizeof(double)), t=buf->ptr;
    amg_fread(t,2*n+1,f);
    for(i=0;i<n;++i) data->cheb_m[i]   = *t++;
    for(i=0;i<n;++i) data->cheb_rho[i] = *t++;
    data->tn = *t++;
  }
#ifdef MPI
  MPI_Bcast(data->cheb_m   ,n,MPI_UNSIGNED,0,data->comm);
  MPI_Bcast(data->cheb_rho ,n,REAL_MPI    ,0,data->comm);
  MPI_Bcast(&data->tn      ,1,REAL_MPI    ,0,data->comm);
#endif
  if(data->pid==0) { unsigned i;
    printf("AMG Chebyshev smoother data:\n");
    for(i=0;i<n;++i)
        printf("AMG  level %u: %u iterations with rho = %g\n",
          i+1, data->cheb_m[i], (double)data->cheb_rho[i]);
    printf("AMG: %g rows\n", data->tn);
  }
}

#ifndef AMG_BLOCK_ROWS
#  define AMG_BLOCK_ROWS 1200
#endif

static void read_all_files(amg_data *data, tuple_list *cdof,
                           crystal_data *crystal, buffer* buf)
{
  static double row_data[AMG_BLOCK_ROWS][6];
  ulong ntot, r; unsigned passes;
  uint cdof_ib=0, cdof_ie, cdof_n=cdof->n;
  tuple_list dof_owner, dof_data;
  const int dof_data_mi=5,
            dof_data_proc=0,
            dof_data_level=1,
            dof_data_mat_len=2;
  amg_file file_main, file_mat[3];
  struct { uint *row_len; size_t size, size_old; buffer buf; } raw_mat[3];
  raw_mat[0].row_len = tmalloc(uint,3*(cdof_n+1));
  raw_mat[1].row_len = raw_mat[0].row_len + (cdof_n+1);
  raw_mat[2].row_len = raw_mat[1].row_len + (cdof_n+1);
  raw_mat[0].size=raw_mat[1].size=raw_mat[2].size=0;
  buffer_init(&raw_mat[0].buf,1024);
  buffer_init(&raw_mat[1].buf,1024);
  buffer_init(&raw_mat[2].buf,1024);
  if(data->pid==0) {
    amg_fopen(&file_main,"amg.dat");
    amg_fopen(&file_mat[0],"amg_W.dat");
    amg_fopen(&file_mat[1],"amg_AfP.dat");
    amg_fopen(&file_mat[2],"amg_Aff.dat");
  }  
  read_main_data(data,&file_main,buf);
  ntot = data->tn;
  passes = (ntot + AMG_BLOCK_ROWS - 1)/AMG_BLOCK_ROWS;
  tuple_list_init_max(&dof_owner,1,1,0,AMG_BLOCK_ROWS);
  tuple_list_init_max(&dof_data ,dof_data_mi,1,1,AMG_BLOCK_ROWS);
  for(r=0;r<ntot;r+=AMG_BLOCK_ROWS) {
    unsigned nr = ntot-r>AMG_BLOCK_ROWS ? AMG_BLOCK_ROWS : (unsigned)ntot-r;
    unsigned m, mat_size[3];
    ulong id_end; uint i;
    if(data->pid==0) {
      amg_fread(&row_data[0][0],nr*6,&file_main);
      id_end = row_data[nr-1][0];
      printf("AMG: reading through row %g, pass %u/%u\n",
        (double)id_end, (unsigned)(r/AMG_BLOCK_ROWS+1), passes);
    }
#ifdef MPI
    MPI_Bcast(&id_end,1,ULONG_MPI,0,data->comm);
#endif
    for(i=cdof_ib;i<cdof_n;++i) {
      ulong id;
      if(cdof->vi[cdof_mi*i+cdof_proc]) break;
      id = cdof->vl[i];
      if(id>id_end) break;
      dof_owner.vl[i-cdof_ib] = id;
      dof_owner.vi[i-cdof_ib] = 0;
    }
    cdof_ie = i;
    dof_owner.n = cdof_ie-cdof_ib;
#ifdef MPI
    transfer(1,&dof_owner,0,crystal);
#endif
    if(data->pid==0) {
      sint *outi=dof_data.vi, *ini=dof_owner.vi;
      slong *outl=dof_data.vl, *inl=dof_owner.vl;
      real *outr=dof_data.vr;
      if(dof_owner.n!=nr)
        failwith("AMG: nobody claimed ownership of some rows");
      dof_data.n=nr;
      tuple_list_sort(&dof_owner,1,buf);
      mat_size[0]=mat_size[1]=mat_size[2]=0;
      for(i=0;i<nr;++i) {
        if(*inl++ != row_data[i][0])
          failwith("AMG: row discrepancy");
        *outl++ = row_data[i][0]; /* id    */
        *outi++ = *ini++;         /* proc  */
        *outi++ = row_data[i][1]; /* level */
        *outi++ = row_data[i][2]; /* W   row len */
        *outi++ = row_data[i][3]; /* AfP row len */
        *outi++ = row_data[i][4]; /* Aff row len */
        *outr++ = row_data[i][5]; /* D_ii */
        mat_size[0] += 2*row_data[i][2];
        mat_size[1] += 2*row_data[i][3];
        mat_size[2] += 2*row_data[i][4];
      }
    } else
      dof_data.n=0;
#ifdef MPI
    transfer(1,&dof_data,0,crystal);
#endif
    for(m=0;m<3;++m) raw_mat[m].size_old = raw_mat[m].size;
    for(i=0;i<dof_data.n;++i) {
      if(dof_data.vl[i]!=cdof->vl[cdof_ib+i])
        failwith("AMG: setup problem");
      cdof->vi[cdof_mi*(cdof_ib+i)+cdof_level]
        = dof_data.vi[dof_data_mi*i+dof_data_level]-1;
      cdof->vr[cdof_ib+i]=dof_data.vr[i];
      for(m=0;m<3;++m) {
        unsigned l = dof_data.vi[dof_data_mi*i+dof_data_mat_len+m];
        raw_mat[m].size += 2*l*sizeof(double);
        raw_mat[m].row_len[cdof_ib+i] = l;
      }
    }
    for(m=0;m<3;++m) buffer_reserve(&raw_mat[m].buf,raw_mat[m].size);
    if(data->pid==0) {
      const char *mat_name[3] = {"W","AfP","Aff"};
      for(m=0;m<3;++m) {
        double *ptr, *my_ptr =
          (double*)((char*)raw_mat[m].buf.ptr+raw_mat[m].size_old);
        buffer_reserve(buf,mat_size[m]*sizeof(double)), ptr=buf->ptr;
        printf("AMG:   reading %g MB of %s\n",
               mat_size[m]*sizeof(double)/(1024*1024.0),mat_name[m]);
        amg_fread(ptr,mat_size[m],&file_mat[m]);
        for(i=0;i<nr;++i) {
          int targ = dof_owner.vi[i];
          int len  = row_data[i][2+m];
          int tag  = dof_owner.vl[i];
#ifdef MPI
          if(targ!=0)
            MPI_Send(ptr,2*len,MPI_DOUBLE,targ,tag,data->comm);
          else
#endif
            memcpy(my_ptr,ptr,2*len*sizeof(double)), my_ptr+=2*len;
          ptr += 2*len;
        }
      }
    } else {
#ifdef MPI
      MPI_Status status;
      for(m=0;m<3;++m) {
        double *ptr = (double*)((char*)raw_mat[m].buf.ptr+raw_mat[m].size_old);
        for(i=0;i<dof_data.n;++i) {
          int len = raw_mat[m].row_len[cdof_ib+i];
          int tag = cdof->vl[cdof_ib+i];
          MPI_Recv(ptr,2*len,MPI_DOUBLE,0,tag,data->comm,&status);
          ptr += 2*len;
        }
      }
#endif
    }
    cdof_ib=cdof_ie;
  }
  tuple_list_free(&dof_owner);
  tuple_list_free(&dof_data);
  if(data->pid==0) {
    amg_fclose(&file_main);
    amg_fclose(&file_mat[0]);
    amg_fclose(&file_mat[1]);
    amg_fclose(&file_mat[2]);
  }
  {
    uint *row_len[3] = { raw_mat[0].row_len,
                         raw_mat[1].row_len,
                         raw_mat[2].row_len };
    double *mat_data[3] = { raw_mat[0].buf.ptr,
                            raw_mat[1].buf.ptr,
                            raw_mat[2].buf.ptr };
    organize_data(data,cdof,row_len,mat_data,crystal,buf);
  }
  free(raw_mat[0].row_len);
  buffer_free(&raw_mat[0].buf);
  buffer_free(&raw_mat[1].buf);
  buffer_free(&raw_mat[2].buf);
}

static void separate_cdof(tuple_list *cdof, uint pid, buffer *buf)
{
  sint *proc = cdof->vi + cdof_proc;
  sint *proc_end = proc + cdof_mi * cdof->n;
  for(;proc!=proc_end;proc+=cdof_mi) *proc = (*proc==pid ? 0 : 1);
  tuple_list_sort(cdof,cdof_proc,buf);
}

amg_data *amg_setup(uint n, const ulong *id,
                    uint nz, const uint *Ai, const uint *Aj, const real *A,
                    uint null_space, crystal_data *crystal)
{
  amg_data *data = tmalloc(amg_data,1);
  sint *perm; tuple_list cdof;
  uint i, bufn, chbn;
#ifdef MPI
  buffer *buf = &crystal->all->buf;
#else
  buffer buf_static, *buf=&buf_static;
  buffer_init(buf,1024);
#endif

#ifdef MPI
  MPI_Comm_dup(crystal->comm,&data->comm);
  data->pid = crystal->id, data->np = crystal->num;
#else
  data->pid = 0, data->np = 1;
#endif

  data->null_space = null_space;
  
  perm = tmalloc(sint,n);
  setup_dofs(&cdof,perm, n,id, data->pid, crystal,buf);
  separate_cdof(&cdof,data->pid,buf);
  read_all_files(data,&cdof,crystal,buf);

  {
#ifdef MPI
    jl_comm_t comm = {crystal->id,crystal->num,crystal->comm};
#else
    jl_comm_t comm = {0,1,0};
#endif
    data->gs_top = jl_gs_setup((const slong*)cdof.vl,cdof.n,&comm);
  }
  
  bufn = 0, chbn=0;
  for(i=0;i<data->levels-1;++i) {
    uint fn = data->lvl_offset[i+1]-data->lvl_offset[i];
    if(fn>chbn) chbn=fn;
    if(data->W[i].cn > bufn) bufn=data->W[i].cn;
    if(data->AfP[i].cn > bufn) bufn=data->AfP[i].cn;
    if(data->Aff[i].cn > bufn) bufn=data->Aff[i].cn;
  }
  data->b = tmalloc(real,2*cdof.n+3*chbn+bufn);
  data->x = data->b + cdof.n;
  data->c = data->x + cdof.n;
  data->c_old = data->c + chbn;
  data->r = data->c_old + chbn;
  data->buf = data->r + chbn;
  
  tuple_list_sort(&cdof,cdof_mi,buf);
  for(i=0;i<n;++i) {
    sint p = perm[i];
    if(p!=-1) perm[i] = cdof.vi[cdof_mi*p+cdof_index];
  }
  data->un = n;
  data->cn = cdof.n;
  data->perm = perm;

  data->timing = tmalloc(double,6*(data->levels-1));
  for(i=0;i<6*(data->levels-1);++i) data->timing[i]=0;
  data->timing_n = 0;

  tuple_list_free(&cdof);
#ifndef MPI
  buffer_free(buf);
#endif

  if(data->pid==0) printf("AMG: initialized\n"), fflush(stdout);
#ifdef DIAGNOSTICS  
  check_amg_data(data);
#ifdef MPI
  MPI_Barrier(data->comm);
#endif
  if(data->pid==0) printf("AMG: sanity check complete\n"), fflush(stdout);
#endif
  return data;
}

void amg_free(amg_data *data)
{
  uint i,n;
  free(data->perm);
  free(data->cheb_m);
  free(data->cheb_rho);
  free(data->lvl_offset);
  free(data->Dff);
  free(data->b);
  jl_gs_free(data->gs_top);
  n = 3*(data->levels-1);
  for(i=0;i<n;++i) {
    /* gs_data_free(data->Q_W[i].gs); */
    jl_gs_free(data->Q_W[i].jgs);
    free(data->W[i].row_off);
    free(data->W[i].pr);
  }
  free(data->Q_W);
  free(data->W);
  free(data->timing);
  free(data);
}
