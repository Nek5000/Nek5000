static void condense_matrix(tuple_list *mat, buffer *buf)
{
  unsigned mi=mat->mi,ml=mat->ml,mr=mat->mr;
  slong *vl=mat->vl;
  uint si,di,n=mat->n;
  tuple_list_sort(mat,mi+1,buf);
  tuple_list_sort(mat,mi+0,buf);
  for(di=0,si=0; si<n; ++si) {
    if(vl[ml*si+0]==vl[ml*di+0]&&vl[ml*si+1]==vl[ml*di+1]) {
      if(di!=si) mat->vr[mr*di]+=mat->vr[mr*si];
    } else {
      ++di;
      if(di!=si) {
        memcpy(mat->vi+mi*di,mat->vi+mi*si,mi*sizeof(sint));
        memcpy(mat->vl+ml*di,mat->vl+ml*si,ml*sizeof(slong));
        memcpy(mat->vr+mr*di,mat->vr+mr*si,mr*sizeof(real));
      }
    }
  }
  if(n) mat->n=di+1;
}

#ifndef MPI
#  define MPI_Comm int
#endif
static void write_matrix(tuple_list *mat, buffer *buf, MPI_Comm comm,
                         uint pid, uint np)
{
  FILE *fi,*fj,*fp;
  uint i,n=mat->n;
  double *v;
  buffer_reserve(buf,3*n*sizeof(double));
  v = buf->ptr;
  for(i=0;i<n;++i) *v++ = mat->vl[2*i+0];
  for(i=0;i<n;++i) *v++ = mat->vl[2*i+1];
  for(i=0;i<n;++i) *v++ = mat->vr[i];
  v = buf->ptr;
  if(pid==0) {
    const double magic = 3.14159;
    fi=fopen("amgdmp_i.dat","w");
    fj=fopen("amgdmp_j.dat","w");
    fp=fopen("amgdmp_p.dat","w");
    fwrite(&magic,sizeof(double),1,fi);
    fwrite(&magic,sizeof(double),1,fj);
    fwrite(&magic,sizeof(double),1,fp);
  }
  for(i=0;i<np;++i) {
#ifdef MPI
    MPI_Barrier(comm);
#endif
    if(pid!=0 && i==pid) {
#ifdef MPI
      MPI_Send(&n,1,UINT_MPI,0,i,comm);
      MPI_Send(v,3*n,MPI_DOUBLE,0,np+i,comm);
#endif
    } else if(pid==0) {
      printf("AMG writing data from proc %u\n",i),fflush(stdout);
      if(i!=0) {
#ifdef MPI
        MPI_Status status;
        MPI_Recv(&n,1,UINT_MPI,i,i,comm,&status);
        MPI_Recv(v,3*n,MPI_DOUBLE,i,np+i,comm,&status);
#endif
      }
      fwrite(v    ,sizeof(double),n,fi);
      fwrite(v+  n,sizeof(double),n,fj);
      fwrite(v+2*n,sizeof(double),n,fp);
    }
  }
#ifdef MPI
  MPI_Barrier(comm);
#endif
  if(pid==0) {
    fclose(fi);
    fclose(fj);
    fclose(fp);
  }
}
#ifndef MPI
#  undef MPI_Comm
#endif

typedef void amg_data;

amg_data *amg_setup(uint n, const ulong *id,
                    uint nz, const uint *Ai, const uint *Aj, const real *A,
                    uint null_space, crystal_data *crystal)
{
  uint pid, np;
  sint *perm;
  tuple_list cdof;
  uint i;
  tuple_list mat;
#ifdef MPI
  MPI_Comm comm;
  buffer *buf = &crystal->all->buf;
#else
  int comm=0;
  buffer buf_static, *buf=&buf_static;
#endif
#ifdef MPI
  MPI_Comm_dup(crystal->comm,&comm);
  pid = crystal->id, np = crystal->num;
#else
  pid = 0, np  = 1;
  buffer_init(buf,1024);
#endif

  perm = tmalloc(sint,n);
  setup_dofs(&cdof,perm, n,id, pid, crystal,buf);
  tuple_list_init_max(&mat,1,2,1,nz), mat.n=0;
  for(i=0;i<nz;++i) {
    sint ai = perm[Ai[i]], aj = perm[Aj[i]];
    if(ai==-1 || aj==-1 || fabsr(A[i])==0) continue;
    mat.vl[2*mat.n+0] = cdof.vl[ai];
    mat.vl[2*mat.n+1] = cdof.vl[aj];
    mat.vi[mat.n] = cdof.vi[cdof_mi*ai+cdof_proc];
    mat.vr[mat.n] = A[i];
    ++mat.n;
  }
  condense_matrix(&mat, buf);
#ifdef MPI
  transfer(1,&mat,0,crystal);
  condense_matrix(&mat, buf);
#endif
  write_matrix(&mat,buf,comm,pid,np);
  tuple_list_free(&mat);
  tuple_list_free(&cdof);
  free(perm);
#ifndef MPI
  buffer_free(buf);
#endif
  
  if(pid==0) printf("AMG dump successful\n"), fflush(stdout);
#ifdef MPI
  MPI_Barrier(comm);
#endif
  exit(0);
  
  return 0;
}

void amg_solve(real *x, amg_data *data, const real *b)
{
}

void amg_free(amg_data *data)
{
}

void amg_stats(amg_data *data)
{
}
