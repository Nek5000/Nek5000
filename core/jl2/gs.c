#include <stdio.h>

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"

#define gs_op gs_op_t   /* fix conflict with fortran */

#include "gs_defs.h"
#include "gs_local.h"
#include "comm.h"
#include "mem.h"
#include "sort.h"
#include "crystal.h"
#include "sarray_sort.h"
#include "sarray_transfer.h"

#define gs               PREFIXED_NAME(gs              )
#define gs_vec           PREFIXED_NAME(gs_vec          )
#define gs_many          PREFIXED_NAME(gs_many         )
#define gs_setup         PREFIXED_NAME(gs_setup        )
#define gs_setup_crystal PREFIXED_NAME(gs_setup_crystal)
#define gs_free          PREFIXED_NAME(gs_free         )
#define gs_unique        PREFIXED_NAME(gs_unique       )

GS_DEFINE_DOM_SIZES()

typedef enum { mode_plain, mode_vec, mode_many,
               mode_dry_run } gs_mode;

static buffer static_buffer = null_buffer;

static void gather_noop(
  void *out, const void *in, const unsigned vn,
  const uint *map, gs_dom dom, gs_op op)
{}

static void scatter_noop(
  void *out, const void *in, const unsigned vn,
  const uint *map, gs_dom dom)
{}

static void init_noop(
  void *out, const unsigned vn,
  const uint *map, gs_dom dom, gs_op op)
{}

typedef void exec_fun(
  void *data, gs_mode mode, unsigned vn, gs_dom dom, gs_op op,
  unsigned transpose, const void *execdata, const struct comm *comm, char *buf);

/*------------------------------------------------------------------------------
  Topology Discovery
------------------------------------------------------------------------------*/

/* nonzero_ids    (local part)

   Creates an array of s_nonzeros, one per nonzero in user id array. The
   output array is grouped by id. Within each group, non-flagged entries come
   first; otherwise the entries within the group are sorted by the index into
   the user id array. The first index in each group is the primary index, and
   is stored along with each entry. The groups themselves are ordered in
   increasing order of the primary index associated with the group (as opposed
   to the user id). */

typedef struct {
  ulong id; uint i, flag, primary;
} s_nonzero;

static void nonzero_ids(array *nz, const slong *id, const uint n, buffer *buf)
{
  ulong last_id = -(ulong)1;
  uint i, primary = -(uint)1;
  s_nonzero *row, *end;
  array_init(s_nonzero,nz,n), end=row=nz->ptr;
  for(i=0;i<n;++i) {
    slong id_i = id[i], abs_id = iabsl(id_i);
    if(id_i==0) continue;
    end->i = i;
    end->id = abs_id;
    end->flag = id_i!=abs_id;
    ++end;
  }
  nz->n = end-row;
  array_resize(s_nonzero,nz,nz->n);
  sarray_sort_two(s_nonzero,nz->ptr,nz->n, id,1, flag,0, buf);
  for(row=nz->ptr,end=row+nz->n;row!=end;++row) {
    ulong this_id = row->id;
    if(this_id!=last_id) primary = row->i;
    row->primary = primary;
    last_id = this_id;
  }
  sarray_sort(s_nonzero,nz->ptr,nz->n, primary,0, buf);
}

/* shared_ids    (global part)

   Creates an array of s_shareds from an array of s_nonzeros. Each entry in the
   output identifies one id shared with one other processor rp. The primary
   index is i locally and ri remotely. Bit 1 of flags indicates the local flag,
   bit 2 indicates the remote flag. The output has no particular ordering. */

/* construct list of all unique id's on this proc */
typedef struct { ulong id; uint work_proc, src_if; } s_unique;
static void unique_ids(array *un, const array *nz, const uint np)
{
  s_unique *un_row;
  const s_nonzero *nz_row, *nz_end;
  array_init(s_unique,un,nz->n), un_row=un->ptr;
  for(nz_row=nz->ptr,nz_end=nz_row+nz->n;nz_row!=nz_end;++nz_row) {
    if(nz_row->i != nz_row->primary) continue;
    un_row->id = nz_row->id;
    un_row->work_proc = nz_row->id%np;
    un_row->src_if = nz_row->flag ? ~nz_row->i : nz_row->i;
    ++un_row;
  }
  un->n = un_row - (s_unique*)un->ptr;
}

#define FLAGS_LOCAL  1
#define FLAGS_REMOTE 2

typedef struct {
  ulong id; uint i, p, ri, bi; unsigned flags;
} s_shared;

typedef struct {
  ulong id, ord; uint i; unsigned flag;
} s_primary;

typedef struct { ulong id,ord; uint p1, p2, i1f, i2f; } s_shared_work;

static void shared_ids_aux(array *sh, array *pr, uint pr_n,
                           array *wa, buffer *buf)
{
  const s_shared_work *w, *we;
  s_shared *s;
  s_primary *p;
  ulong last_id = -(ulong)1;
  /* translate work array to output arrays */
  sarray_sort(s_shared_work,wa->ptr,wa->n, id,1, buf);
  array_init(s_shared,sh,wa->n), sh->n=wa->n, s=sh->ptr;
  array_init(s_primary,pr,pr_n), p=pr->ptr;
  for(w=wa->ptr,we=w+wa->n;w!=we;++w) {
    uint i1f = w->i1f, i2f = w->i2f;
    uint i1 = ~i1f<i1f?~i1f:i1f, i2 = ~i2f<i2f?~i2f:i2f;
    s->id=w->id, s->i=i1, s->p=w->p2, s->ri=i2;
    s->flags = ((i2f^i2)&FLAGS_REMOTE) | ((i1f^i1)&FLAGS_LOCAL);
    ++s;
    if(w->id!=last_id) {
      last_id=w->id;
      p->id=last_id, p->ord=w->ord, p->i=i1, p->flag=(i1f^i1)&FLAGS_LOCAL;
      ++p;
    }
  }
  pr->n = p-(s_primary*)pr->ptr;
  sarray_sort(s_primary,pr->ptr,pr->n, i,0, buf);
}

static ulong shared_ids(array *sh, array *pr, const array *nz, crystal_data *cr)
{
  array un; s_unique *un_row, *un_end, *other;
  ulong last_id = -(ulong)1;
  ulong ordinal[2], n_shared=0, scan_buf[2];
  array wa; s_shared_work *w;
  uint n_unique;
  /* construct list of all unique id's on this proc */
  unique_ids(&un,nz,cr->comm.np);
  n_unique = un.n;
  /* transfer list to work procs */
  sarray_transfer(s_unique,&un, work_proc, cr);
  /* group by id, put flagged entries after unflagged (within each group) */
  sarray_sort_two(s_unique,un.ptr,un.n, id,1, src_if,0, &cr->data);
  /* count shared id's */
  for(un_row=un.ptr,un_end=un_row+un.n;un_row!=un_end;++un_row) {
    ulong id = un_row->id;
    if(~un_row->src_if<un_row->src_if) continue;
    if(id==last_id) continue;
    other=un_row+1;
    if(other!=un_end&&other->id==id) last_id=id, ++n_shared;
  }
  comm_scan(ordinal, &cr->comm,gs_slong,gs_add, &n_shared,1, scan_buf);
  /* construct list of shared ids */
  last_id = -(ulong)1;
  array_init(s_shared_work,&wa,un.n), wa.n=0, w=wa.ptr;
  for(un_row=un.ptr,un_end=un_row+un.n;un_row!=un_end;++un_row) {
    ulong id = un_row->id;
    uint p1 = un_row->work_proc, i1f = un_row->src_if;
    if(~i1f<i1f) continue;
    for(other=un_row+1;other!=un_end&&other->id==id;++other) {
      uint p2 = other->work_proc, i2f = other->src_if;
      ulong ord;
      if(id!=last_id) last_id=id, ++ordinal[0];
      ord=ordinal[0]-1;
      if(wa.n+2>wa.max)
        array_reserve(s_shared_work,&wa,wa.n+2), w=(s_shared_work*)wa.ptr+wa.n;
      w->id=id, w->ord=ord, w->p1=p1, w->p2=p2, w->i1f=i1f, w->i2f=i2f, ++w;
      w->id=id, w->ord=ord, w->p1=p2, w->p2=p1, w->i1f=i2f, w->i2f=i1f, ++w;
      wa.n+=2;
    }
  }
  /* transfer shared list to source procs */
  sarray_transfer(s_shared_work,&wa, p1, cr);
  /* fill output arrays from work array */
  shared_ids_aux(sh,pr,n_unique,&wa,&cr->data);
  array_free(&un);
  array_free(&wa);
  return ordinal[1];
}

static ulong gs_topology(array *nz, array *sh, array *pr,
                         const slong *id, uint n,
                         crystal_data *cr)
{
  nonzero_ids(nz,id,n,&cr->data);
  return shared_ids(sh,pr, nz,cr);
}

static void make_nz_unique(array *nz, slong *id, buffer *buf)
{
  s_nonzero *p,*e;
  sarray_sort(s_nonzero,nz->ptr,nz->n, i,0, buf);
  for(p=nz->ptr,e=p+nz->n;p!=e;++p)
    if(p->i != p->primary) id[p->i] = -(slong)p->id;
  sarray_sort(s_nonzero,nz->ptr,nz->n, primary,0, buf);
}

static void make_sh_unique(array *sh, slong *id, uint pid, buffer *buf)
{
  s_shared *pb, *pm, *pe, *e;
  /* create sentinel with i = -1 */
  array_reserve(s_shared,sh,sh->n+1);
  ((s_shared*)sh->ptr)[sh->n].i = -(uint)1;
  /* in the sorted list of procs sharing a given id,
     the owner is chosen to be the j^th proc,
     where j = id mod (length of list) */
  sarray_sort_two(s_shared,sh->ptr,sh->n, i,0, p,0, buf);
  for(pb=sh->ptr,e=pb+sh->n;pb!=e;pb=pe) {
    uint i = pb->i;
    /* note: current proc not in list */
    pm = pb; while(pm!=e && pm->i==i && pm->p < pid) ++pm;
    pe = pm; while(pe!=e && pe->i==i) ++pe;
    /* position of current proc in list = pm - pb
       length of list                   = pe - pb + 1 */
    if((uint)(pm-pb) != pb->id % (pe-pb+1)) id[i] = -(slong)pb->id;
  }
}

/*------------------------------------------------------------------------------
  Local setup
------------------------------------------------------------------------------*/

/* assumes nz is sorted by primary, then flag, then index */
static const uint *local_map(const array *nz, const int ignore_flagged)
{
  uint *map, *p, count = 1;
  const s_nonzero *row, *other, *end;
#define DO_COUNT(cond) do \
    for(row=nz->ptr,end=row+nz->n;row!=end;) {                     \
      ulong row_id = row->id; int any=0;                           \
      for(other=row+1;other!=end&&other->id==row_id&&cond;++other) \
        any=2, ++count;                                            \
      count+=any, row=other;                                       \
    } while(0)
  if(ignore_flagged) DO_COUNT(other->flag==0); else DO_COUNT(1);
#undef DO_COUNT
  p = map = tmalloc(uint,count);
#define DO_SET(cond) do \
    for(row=nz->ptr,end=row+nz->n;row!=end;) {                     \
      ulong row_id = row->id; int any=0;                           \
      *p++ = row->i;                                               \
      for(other=row+1;other!=end&&other->id==row_id&&cond;++other) \
        any=1, *p++ = other->i;                                    \
      if(any) *p++ = -(uint)1; else --p;                           \
      row=other;                                                   \
    } while(0)
  if(ignore_flagged) DO_SET(other->flag==0); else DO_SET(1);
#undef DO_SET
  *p = -(uint)1;
  return map;
}

static const uint *flagged_primaries_map(const array *nz)
{
  uint *map, *p, count=1;
  const s_nonzero *row, *end;
  for(row=nz->ptr,end=row+nz->n;row!=end;++row)
    if(row->i==row->primary && row->flag==1) ++count;
  p = map = tmalloc(uint,count);
  for(row=nz->ptr,end=row+nz->n;row!=end;++row)
    if(row->i==row->primary && row->flag==1) *p++ = row->i;
  *p = -(uint)1;
  return map;
}

/*------------------------------------------------------------------------------
  Pairwise Execution
------------------------------------------------------------------------------*/
typedef struct {
  uint n;      /* number of messages */
  uint *p;     /* message source/dest proc */
  uint *size;  /* size of message */
  uint total;  /* sum of message sizes */
} pw_comm_data;

typedef struct {
  pw_comm_data comm[2];
  const uint *map[2];
  comm_req *req;
  uint buffer_size;
} pw_data;

static char *pw_exec_recvs(char *buf, const unsigned unit_size,
                           const struct comm *comm, const pw_comm_data *c,
                           comm_req *req)
{
  const uint *p, *pe, *size=c->size;
  for(p=c->p,pe=p+c->n;p!=pe;++p) {
    size_t len = *(size++)*unit_size;
    comm_irecv(req++,comm,buf,len,*p,*p);
    buf += len;
  }
  return buf;
}

static char *pw_exec_sends(char *buf, const unsigned unit_size,
                           const struct comm *comm, const pw_comm_data *c,
                           comm_req *req)
{
  const uint *p, *pe, *size=c->size;
  for(p=c->p,pe=p+c->n;p!=pe;++p) {
    size_t len = *(size++)*unit_size;
    comm_isend(req++,comm,buf,len,*p,comm->id);
    buf += len;
  }
  return buf;
}

static void pw_exec(
  void *data, gs_mode mode, unsigned vn, gs_dom dom, gs_op op,
  unsigned transpose, const void *execdata, const struct comm *comm, char *buf)
{
  const pw_data *pwd = execdata;
  static gs_scatter_fun *const scatter_to_buf[] =
    { &gs_scatter, &gs_scatter_vec, &gs_scatter_many_to_vec, &scatter_noop };
  static gs_gather_fun *const gather_from_buf[] =
    { &gs_gather, &gs_gather_vec, &gs_gather_vec_to_many, &gather_noop };
  const unsigned recv = 0^transpose, send = 1^transpose;
  unsigned unit_size = vn*gs_dom_size[dom];
  char *sendbuf;
  /* post receives */
  sendbuf = pw_exec_recvs(buf,unit_size,comm,&pwd->comm[recv],pwd->req);
  /* fill send buffer */
  scatter_to_buf[mode](sendbuf,data,vn,pwd->map[send],dom);
  /* post sends */
  pw_exec_sends(sendbuf,unit_size,comm,&pwd->comm[send],
                &pwd->req[pwd->comm[recv].n]);
  comm_wait(pwd->req,pwd->comm[0].n+pwd->comm[1].n);
  /* gather using recv buffer */
  gather_from_buf[mode](data,buf,vn,pwd->map[recv],dom,op);
}

/*------------------------------------------------------------------------------
  Pairwise setup
------------------------------------------------------------------------------*/
static void pw_comm_setup(pw_comm_data *data, array *sh,
                          const unsigned flags_mask, buffer *buf)
{
  uint n=0,count=0, lp=-(uint)1;
  s_shared *s, *se;
  /* sort by remote processor and id (a globally consistent ordering) */
  sarray_sort_two(s_shared,sh->ptr,sh->n, p,0, id,1, buf);
  /* assign index into buffer */
  for(s=sh->ptr,se=s+sh->n;s!=se;++s) {
    if(s->flags&flags_mask) { s->bi = -(uint)1; continue; }
    s->bi = count++;
    if(s->p!=lp) lp=s->p, ++n;
  }
  data->n = n;
  data->p = tmalloc(uint,2*n);
  data->size = data->p + n;
  data->total = count;
  n = 0, lp=-(uint)1;
  for(s=sh->ptr,se=s+sh->n;s!=se;++s) {
    if(s->flags&flags_mask) continue;
    if(s->p!=lp) {
      lp=s->p;
      if(n!=0) data->size[n-1] = count;
      count=0, data->p[n++]=lp;
    }
    ++count;
  }
  if(n!=0) data->size[n-1] = count;
}

static void pw_comm_free(pw_comm_data *data) { free(data->p); }

/* assumes that the bi field of sh is set */
static const uint *pw_map_setup(array *sh, buffer *buf)
{
  uint count=0, *map, *p;
  s_shared *s, *se;
  sarray_sort(s_shared,sh->ptr,sh->n, i,0, buf);
  /* calculate map size */
  count=1;
  for(s=sh->ptr,se=s+sh->n;s!=se;) {
    uint i=s->i;
    if(s->bi==-(uint)1) { ++s; continue; }
    count+=3;
    for(++s;s!=se&&s->i==i;++s) if(s->bi!=-(uint)1) ++count;
  }
  /* write map */
  p = map = tmalloc(uint,count);
  for(s=sh->ptr,se=s+sh->n;s!=se;) {
    uint i=s->i;
    if(s->bi==-(uint)1) { ++s; continue; }
    *p++ = i, *p++ = s->bi;
    for(++s;s!=se&&s->i==i;++s) if(s->bi!=-(uint)1) *p++ = s->bi;
    *p++ = -(uint)1;
  }
  *p = -(uint)1;
  return map;
}

static pw_data *pw_setup(array *sh, buffer *buf)
{
  pw_data *pwd = tmalloc(pw_data,1);
  
  /* default behavior: receive only remotely unflagged data */
  pw_comm_setup(&pwd->comm[0],sh, FLAGS_REMOTE, buf);
  pwd->map[0] = pw_map_setup(sh, buf);

  /* default behavior: send only locally unflagged data */
  pw_comm_setup(&pwd->comm[1],sh, FLAGS_LOCAL, buf);
  pwd->map[1] = pw_map_setup(sh, buf);
  
  pwd->req = tmalloc(comm_req,pwd->comm[0].n+pwd->comm[1].n);
  pwd->buffer_size = pwd->comm[0].total + pwd->comm[1].total;
  return pwd;
}

static void pw_free(pw_data *data)
{
  pw_comm_free(&data->comm[0]);
  pw_comm_free(&data->comm[1]);
  free((uint*)data->map[0]);
  free((uint*)data->map[1]);
  free(data->req);
  free(data);
}

/*------------------------------------------------------------------------------
  Crystal-Router Execution
------------------------------------------------------------------------------*/
typedef struct {
  const uint *scatter_map, *gather_map;
  uint size_r, size_r1, size_r2;
  uint size_sk, size_s, size_total;
  uint p1, p2;
  unsigned nrecvn;
} cr_stage_data;

typedef struct {
  cr_stage_data *stage[2];
  unsigned nstages;
  uint buffer_size, stage_buffer_size;
} cr_data;

static void cr_exec(
  void *data, gs_mode mode, unsigned vn, gs_dom dom, gs_op op,
  unsigned transpose, const void *execdata, const struct comm *comm, char *buf)
{
  const cr_data *crd = execdata;
  static gs_scatter_fun *const scatter_user_to_buf[] =
    { &gs_scatter, &gs_scatter_vec, &gs_scatter_many_to_vec, &scatter_noop };
  static gs_scatter_fun *const scatter_buf_to_buf[] =
    { &gs_scatter, &gs_scatter_vec, &gs_scatter_vec, &gs_scatter };
  static gs_scatter_fun *const scatter_buf_to_user[] =
    { &gs_scatter, &gs_scatter_vec, &gs_scatter_vec_to_many, &scatter_noop };
  static gs_gather_fun *const gather_buf_to_user[] =
    { &gs_gather, &gs_gather_vec, &gs_gather_vec_to_many, &gather_noop };
  static gs_gather_fun *const gather_buf_to_buf[] =
    { &gs_gather, &gs_gather_vec, &gs_gather_vec, &gs_gather };
  const unsigned unit_size = vn*gs_dom_size[dom], nstages=crd->nstages;
  unsigned k;
  char *sendbuf, *buf_old, *buf_new;
  const cr_stage_data *stage = crd->stage[transpose];
  buf_old = buf;
  buf_new = buf_old + unit_size*crd->stage_buffer_size;
  /* crystal router */
  for(k=0;k<nstages;++k) {
    comm_req req[3];
    if(stage[k].nrecvn)
      comm_irecv(&req[1],comm,buf_new,unit_size*stage[k].size_r1,
               stage[k].p1, comm->np+k);
    if(stage[k].nrecvn==2)
      comm_irecv(&req[2],comm,buf_new+unit_size*stage[k].size_r1,
               unit_size*stage[k].size_r2, stage[k].p2, comm->np+k);
    sendbuf = buf_new+unit_size*stage[k].size_r;
    if(k==0)
      scatter_user_to_buf[mode](sendbuf,data,vn,stage[0].scatter_map,dom);
    else
      scatter_buf_to_buf[mode](sendbuf,buf_old,vn,stage[k].scatter_map,dom),
      gather_buf_to_buf [mode](sendbuf,buf_old,vn,stage[k].gather_map ,dom,op);

    comm_isend(&req[0],comm,sendbuf,unit_size*stage[k].size_s,
               stage[k].p1, comm->np+k);
    comm_wait(&req[0],1+stage[k].nrecvn);
    { char *t = buf_old; buf_old=buf_new; buf_new=t; }
  }
  scatter_buf_to_user[mode](data,buf_old,vn,stage[k].scatter_map,dom);
  gather_buf_to_user [mode](data,buf_old,vn,stage[k].gather_map ,dom,op);
}

/*------------------------------------------------------------------------------
  Crystal-Router setup
------------------------------------------------------------------------------*/
static void cr_schedule(cr_data *data, const struct comm *comm)
{
  const uint id = comm->id;
  uint bl=0, n=comm->np;
  unsigned k = 0;
  while(n>1) {
    uint nl = (n+1)/2, bh = bl+nl;
    if(id<bh) n=nl; else n-=nl,bl=bh;
    ++k;
  }
  data->nstages = k;
  data->stage[0] = tmalloc(cr_stage_data,2*(k+1));
  data->stage[1] = data->stage[0] + (k+1);
  bl=0, n=comm->np, k=0;
  while(n>1) {
    uint nl = (n+1)/2, bh = bl+nl;
    uint targ; unsigned recvn;
    recvn = 1, targ = n-1-(id-bl)+bl;
    if(id==targ) targ=bh, recvn=0;
    if(n&1 && id==bh) recvn=2;
    data->stage[1][k].nrecvn=data->stage[0][k].nrecvn=recvn;
    data->stage[1][k].p1    =data->stage[0][k].p1    =targ;
    data->stage[1][k].p2    =data->stage[0][k].p2    =comm->id-1;
    if(id<bh) n=nl; else n-=nl,bl=bh;
    ++k;
  }
}

typedef struct {
  ulong id; uint p, ri, si, bi, send;
} s_crl_work;

/* assumes sh is grouped by i (e.g., sorted by i or by id) */
static void crl_work_init(array *cw, array *sh, const unsigned send_mask,
                          uint this_p)
{
  const unsigned recv_mask = send_mask^(FLAGS_REMOTE|FLAGS_LOCAL);
  uint last_i=-(uint)1; int added_myself;
  uint cw_n = 0, cw_max = cw->max;
  s_crl_work *w = cw->ptr;
  s_shared *s, *se;

#define CW_ADD(aid,ap,ari,asi) do { \
    if(cw_n==cw_max)                                      \
      array_reserve(s_crl_work,cw,cw_n+1),cw_max=cw->max, \
      w=(s_crl_work*)cw->ptr+cw_n;                        \
    w->id=aid, w->p=ap, w->ri=ari, w->si=asi;             \
    ++w, ++cw_n;                                          \
  } while(0)
  
  for(s=sh->ptr,se=s+sh->n;s!=se;++s) {
    int send = (s->flags&send_mask)==0;
    int recv = (s->flags&recv_mask)==0;
    if(s->i!=last_i) last_i=s->i, added_myself=0;
    if(!added_myself && recv && (s->flags&FLAGS_LOCAL)==0) {
      added_myself=1;
      CW_ADD(s->id,this_p,s->i,s->i);
    }
    if(send) CW_ADD(s->id,s->p,s->ri,s->i);
  }
  cw->n=cw_n;
#undef CW_ADD  
}

static void crl_maps(cr_stage_data *stage, array *cw, buffer *buf)
{
  s_crl_work *w, *we, *other;
  uint scount=1, gcount=1, *sp, *gp;
  sarray_sort_two(s_crl_work,cw->ptr,cw->n, bi,0, si,0, buf);
  for(w=cw->ptr,we=w+cw->n;w!=we;w=other) {
    uint bi=w->bi,any=0,si=w->si;
    scount+=3;
    for(other=w+1;other!=we&&other->bi==bi;++other)
      if(other->si!=si) si=other->si, any=2, ++gcount;
    gcount+=any;
  }
  stage->scatter_map = sp = tmalloc(uint,scount+gcount);
  stage->gather_map  = gp = sp + scount;
  for(w=cw->ptr,we=w+cw->n;w!=we;w=other) {
    uint bi=w->bi,any=0,si=w->si;
    *sp++ = w->si, *sp++ = bi;
    *gp++ = bi;
    for(other=w+1;other!=we&&other->bi==bi;++other)
      if(other->si!=si) si=other->si, any=1, *gp++ = si;
    if(any) *gp++ = -(uint)1; else --gp;
    *sp++ = -(uint)1;
  }
  *sp=-(uint)1, *gp=-(uint)1;
}

static uint crl_work_label(array *cw, cr_stage_data *stage,
                           uint cutoff, int send_hi, buffer *buf)
{
  s_crl_work *w, *we, *start;
  uint nsend, nkeep = 0, nks = 0, bi=0;
  /* here w->send has a reverse meaning */
  if(send_hi) for(w=cw->ptr,we=w+cw->n;w!=we;++w) w->send = w->p< cutoff;
         else for(w=cw->ptr,we=w+cw->n;w!=we;++w) w->send = w->p>=cutoff;
  sarray_sort_two(s_crl_work,cw->ptr,cw->n, id,1, send,0, buf);
  for(start=cw->ptr,w=start,we=w+cw->n;w!=we;++w) {
    nkeep += w->send;
    if(w->id!=start->id) start=w;
    if(w->send!=start->send) w->send=0,w->bi=1, ++nks; else w->bi=0;
  }
  nsend = cw->n-nkeep;
  /* assign indices; sent ids have priority (hence w->send is reversed) */
  sarray_sort(s_crl_work,cw->ptr,cw->n, send,0, buf);
  for(start=cw->ptr,w=start,we=w+nsend+nks;w!=we;++w) {
    if(w->id!=start->id) start=w, ++bi;
    if(w->bi!=1) w->send=1;   /* switch back to the usual semantics */
    w->bi = bi;
  }
  stage->size_s = nsend+nks==0 ? 0 : bi+1;
  for(we=(s_crl_work*)cw->ptr+cw->n;w!=we;++w) {
    if(w->id!=start->id) start=w, ++bi;
    w->send = 0;              /* switch back to the usual semantics */
    w->bi = bi;
  }
  stage->size_sk = cw->n==0 ? 0 : bi+1;
  crl_maps(stage,cw,buf);
  return nsend;
}

static void crl_bi_to_si(s_crl_work *w, uint n, uint v) {
  for(;n;--n) w->si=w->bi+v, ++w;
}

static void crl_ri_to_bi(s_crl_work *w, uint n) {
  for(;n;--n) w->bi=w->ri, ++w;
}

static uint cr_learn(array *cw, cr_stage_data *stage, const struct comm *comm,
                     buffer *buf)
{
  comm_req req[3];
  const uint id = comm->id;
  uint bl=0, n=comm->np;
  uint size_max=0;
  uint tag = comm->np;
  while(n>1) {
    uint nl = (n+1)/2, bh = bl+nl;
    uint nkeep, nsend[2], nrecv[2][2] = {{0,0},{0,0}};
    s_crl_work *wrecv[2], *wsend;
    nsend[0] = crl_work_label(cw,stage,bh,id<bh,buf);
    nsend[1] = stage->size_s;
    nkeep = cw->n - nsend[0];

    if(stage->nrecvn   ) comm_irecv(&req[1],comm,nrecv[0],2*sizeof(uint),
                                    stage->p1,tag);
    if(stage->nrecvn==2) comm_irecv(&req[2],comm,nrecv[1],2*sizeof(uint),
                                    stage->p2,tag);
    comm_isend(&req[0],comm,nsend,2*sizeof(uint),stage->p1,tag);
    comm_wait(req,1+stage->nrecvn),++tag;
    
    stage->size_r1 = nrecv[0][1], stage->size_r2 = nrecv[1][1];
    stage->size_r = stage->size_r1 + stage->size_r2;
    stage->size_total = stage->size_r + stage->size_sk;
    if(stage->size_total>size_max) size_max=stage->size_total;
    
    array_reserve(s_crl_work,cw,cw->n+nrecv[0][0]+nrecv[1][0]);
    wrecv[0] = cw->ptr, wrecv[0] += cw->n, wrecv[1] = wrecv[0]+nrecv[0][0];
    wsend = cw->ptr, wsend += nkeep;
    if(stage->nrecvn   )
      comm_irecv(&req[1],comm,wrecv[0],nrecv[0][0]*sizeof(s_crl_work),
                 stage->p1,tag);
    if(stage->nrecvn==2)
      comm_irecv(&req[2],comm,wrecv[1],nrecv[1][0]*sizeof(s_crl_work),
                 stage->p2,tag);
    sarray_sort_two(s_crl_work,cw->ptr,cw->n, send,0, bi,0, buf);
    comm_isend(&req[0],comm,wsend,nsend[0]*sizeof(s_crl_work),stage->p1,tag);
    comm_wait(req,1+stage->nrecvn),++tag;

    crl_bi_to_si(cw->ptr,nkeep,stage->size_r);
    if(stage->nrecvn)    crl_bi_to_si(wrecv[0],nrecv[0][0],0);
    if(stage->nrecvn==2) crl_bi_to_si(wrecv[1],nrecv[1][0],stage->size_r1);
    memmove(wsend,wrecv[0],(nrecv[0][0]+nrecv[1][0])*sizeof(s_crl_work));
    cw->n += nrecv[0][0] + nrecv[1][0];
    cw->n -= nsend[0];
    
    if(id<bh) n=nl; else n-=nl,bl=bh;
    ++stage;
  }
  crl_ri_to_bi(cw->ptr,cw->n);
  crl_maps(stage,cw,buf);
  return size_max;
}

static cr_data *cr_setup(array *sh, const struct comm *comm, buffer *buf)
{
  uint size_max[2];
  array cw = null_array;
  cr_data *crd = tmalloc(cr_data,1);
  
  /* default behavior: receive only remotely unflagged data */
  /* default behavior: send only locally unflagged data */
  
  cr_schedule(crd,comm);

  sarray_sort(s_shared,sh->ptr,sh->n, i,0, buf);
  crl_work_init(&cw,sh, FLAGS_LOCAL , comm->id);
  size_max[0]=cr_learn(&cw,crd->stage[0],comm,buf);
  crl_work_init(&cw,sh, FLAGS_REMOTE, comm->id);
  size_max[1]=cr_learn(&cw,crd->stage[1],comm,buf);
  
  crd->stage_buffer_size = size_max[1]>size_max[0]?size_max[1]:size_max[0];

  array_free(&cw);
  
  crd->buffer_size = 2*crd->stage_buffer_size;
  return crd;
}

static void cr_free_stage_maps(cr_stage_data *stage, unsigned kmax)
{
  unsigned k;
  for(k=0; k<kmax; ++k) {
    free((uint*)stage->scatter_map);
    ++stage;
  }
  free((uint*)stage->scatter_map);
}

static void cr_free(cr_data *data)
{
  cr_free_stage_maps(data->stage[0],data->nstages);
  cr_free_stage_maps(data->stage[1],data->nstages);
  free(data->stage[0]);
  free(data);
}

/*------------------------------------------------------------------------------
  All-reduce Execution
------------------------------------------------------------------------------*/
typedef struct {
  const uint *map_to_buf[2], *map_from_buf[2];
  uint buffer_size;
} allreduce_data;

static void allreduce_exec(
  void *data, gs_mode mode, unsigned vn, gs_dom dom, gs_op op,
  unsigned transpose, const void *execdata, const struct comm *comm, char *buf)
{
  const allreduce_data *ard = execdata;
  static gs_scatter_fun *const scatter_to_buf[] =
    { &gs_scatter, &gs_scatter_vec, &gs_scatter_many_to_vec, &scatter_noop };
  static gs_scatter_fun *const scatter_from_buf[] =
    { &gs_scatter, &gs_scatter_vec, &gs_scatter_vec_to_many, &scatter_noop };
  uint gvn = vn*(ard->buffer_size/2);
  unsigned unit_size = gs_dom_size[dom];
  char *ardbuf;
  ardbuf = buf+unit_size*gvn;
  /* user array -> buffer */
  gs_init_array(buf,gvn,dom,op);
  scatter_to_buf[mode](buf,data,vn,ard->map_to_buf[transpose],dom);
  /* all reduce */
  comm_allreduce(comm,dom,op, buf,gvn, ardbuf);
  /* buffer -> user array */
  scatter_from_buf[mode](data,buf,vn,ard->map_from_buf[transpose],dom);
}

/*------------------------------------------------------------------------------
  All-reduce setup
------------------------------------------------------------------------------*/
static const uint *allreduce_map_setup(array *pr, const unsigned flags_mask,
                                       int to_buf)
{
  s_primary *p, *pe;
  uint count=1, *map, *m;
  for(p=pr->ptr,pe=p+pr->n;p!=pe;++p)
    if((p->flag&flags_mask)==0) count+=3;
  m=map=tmalloc(uint,count);
  if(to_buf) {
    for(p=pr->ptr,pe=p+pr->n;p!=pe;++p)
      if((p->flag&flags_mask)==0)
        *m++ = p->i, *m++ = p->ord, *m++ = -(uint)1;
  } else {
    for(p=pr->ptr,pe=p+pr->n;p!=pe;++p)
      if((p->flag&flags_mask)==0)
        *m++ = p->ord, *m++ = p->i, *m++ = -(uint)1;
  }
  *m=-(uint)1;
  return map;
}

static allreduce_data *allreduce_setup(array *pr, ulong total_shared)
{
  allreduce_data *ard = tmalloc(allreduce_data,1);
  
  /* default behavior: reduce only unflagged data, copy to all */
  ard->map_to_buf  [0] = allreduce_map_setup(pr,1,1);
  ard->map_from_buf[0] = allreduce_map_setup(pr,0,0);

  /* transpose behavior: reduce all data, copy to unflagged */
  ard->map_to_buf  [1] = allreduce_map_setup(pr,0,1);
  ard->map_from_buf[1] = allreduce_map_setup(pr,1,0);
  
  ard->buffer_size = total_shared*2;
  return ard;
}

static void allreduce_free(allreduce_data *ard)
{
  free((uint*)ard->map_to_buf[0]);
  free((uint*)ard->map_to_buf[1]);
  free((uint*)ard->map_from_buf[0]);
  free((uint*)ard->map_from_buf[1]);
  free(ard);
}

/*------------------------------------------------------------------------------
  Main Execution
------------------------------------------------------------------------------*/
typedef struct {
  struct comm comm;
  const uint *map_local[2]; /* 0=unflagged, 1=all */
  const uint *flagged_primaries;
  pw_data *pwd;
  cr_data *crd;
  allreduce_data *ard;
  uint buffer_size;
  void *execdata;
  exec_fun *exec;
} gs_data;

static void gs_aux(
  void *u, gs_mode mode, unsigned vn, gs_dom dom, gs_op op, unsigned transpose,
  gs_data *gsh, buffer *buf)
{
  static gs_scatter_fun *const local_scatter[] =
    { &gs_scatter, &gs_scatter_vec, &gs_scatter_many, &scatter_noop };
  static gs_gather_fun  *const local_gather [] =
    { &gs_gather,  &gs_gather_vec,  &gs_gather_many, &gather_noop  };
  static gs_init_fun *const init[] =
    { &gs_init, &gs_init_vec, &gs_init_many, &init_noop };
  if(!buf) buf = &static_buffer;
  buffer_reserve(buf,vn*gs_dom_size[dom]*gsh->buffer_size);
  local_gather [mode](u,u,vn,gsh->map_local[0^transpose],dom,op);
  if(transpose==0) init[mode](u,vn,gsh->flagged_primaries,dom,op);
  gsh->exec(u,mode,vn,dom,op,transpose,gsh->execdata,&gsh->comm,buf->ptr);
  local_scatter[mode](u,u,vn,gsh->map_local[1^transpose],dom);
}

void gs(void *u, gs_dom dom, gs_op op, unsigned transpose,
        gs_data *gsh, buffer *buf)
{
  gs_aux(u,mode_plain,1,dom,op,transpose,gsh,buf);
}

void gs_vec(void *u, unsigned vn, gs_dom dom, gs_op op,
            unsigned transpose, gs_data *gsh, buffer *buf)
{
  gs_aux(u,mode_vec,vn,dom,op,transpose,gsh,buf);
}

void gs_many(void *const*u, unsigned vn, gs_dom dom, gs_op op,
             unsigned transpose, gs_data *gsh, buffer *buf)
{
  gs_aux((void*)u,mode_many,vn,dom,op,transpose,gsh,buf);
}

/*------------------------------------------------------------------------------
  Main Setup
------------------------------------------------------------------------------*/
static void dry_run_time(double *times, exec_fun *exec, void *data,
                         const struct comm* comm, buffer *buf, uint buf_size)
{
  int i; double t;
  buffer_reserve(buf,gs_dom_size[gs_double]*buf_size);
  for(i= 2;i;--i) exec(0,mode_dry_run,1,gs_double,gs_add,0,data,comm,buf->ptr);
  comm_barrier(comm);
  t = comm_time();
  for(i=10;i;--i) exec(0,mode_dry_run,1,gs_double,gs_add,0,data,comm,buf->ptr);
  t = (comm_time() - t)/10;
  times[0] = t/comm->np, times[1] = t, times[2] = t;
  comm_allreduce(comm,gs_double,gs_add, &times[0],1, &t);
  comm_allreduce(comm,gs_double,gs_min, &times[1],1, &t);
  comm_allreduce(comm,gs_double,gs_max, &times[2],1, &t);
}

static void gs_setup_aux(gs_data *gsh, const slong *id, uint n)
{
  array nz, sh, pr;
  ulong total_shared;
  crystal_data cr;
  crystal_init(&cr,&gsh->comm);
  total_shared = gs_topology(&nz,&sh,&pr, id,n, &cr);
  gsh->map_local[0] = local_map(&nz,1);
  gsh->map_local[1] = local_map(&nz,0);
  gsh->flagged_primaries = flagged_primaries_map(&nz);
  array_free(&nz);

  if(gsh->comm.id==0)
    printf("gs_setup: %ld unique labels shared\n",(long)total_shared);

  gsh->pwd = pw_setup(&sh,&cr.data);
  gsh->buffer_size = gsh->pwd->buffer_size;
  gsh->execdata = gsh->pwd;
  gsh->exec = &pw_exec;
  
  gsh->crd = 0; gsh->ard = 0;
  
  if(gsh->comm.np>1) {
    double time[2][3];

    gsh->crd = cr_setup(&sh,&gsh->comm,&cr.data);
      
    if(gsh->comm.id==0) printf("   pairwise times (avg, min, max): ");
    dry_run_time(time[0],&pw_exec,gsh->pwd,&gsh->comm,
                 &cr.data,gsh->pwd->buffer_size);
    if(gsh->comm.id==0) printf("%g %g %g\n",time[0][0],time[0][1],time[0][2]);

    if(gsh->comm.id==0) printf("   crystal router                : ");
    dry_run_time(time[1],&cr_exec,gsh->crd,&gsh->comm,
                 &cr.data,gsh->crd->buffer_size);
    if(gsh->comm.id==0) printf("%g %g %g\n",time[1][0],time[1][1],time[1][2]);
    if(time[1][2]<time[0][2]) {
      time[0][2] = time[1][2];
      pw_free(gsh->pwd), gsh->pwd=0;
      gsh->buffer_size = gsh->crd->buffer_size;
      gsh->execdata = gsh->crd, gsh->exec = &cr_exec;
    } else
      cr_free(gsh->crd), gsh->crd=0;
    
    if(total_shared<100000) {
      gsh->ard = allreduce_setup(&pr,total_shared);
      if(gsh->comm.id==0) printf("   all reduce                    : ");
      dry_run_time(time[1],&allreduce_exec,gsh->ard,&gsh->comm,
                   &cr.data,gsh->ard->buffer_size);
      if(gsh->comm.id==0)
         printf("%g %g %g\n",time[1][0],time[1][1],time[1][2]);
      if(time[1][2]<time[0][2]) {
        if(gsh->pwd) pw_free(gsh->pwd), gsh->pwd=0;
        if(gsh->crd) cr_free(gsh->crd), gsh->crd=0;
        gsh->buffer_size = gsh->ard->buffer_size;
        gsh->execdata = gsh->ard, gsh->exec = &allreduce_exec;
      } else
        allreduce_free(gsh->ard), gsh->ard=0;
    }

    if(gsh->comm.id==0) {
      if(gsh->pwd) printf("   used all_to_all method: pairwise\n");
      if(gsh->crd) printf("   used all_to_all method: crystal router\n");
      if(gsh->ard) printf("   used all_to_all method: allreduce\n");
    }

  }

  array_free(&pr);
  array_free(&sh);
  crystal_free(&cr);
}

gs_data *gs_setup(const slong *id, uint n, const struct comm *comm)
{
  gs_data *gsh = tmalloc(gs_data,1);
  comm_dup(&gsh->comm,comm);
  gs_setup_aux(gsh,id,n);
  return gsh;
}

gs_data *gs_setup_crystal(const slong *id, uint n, const struct comm *comm)
{
  gs_data *gsh = tmalloc(gs_data,1);
  array nz, sh, pr;
  ulong total_shared;
  crystal_data cr;
  comm_dup(&gsh->comm,comm);
  crystal_init(&cr,&gsh->comm);
  total_shared = gs_topology(&nz,&sh,&pr, id,n, &cr);
  gsh->map_local[0] = local_map(&nz,1);
  gsh->map_local[1] = local_map(&nz,0);
  gsh->flagged_primaries = flagged_primaries_map(&nz);
  array_free(&nz);
  
  gsh->pwd = 0; gsh->ard = 0;
  gsh->crd = cr_setup(&sh,&gsh->comm,&cr.data);
  gsh->buffer_size = gsh->crd->buffer_size;
  gsh->execdata = gsh->crd, gsh->exec = &cr_exec;

  array_free(&pr);
  array_free(&sh);
  crystal_free(&cr);
  return gsh;
}

void gs_free(gs_data *gsh)
{
  comm_free(&gsh->comm);
  free((uint*)gsh->map_local[0]), free((uint*)gsh->map_local[1]);
  free((uint*)gsh->flagged_primaries);
  if(gsh->pwd) pw_free(gsh->pwd);
  if(gsh->crd) cr_free(gsh->crd);
  if(gsh->ard) allreduce_free(gsh->ard);
  free(gsh);
}

void gs_unique(slong *id, uint n, const struct comm *comm)
{
  array nz, sh, pr;
  crystal_data cr;
  crystal_init(&cr,comm);
  gs_topology(&nz,&sh,&pr, id,n, &cr);
  array_free(&pr);
  make_nz_unique(&nz,id,&cr.data);
  make_sh_unique(&sh,id,comm->id,&cr.data);
  array_free(&sh);
  array_free(&nz);
  crystal_free(&cr);
}

/*------------------------------------------------------------------------------
  FORTRAN interface
------------------------------------------------------------------------------*/

#undef gs_op

#undef gs_free
#undef gs_setup
#undef gs_many
#undef gs_vec
#undef gs
#define cgs       PREFIXED_NAME(gs      )
#define cgs_vec   PREFIXED_NAME(gs_vec  )
#define cgs_many  PREFIXED_NAME(gs_many )
#define cgs_setup PREFIXED_NAME(gs_setup)
#define cgs_free  PREFIXED_NAME(gs_free )

#define fgs_setup  FORTRAN_NAME(gs_setup    ,GS_SETUP    )
#define fgs        FORTRAN_NAME(gs_op       ,GS_OP       )
#define fgs_vec    FORTRAN_NAME(gs_op_vec   ,GS_OP_VEC   )
#define fgs_many   FORTRAN_NAME(gs_op_many  ,GS_OP_MANY  )
#define fgs_fields FORTRAN_NAME(gs_op_fields,GS_OP_FIELDS)
#define fgs_free   FORTRAN_NAME(gs_free     ,GS_FREE     )

static gs_data **fgs_info = 0;
static int fgs_max = 0;
static int fgs_n = 0;

void fgs_setup(sint *handle, const slong id[], const sint *n,
               const MPI_Fint *comm, const sint *np)
{
  gs_data *gsh;
  if(fgs_n==fgs_max) fgs_max+=fgs_max/2+1,
                     fgs_info=trealloc(gs_data*,fgs_info,fgs_max);
  gsh=fgs_info[fgs_n]=tmalloc(gs_data,1);
  comm_init_check(&gsh->comm,*comm,*np);
  gs_setup_aux(gsh,id,*n);
  *handle = fgs_n++;
}

static void fgs_check_handle(sint handle, const char *func, unsigned line)
{
  if(handle<0 || handle>=fgs_n || !fgs_info[handle])
    fail(1,__FILE__,line,"%s: invalid handle", func);
}

static void fgs_check_parms(sint handle, sint dom, sint op,
                            const char *func, unsigned line)
{
  if(dom<1 || dom>4)
    fail(1,__FILE__,line,"%s: datatype %d not in valid range 1-4",func,dom);
  if(op <1 || op >4)
    fail(1,__FILE__,line,"%s: op %d not in valid range 1-4",func,op);
  fgs_check_handle(handle,func,line);
}

void fgs(const sint *handle, void *u, const sint *dom, const sint *op,
         const sint *transpose)
{
  fgs_check_parms(*handle,*dom,*op,"gs_op",__LINE__);
  cgs(u,(gs_dom)(*dom-1),(gs_op_t)(*op-1),*transpose!=0,fgs_info[*handle],0);
}

void fgs_vec(const sint *handle, void *u, const sint *n,
             const sint *dom, const sint *op, const sint *transpose)
{
  fgs_check_parms(*handle,*dom,*op,"gs_op_vec",__LINE__);
  cgs_vec(u,*n,(gs_dom)(*dom-1),(gs_op_t)(*op-1),*transpose!=0,
          fgs_info[*handle],0);
}

void fgs_many(const sint *handle, void *u1, void *u2, void *u3,
              void *u4, void *u5, void *u6, const sint *n,
              const sint *dom, const sint *op, const sint *transpose)
{
  void *uu[6];
  uu[0]=u1,uu[1]=u2,uu[2]=u3,uu[3]=u4,uu[4]=u5,uu[5]=u6;
  fgs_check_parms(*handle,*dom,*op,"gs_op_many",__LINE__);
  cgs_many((void *const*)uu,*n,(gs_dom)(*dom-1),(gs_op_t)(*op-1),*transpose!=0,
           fgs_info[*handle],0);
}

static array fgs_fields_array = null_array;

void fgs_fields(const sint *handle,
                void *u, const sint *stride, const sint *n,
                const sint *dom, const sint *op, const sint *transpose)
{
  size_t offset;
  void **p;
  uint i;
  
  fgs_check_parms(*handle,*dom,*op,"gs_op_fields",__LINE__);
  if(*n<0) return;

  array_reserve(void*,&fgs_fields_array,*n);
  p = fgs_fields_array.ptr;
  offset = *stride * gs_dom_size[*dom-1];
  for(i=*n;i;--i) *p++ = u, u = (char*)u + offset;

  cgs_many((void *const*)fgs_fields_array.ptr,*n,
           (gs_dom)(*dom-1),(gs_op_t)(*op-1),
           *transpose!=0, fgs_info[*handle],0);
}

void fgs_free(const sint *handle)
{
  fgs_check_handle(*handle,"gs_free",__LINE__);
  cgs_free(fgs_info[*handle]);
  fgs_info[*handle] = 0;
}

