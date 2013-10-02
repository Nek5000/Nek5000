#ifndef FINDPTS_LOCAL_H
#define FINDPTS_LOCAL_H

#define TOKEN_PASTE_(a,b) a##b
#define TOKEN_PASTE(a,b) TOKEN_PASTE_(a,b)

#ifdef PREFIX
#  define PREFIXED_NAME(x) TOKEN_PASTE(PREFIX,x)
#else
#  define PREFIXED_NAME(x) x
#endif

#define uint unsigned long

typedef struct { void *ptr; size_t n,max; } buffer;
struct dbl_range { double min, max; };

typedef void lagrange_fun(double *p,
  double *data, unsigned n, int d, double x);

struct findpts_el_pt_2 {
  double x[2],r[2],oldr[2],dist2,dist2p,tr;
  unsigned index,flags;
};

struct findpts_el_gedge_2 { const double *x[2], *dxdn[2]; };
struct findpts_el_gpt_2   { double x[2], jac[4], hes[4]; };

struct findpts_el_data_2 {
  unsigned npt_max;
  struct findpts_el_pt_2 *p;

  unsigned n[2];
  double *z[2];
  lagrange_fun *lag[2];
  double *lag_data[2];
  double *wtend[2];
  
  const double *x[2];
  
  unsigned side_init;
  double *sides;
  struct findpts_el_gedge_2 edge[4]; /* R S=-1; R S=1; ... */
  struct findpts_el_gpt_2 pt[4];

  double *work;
};

struct findpts_el_pt_3 {
  double x[3],r[3],oldr[3],dist2,dist2p,tr;
  unsigned index,flags;
};

struct findpts_el_gface_3 { const double *x[3], *dxdn[3]; };
struct findpts_el_gedge_3 { const double *x[3], *dxdn1[3], *dxdn2[3],
                                         *d2xdn1[3], *d2xdn2[3]; };
struct findpts_el_gpt_3   { double x[3], jac[9], hes[18]; };

struct findpts_el_data_3 {
  unsigned npt_max;
  struct findpts_el_pt_3 *p;

  unsigned n[3];
  double *z[3];
  lagrange_fun *lag[3];
  double *lag_data[3];
  double *wtend[3];
  
  const double *x[3];
  
  unsigned side_init;
  double *sides;
  struct findpts_el_gface_3 face[6]; /* ST R=-1,R=+1; TR S=-1,S=+1; ... */
  struct findpts_el_gedge_3 edge[12]; /* R S=-1,T=-1; R S=1,T=-1; ... */
  struct findpts_el_gpt_3 pt[8];

  double *work;
};



#define findpts_local_setup_2   PREFIXED_NAME(findpts_local_setup_2)
#define findpts_local_free_2    PREFIXED_NAME(findpts_local_free_2 )
#define findpts_local_2         PREFIXED_NAME(findpts_local_2      )
#define findpts_local_eval_2    PREFIXED_NAME(findpts_local_eval_2 )

struct findpts_local_hash_data_2 {
  uint hash_n;
  struct dbl_range bnd[2];
  double fac[2];
  uint *offset;
  uint max;
};

struct findpts_local_data_2 {
  unsigned ntot;
  const double *elx[2];
  struct obbox_2 *obb;
  struct findpts_local_hash_data_2 hd;
  struct findpts_el_data_2 fed;
  double tol;
};

void findpts_local_setup_2(struct findpts_local_data_2 *const fd,
                           const double *const elx[2],
                           const unsigned n[2], const uint nel,
                           const unsigned m[2], const double bbox_tol,
                           const uint max_hash_size,
                           const unsigned npt_max, const double newt_tol);
void findpts_local_free_2(struct findpts_local_data_2 *const fd);
void findpts_local_2(
        uint   *const  code_base   , const unsigned  code_stride   ,
        uint   *const    el_base   , const unsigned    el_stride   ,
        double *const     r_base   , const unsigned     r_stride   ,
        double *const dist2_base   , const unsigned dist2_stride   ,
  const double *const     x_base[2], const unsigned     x_stride[2],
  const uint npt, struct findpts_local_data_2 *const fd,
  buffer *buf);
void findpts_local_eval_2(
        double *const out_base, const unsigned out_stride,
  const uint   *const  el_base, const unsigned  el_stride,
  const double *const   r_base, const unsigned   r_stride,
  const uint npt,
  const double *const in, struct findpts_local_data_2 *const fd);

#define findpts_local_setup_3   PREFIXED_NAME(findpts_local_setup_3)
#define findpts_local_free_3    PREFIXED_NAME(findpts_local_free_3 )
#define findpts_local_3         PREFIXED_NAME(findpts_local_3      )
#define findpts_local_eval_3    PREFIXED_NAME(findpts_local_eval_3 )

struct findpts_local_hash_data_3 {
  uint hash_n;
  struct dbl_range bnd[3];
  double fac[3];
  uint *offset;
  uint max;
};

struct findpts_local_data_3 {
  unsigned ntot;
  const double *elx[3];
  struct obbox_3 *obb;
  struct findpts_local_hash_data_3 hd;
  struct findpts_el_data_3 fed;
  double tol;
};

void findpts_local_setup_3(struct findpts_local_data_3 *const fd,
                           const double *const elx[3],
                           const unsigned n[3], const uint nel,
                           const unsigned m[3], const double bbox_tol,
                           const uint max_hash_size,
                           const unsigned npt_max, const double newt_tol);
void findpts_local_free_3(struct findpts_local_data_3 *const fd);
void findpts_local_3(
        uint   *const  code_base   , const unsigned  code_stride   ,
        uint   *const    el_base   , const unsigned    el_stride   ,
        double *const     r_base   , const unsigned     r_stride   ,
        double *const dist2_base   , const unsigned dist2_stride   ,
  const double *const     x_base[3], const unsigned     x_stride[3],
  const uint npt, struct findpts_local_data_3 *const fd,
  buffer *buf);
void findpts_local_eval_3(
        double *const out_base, const unsigned out_stride,
  const uint   *const  el_base, const unsigned  el_stride,
  const double *const   r_base, const unsigned   r_stride,
  const uint npt,
  const double *const in, struct findpts_local_data_3 *const fd);

#endif
