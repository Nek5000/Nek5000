#ifndef FINDPTS_EL_H
#define FINDPTS_EL_H

#if !defined(TYPES_H) || !defined(NAME_H)
#warning "findpts_el.h" requires "types.h" and "name.h"
#endif

#define findpts_el_3_setup   PREFIXED_NAME(findpts_el_3_setup )
#define findpts_el_3_free    PREFIXED_NAME(findpts_el_3_free  )
#define findpts_el_3         PREFIXED_NAME(findpts_el_3       )

struct findpts_el_3_pt { double x[3], r[3]; double dist2; unsigned flags; };

struct findpts_el_3_gface { const double *x[3], *dxdn[3]; };
struct findpts_el_3_gedge { const double *x[3], *dxdn1[3], *dxdn2[3],
                                         *d2xdn1[3], *d2xdn2[3]; };
struct findpts_el_3_gpt   { double x[3], jac[9], hes[18]; };

struct findpts_el_3_data {
  struct findpts_el_3_pt *p;

  unsigned n[3];
  lagrange_fun *lag[3];
  double *lag_data[3];
  double *wtend[3];
  
  const double *x[3];
  
  unsigned side_init;
  double *sides;
  struct findpts_el_3_gface face[6]; /* ST R=-1,R=+1; TR S=-1,S=+1; ... */
  struct findpts_el_3_gedge edge[12]; /* R S=-1,T=-1; R S=1,T=-1; ... */
  struct findpts_el_3_gpt pt[8];

  double *work;
};

struct findpts_el_3_data *findpts_el_3_setup(
  const unsigned nr, const unsigned ns, const unsigned nt,
  const unsigned npt_max);
void findpts_el_3_free(struct findpts_el_3_data *fd);
void findpts_el_3(struct findpts_el_3_data *fd, uint npt, const double tol);

static void findpts_el_3_start(struct findpts_el_3_data *fd)
{
  fd->side_init=0;
}

static struct findpts_el_3_pt *findpts_el_3_points(struct findpts_el_3_data *fd)
{
  return fd->p;
}

#endif
