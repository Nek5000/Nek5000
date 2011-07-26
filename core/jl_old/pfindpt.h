#ifndef PFINDPT_H
#define PFINDPT_H

/* requires "types.h", "tuple_list.h", and, 
   when MPI is defined, "crystal.h" */
#if !defined(TYPES_H) || !defined(TUPLE_LIST_H) || \
    ( defined(MPI) && !defined(CRYSTAL_H) )
#warning "pfindpt.h" requires "types.h", "tuple_list.h", and "crystal.h"
#endif

typedef struct pfindpt_data_ pfindpt_data;

#ifndef MPI
#  define crystal_data void
#endif

pfindpt_data *pfindpt_setup(unsigned ndim, const real *const*xw,
                            const unsigned *n, uint nel,
                            uint max_hash_size, real bbox_tol,
                            crystal_data *crystal);

#ifndef MPI
#  undef crystal_data
#endif
                            
void pfindpt_free(pfindpt_data *p);
void pfindpt_transfer(pfindpt_data *p, tuple_list *list, int dynamic);
void pfindpt(pfindpt_data *p, tuple_list *list, int guess);
void pfindpt_weights(pfindpt_data *p, const real *r);
real pfindpt_eval(pfindpt_data *p, const real *u);

#endif

