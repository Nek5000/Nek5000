#ifndef GS_H
#define GS_H

/* requires "types.h", and, when MPI is defined, "crystal.h" */
#if !defined(TYPES_H) || ( defined(MPI) && !defined(CRYSTAL_H) )
#warning "gs.h" requires "types.h" and "crystal.h"
#endif

typedef struct gs_data_ gs_data;

#ifndef MPI
#  define crystal_data void
#endif

gs_data *gs_data_setup(uint n, const ulong *label,
                       uint maxv, crystal_data *crystal);

void gs_data_stats(double stats[3], const gs_data *data);

#ifndef MPI
#  undef crystal_data
#endif

void gs_data_free(gs_data *data);
double gs_op(real *u, int op, const gs_data *data);
double gs_op_vec(real *u, uint n, int op, const gs_data *data);
double gs_op_many(real **u, uint n, int op, const gs_data *data);

#define GS_OP_ADD 1
#define GS_OP_MUL 2
#define GS_OP_MIN 3
#define GS_OP_MAX 4
#define GS_OP_BPR 5

#if defined(TUPLE_LIST_H)
  void gs_data_dump(tuple_list *dump, const gs_data *data);
#endif

#endif

