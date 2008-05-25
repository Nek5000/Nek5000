
#include "gs_defs.h"

typedef struct {
  uint id, np;
#ifdef MPI
  MPI_Comm comm;
#else
  int comm;
#endif
} jl_comm_t;

typedef struct jl_gs_data_ jl_gs_data;

void jl_gs_op(void *u, gs_dom_t dom, gs_op_t op, unsigned transpose,
              jl_gs_data *gs, void *buf);
jl_gs_data *jl_gs_setup(const slong *id, uint n, const jl_comm_t *comm);
void jl_gs_free(jl_gs_data *gs);

