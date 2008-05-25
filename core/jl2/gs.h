#ifndef GS_H
#define GS_H

#if !defined(ERRMEM_H) || !defined(TYPES_H) || !defined(COMM_H)
#warning "gs.h" requires "errmem.h", "types.h", and "comm.h"
#endif

#include "gs_defs.h"

#ifdef PREFIX
#  define gs_op      TOKEN_PASTE(PREFIX,gs_op     )
#  define gs_op_vec  TOKEN_PASTE(PREFIX,gs_op_vec )
#  define gs_op_many TOKEN_PASTE(PREFIX,gs_op_many)
#  define gs_setup   TOKEN_PASTE(PREFIX,gs_setup  )
#  define gs_free    TOKEN_PASTE(PREFIX,gs_free   )
#  define gs_unique  TOKEN_PASTE(PREFIX,gs_unique )
#endif

typedef struct gs_data_ gs_data;

void gs_op(void *u, gs_dom_t dom, gs_op_t op, unsigned transpose,
           gs_data *gs, void *buf);
void gs_op_vec(void *u, unsigned vn, gs_dom_t dom, gs_op_t op,
               unsigned transpose, gs_data *gs, void *buf);
void gs_op_many(void *const*u, unsigned vn, gs_dom_t dom, gs_op_t op,
                unsigned transpose, gs_data *gs, void *buf);
gs_data *gs_setup(const slong *id, uint n, const comm_t *comm);
void gs_free(gs_data *gs);
void gs_unique(slong *id, uint n, const comm_t *comm);

#endif
