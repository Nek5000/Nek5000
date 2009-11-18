#ifndef GS_H
#define GS_H

#if !defined(COMM_H) || !defined(GS_DEFS_H) || !defined(MEM_H)
#warning "gs.h" requires "comm.h", "gs_defs.h", and "mem.h"
#endif

#define gs               PREFIXED_NAME(gs              )
#define gs_vec           PREFIXED_NAME(gs_vec          )
#define gs_many          PREFIXED_NAME(gs_many         )
#define gs_setup         PREFIXED_NAME(gs_setup        )
#define gs_setup_crystal PREFIXED_NAME(gs_setup_crystal)
#define gs_free          PREFIXED_NAME(gs_free         )
#define gs_unique        PREFIXED_NAME(gs_unique       )

typedef struct gs_data_ gs_data;

void gs(void *u, gs_dom dom, gs_op op, unsigned transpose,
        gs_data *gsh, buffer *buf);
void gs_vec(void *u, unsigned vn, gs_dom dom, gs_op op,
            unsigned transpose, gs_data *gsh, buffer *buf);
void gs_many(void *const*u, unsigned vn, gs_dom dom, gs_op op,
             unsigned transpose, gs_data *gsh, buffer *buf);
gs_data *gs_setup(const slong *id, uint n, const struct comm *comm);
gs_data *gs_setup_crystal(const slong *id, uint n, const struct comm *comm);
void gs_free(gs_data *gsh);
void gs_unique(slong *id, uint n, const struct comm *comm);

#endif
