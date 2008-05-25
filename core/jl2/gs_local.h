#ifndef GS_LOCAL_H
#define GS_LOCAL_H

#if !defined(NAME_H) || !defined(TYPES_H) || !defined(GS_DEFS_H)
#warning "gs_local.h" requires "name.h", "types.h", and "comm.h"
#endif

#ifdef PREFIX
#  define gs_gather_array        TOKEN_PASTE(PREFIX,gs_gather_array       )
#  define gs_init_array          TOKEN_PASTE(PREFIX,gs_init_array         )
#  define gs_gather              TOKEN_PASTE(PREFIX,gs_gather             )
#  define gs_scatter             TOKEN_PASTE(PREFIX,gs_scatter            )
#  define gs_init                TOKEN_PASTE(PREFIX,gs_init               )
#  define gs_gather_vec          TOKEN_PASTE(PREFIX,gs_gather_vec         )
#  define gs_scatter_vec         TOKEN_PASTE(PREFIX,gs_scatter_vec        )
#  define gs_init_vec            TOKEN_PASTE(PREFIX,gs_init_vec           )
#  define gs_gather_many         TOKEN_PASTE(PREFIX,gs_gather_many        )
#  define gs_scatter_many        TOKEN_PASTE(PREFIX,gs_scatter_many       )
#  define gs_init_many           TOKEN_PASTE(PREFIX,gs_init_many          )
#  define gs_gather_vec_to_many  TOKEN_PASTE(PREFIX,gs_gather_vec_to_many )
#  define gs_scatter_many_to_vec TOKEN_PASTE(PREFIX,gs_scatter_many_to_vec)
#  define gs_scatter_vec_to_many TOKEN_PASTE(PREFIX,gs_scatter_vec_to_many)
#endif

void gs_gather_array(void *out, const void *in, uint n,
                     gs_dom_t dom, gs_op_t op);
void gs_init_array(void *out, uint n, gs_dom_t dom, gs_op_t op);

typedef void gs_gather_func_t(
  void *out, const void *in, const unsigned vn,
  const uint *map, gs_dom_t dom, gs_op_t op);
typedef void gs_scatter_func_t(
  void *out, const void *in, const unsigned vn,
  const uint *map, gs_dom_t dom);
typedef void gs_init_func_t(
  void *out, const unsigned vn,
  const uint *map, gs_dom_t dom, gs_op_t op);

extern gs_gather_func_t gs_gather, gs_gather_vec, gs_gather_many,
                        gs_gather_vec_to_many;
extern gs_scatter_func_t gs_scatter, gs_scatter_vec, gs_scatter_many,
                         gs_scatter_many_to_vec, gs_scatter_vec_to_many;
extern gs_init_func_t gs_init, gs_init_vec, gs_init_many;


#endif
