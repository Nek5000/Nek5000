
#if defined(USE_LONG_LONG) || defined(GLOBAL_LONG_LONG)
typedef long long long_long;
#  define WHEN_LONG_LONG(x) x
#else
#  define WHEN_LONG_LONG(x)
#endif

/* the supported domains */
#define GS_FOR_EACH_DOMAIN(macro) \
  macro(double) \
  macro(float ) \
  macro(int   ) \
  macro(long  ) \
  WHEN_LONG_LONG(macro(long_long))
  
/* the supported ops */
#define GS_FOR_EACH_OP(T,macro) \
  macro(T,add) \
  macro(T,mul) \
  macro(T,min) \
  macro(T,max) \
  macro(T,bpr)

/* domain enum */
#define LIST GS_FOR_EACH_DOMAIN(ITEM) gs_dom_n
#define ITEM(T) gs_##T,
typedef enum { LIST } gs_dom_t;
#undef ITEM
#undef LIST

#define gs_sint   TYPE_LOCAL(gs_int,gs_long,gs_long_long)
#define gs_slong TYPE_GLOBAL(gs_int,gs_long,gs_long_long)

/* operation enum */
#define LIST GS_FOR_EACH_OP(T,ITEM) gs_op_n
#define ITEM(T,op) gs_##op,
typedef enum { LIST } gs_op_t;
#undef ITEM
#undef LIST


typedef struct {
  uint id, np;
#ifdef MPI
  MPI_Comm comm;
#else
  int comm;
#endif
} jl_comm_t;

typedef struct jl_gs_data_ jl_gs_data;

typedef enum { gs_pairwise, gs_crystal_router, gs_all_reduce,
               gs_auto } gs_method_t;

void jl_gs(void *u, gs_dom_t dom, gs_op_t op, unsigned transpose,
           jl_gs_data *gsh, void *buf);
jl_gs_data *jl_gs_setup(const slong *id, uint n, const jl_comm_t *comm,
                        int unique, gs_method_t method, int verbose);
void jl_gs_free(jl_gs_data *gsh);

