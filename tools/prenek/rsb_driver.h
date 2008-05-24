
#ifndef _rsb_driver_h
#define _rsb_driver_h


#define _2D                  0
#define _3D                  1
#define _MHC                 24
#define _2D_DEFAULT_RES_LIMIT    1.0e-1
#define _3D_DEFAULT_RES_LIMIT    1.0e-1
#define DEFAULT_MAX_REPEATS      1


extern int max_repeats;
extern REAL res_limit;
extern int rsb_dim;
extern REAL *xc, *yc, *zc;
extern int *pre_nek_color;
#endif
