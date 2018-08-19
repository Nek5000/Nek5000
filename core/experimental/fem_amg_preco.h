#include <stdbool.h>

// Datatypes
typedef double (*Basis)(double*);
typedef void (*DBasis)(double**, double*);

// Interface declaration
void fem_amg_setup(int *precond_type_, int *meshing_type_,
                   int *n_x_, int *n_y_, int *n_z_, int *n_elem_, int *n_dim_,
                   double *x_m_, double *y_m_, double *z_m_,
                   long long *glo_num_, double *pmask_, double* binv_,
                   const MPI_Fint *ce, const int *gsh);
void fem_amg_solve(double *z, double *w);

// FEM assembly declaration
void fem_assembly();
void quadrature_rule(double***, double**, int, int);
void mesh_connectivity(int***, int***, int, int);
void x_map(double**, double*, double**, int, Basis*);
void J_xr_map(double***, double*, double**, int, DBasis*);

// Basis functions and derivatives in 2D
double phi_2D_1(double*);
double phi_2D_2(double*);
double phi_2D_3(double*);
void dphi_2D_1(double**, double*);
void dphi_2D_2(double**, double*);
void dphi_2D_3(double**, double*);

// Basis functions and derivatives in 3D
double phi_3D_1(double*);
double phi_3D_2(double*);
double phi_3D_3(double*);
double phi_3D_4(double*);
void dphi_3D_1(double**, double*);
void dphi_3D_2(double**, double*);
void dphi_3D_3(double**, double*);
void dphi_3D_4(double**, double*);

// Math functions
double determinant(double**, int);
void inverse(double***, double**, int);
HYPRE_Int maximum(HYPRE_Int, HYPRE_Int);

// Memory management
int* mem_alloc_1D_int(int);
HYPRE_Int* mem_alloc_1D_long(int);
double* mem_alloc_1D_double(int);
int** mem_alloc_2D_int(int, int);
double** mem_alloc_2D_double(int, int);
void mem_free_1D_int(int**, int);
void mem_free_1D_long(HYPRE_Int**, int);
void mem_free_1D_double(double**, int);
void mem_free_2D_int(int***, int, int);
void mem_free_2D_double(double***, int, int);
