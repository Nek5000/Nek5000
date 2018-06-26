/*
 * Low-Order finite element preconditioner computed with AMG and Hypre
 *
 * Author: Pedro D. Bello-Maldonado (belloma2@illinois.edu)
 */

// Headers
#include <stdbool.h>
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

// Global variables
extern int n_x, n_y, n_z, n_elem, n_dim;
extern int n_xyz, n_xyze;
extern double *x_m, *y_m, *z_m;
extern long *glo_num;
extern double *pmask;
extern double *binv;
extern int num_procs;
extern int proc_id;
extern int num_loc_dofs;
extern long *dof_map;
extern bool amg_ready;
extern long row_start;
extern long row_end;
extern HYPRE_IJMatrix A_bc;
extern HYPRE_ParCSRMatrix A_fem;
extern HYPRE_IJMatrix B_bc;
extern HYPRE_ParCSRMatrix B_fem;
extern HYPRE_IJVector Bd_bc;
extern HYPRE_ParVector Bd_fem;
extern HYPRE_IJVector f_bc;
extern HYPRE_ParVector f_fem;
extern HYPRE_IJVector u_bc;
extern HYPRE_ParVector u_fem;
extern HYPRE_IJVector Bf_bc;
extern HYPRE_ParVector Bf_fem;
extern HYPRE_IJVector Binv_sem_bc;
extern HYPRE_ParVector Binv_sem;
extern HYPRE_Solver amg_preconditioner;

// Datatypes
typedef double (*Basis)(double*);
typedef void (*DBasis)(double**, double*);

// Fortran functions
long i8glmax_(long*, int*);
void dssum_(double*, int*, int*, int*);

// Interface declaration
void amg_setup_(double*, double*, int*, int*, int*, int*, int*, double*, double*, double*, long*, double*, double*);
void amg_solve_(double*, double*);

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
long maximum(long, long);

// Utilities declaration
void breakpoint_();
void enamge_mpi_output_();

// Memory management
int* mem_alloc_1D_int(int);
long* mem_alloc_1D_long(int);
double* mem_alloc_1D_double(int);
int** mem_alloc_2D_int(int, int);
long** mem_alloc_2D_long(int, int);
double** mem_alloc_2D_double(int, int);
void mem_free_1D_int(int**, int);
void mem_free_1D_long(long**, int);
void mem_free_1D_double(double**, int);
void mem_free_2D_int(int***, int, int);
void mem_free_2D_long(long***, int, int);
void mem_free_2D_double(double***, int, int);
