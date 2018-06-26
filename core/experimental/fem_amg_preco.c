#include "name.h"

#define fem_amg_setup FORTRAN_UNPREFIXED(fem_amg_setup, FEM_AMG_SETUP)
#define fem_amg_solve FORTRAN_UNPREFIXED(fem_amg_solve, FEM_AMG_SOLVE)

#ifdef HYPRE
/*
 * Low-Order finite element preconditioner computed with AMG and Hypre
 *
 * Author: Pedro D. Bello-Maldonado (belloma2@illinois.edu)
 */

// Headers
#include "fem_amg_preco.h"
#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_ls.h"

// Global variables
int precond_type;
int meshing_type;
int n_x, n_y, n_z, n_elem, n_dim;
int n_xyz, n_xyze;
double *x_m, *y_m, *z_m;
long *glo_num;
double *pmask;
double *binv;
int num_loc_dofs;
long *dof_map;
bool amg_ready = false;
long row_start;
long row_end;
HYPRE_IJMatrix A_bc;
HYPRE_ParCSRMatrix A_fem;
HYPRE_IJMatrix B_bc;
HYPRE_ParCSRMatrix B_fem;
HYPRE_IJVector Bd_bc;
HYPRE_ParVector Bd_fem;
HYPRE_IJVector f_bc;
HYPRE_ParVector f_fem;
HYPRE_IJVector u_bc;
HYPRE_ParVector u_fem;
HYPRE_IJVector Bf_bc;
HYPRE_ParVector Bf_fem;
HYPRE_IJVector Binv_sem_bc;
HYPRE_ParVector Binv_sem;
HYPRE_Solver amg_preconditioner;

MPI_Comm comm;
int num_procs;
int proc_id;

// Interface definition
void fem_amg_setup(int *precond_type_, int *meshing_type_, 
                   int *n_x_, int *n_y_, int *n_z_, int *n_elem_, int *n_dim_, 
                   double *x_m_, double *y_m_, double *z_m_, 
                   long *glo_num_, double *pmask_, double* binv_,
                   const MPI_Fint *ce)
{
    // Parallel run information
    MPI_Comm_dup(MPI_Comm_f2c(*ce),&comm);
    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &proc_id);

    double time0 = MPI_Wtime();
    if (proc_id == 0) printf("fem_amg_setup:");
    setbuf(stdout, NULL);

    // Mesh structure
    precond_type = *precond_type_;
    meshing_type = *meshing_type_;
    n_x = *n_x_;
    n_y = *n_y_;
    n_z = *n_z_;
    n_elem = *n_elem_;
    n_dim = *n_dim_;
    n_xyz = n_x * n_y * n_z;
    n_xyze = n_x * n_y * n_z * n_elem;
    x_m = x_m_;
    y_m = y_m_;
    z_m = z_m_;
    glo_num = glo_num_;
    pmask = pmask_;
    binv = binv_;

    // Assemble FEM matrices and vectors
    fem_assembly();

    // Initialize preconditioning vectors
    int row;
    int row_start = hypre_ParCSRMatrixFirstRowIndex(A_fem);
    int row_end = hypre_ParCSRMatrixLastRowIndex(A_fem);

    HYPRE_IJVectorCreate(comm, row_start, row_end, &u_bc);
    HYPRE_IJVectorSetObjectType(u_bc, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(u_bc);
    HYPRE_IJVectorAssemble(u_bc);
    HYPRE_IJVectorGetObject(u_bc, (void**) &u_fem);

    HYPRE_IJVectorCreate(comm, row_start, row_end, &f_bc);
    HYPRE_IJVectorSetObjectType(f_bc, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(f_bc);
    HYPRE_IJVectorAssemble(f_bc);
    HYPRE_IJVectorGetObject(f_bc, (void**) &f_fem);

    HYPRE_IJVectorCreate(comm, row_start, row_end, &Bf_bc);
    HYPRE_IJVectorSetObjectType(Bf_bc, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(Bf_bc);
    HYPRE_IJVectorAssemble(Bf_bc);
    HYPRE_IJVectorGetObject(Bf_bc, (void**) &Bf_fem);

    HYPRE_IJVectorCreate(comm, row_start, row_end, &Binv_sem_bc);
    HYPRE_IJVectorSetObjectType(Binv_sem_bc, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(Binv_sem_bc);

    for (row = row_start; row <= row_end; row++)
        HYPRE_IJVectorSetValues(Binv_sem_bc, 1, &row, &(binv[dof_map[row - row_start]]));

    HYPRE_IJVectorAssemble(Binv_sem_bc);
    HYPRE_IJVectorGetObject(Binv_sem_bc, (void**) &Binv_sem);

    // AMG parameters
    HYPRE_BoomerAMGCreate(&amg_preconditioner);

    // Set some parameters (See Reference Manual for more parameters)
    HYPRE_BoomerAMGSetMaxRowSum(amg_preconditioner, 1); // Don't check for maximum row sum
    HYPRE_BoomerAMGSetCoarsenType(amg_preconditioner, 10); // 0 for CLJP, 6 for Falgout, 10, HMIS
    HYPRE_BoomerAMGSetInterpType(amg_preconditioner, 6); // Interpolation type, 0 for classical modified interpolation, 6 for extended+i interpolation
    HYPRE_BoomerAMGSetPMaxElmts(amg_preconditioner, 4); // Maximum number of elements per row for interpolation
    HYPRE_BoomerAMGSetAggNumLevels(amg_preconditioner, 0); // 0 for no-aggressive coarsening
    HYPRE_BoomerAMGSetStrongThreshold(amg_preconditioner, 0.25); // Strength threshold
    HYPRE_BoomerAMGSetMaxCoarseSize(amg_preconditioner, 50); // maximum number of rows in coarse level
    HYPRE_BoomerAMGSetRelaxType(amg_preconditioner, 3); // G-S/Jacobi hybrid relaxation, 3 means SOR
    HYPRE_BoomerAMGSetPrintLevel(amg_preconditioner, 0);  // print solve info + parameters
    HYPRE_BoomerAMGSetMaxIter(amg_preconditioner, 1); // maximum number of V-cycles
    HYPRE_BoomerAMGSetTol(amg_preconditioner, 0); // convergence tolerance
    HYPRE_BoomerAMGSetNonGalerkinTol(amg_preconditioner, 0.1); // Non-Galerkin tolerance
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(amg_preconditioner, 0.0, 0); // Skip first level when droping entries

    // Setup preconditioner
    HYPRE_BoomerAMGSetup(amg_preconditioner, A_fem, NULL, NULL);
    amg_ready = true;

    double time1 = MPI_Wtime();
    if (proc_id == 0) printf("done %fs\n",time1-time0);
    setbuf(stdout, NULL);
}

void fem_amg_solve(double *z, double *w)
{
    /*
     * Solves the system $\boldsymbol{M} \boldsymbol{z} = \boldsymbol{w}$ using Algebraic Multigrid
     */

    // Variables
    int row, idx;

    // Check if the AMG matrices have been setup correctly
    if (!amg_ready)
    {
        MPI_Barrier(comm);

        if (proc_id == 0)
            printf("ERROR: AMG hasn't been setup. Call the amg_setup function first before calling this function\n");

        exit(EXIT_FAILURE);
    }

    // Choose preconditioner type
    HYPRE_IJVectorInitialize(f_bc);

    switch (precond_type)
    {
        // No mass operator
        case 1:
            for (row = row_start; row <= row_end; row++)
                HYPRE_IJVectorSetValues(f_bc, 1, &row, &(w[dof_map[row - row_start]]));

            break;

        // Diagonal mass operator
        case 2:
            for (row = row_start; row <= row_end; row++)
            {
                double Bd_fem_value;
                double Binv_sem_value;

                HYPRE_IJVectorGetValues(Bd_bc, 1, &row, &Bd_fem_value);
                HYPRE_IJVectorGetValues(Binv_sem_bc, 1, &row, &Binv_sem_value);

                double value = Bd_fem_value * Binv_sem_value * w[dof_map[row - row_start]];

                HYPRE_IJVectorSetValues(f_bc, 1, &row, &value);
            }

            break;

        // Full mass operator
        case 3:
            for (row = row_start; row <= row_end; row++)
            {
                double Binv_sem_value;

                HYPRE_IJVectorGetValues(Binv_sem_bc, 1, &row, &Binv_sem_value);

                double value = Binv_sem_value * w[dof_map[row - row_start]];

                HYPRE_IJVectorSetValues(Bf_bc, 1, &row, &value);
            }

            HYPRE_ParCSRMatrixMatvec(1.0, B_fem, Bf_fem, 0.0, f_fem);

            break;
    }

    HYPRE_IJVectorAssemble(f_bc);

    // Solve preconditioned system
    HYPRE_BoomerAMGSolve(amg_preconditioner, A_fem, f_fem, u_fem);

    // Map data back to SEM mesh
    for (idx = 0; idx < n_xyze; idx++)
        z[idx] = 0.0;

    for (row = row_start; row <= row_end; row++)
        HYPRE_IJVectorGetValues(u_bc, 1, &row, &(z[dof_map[row - row_start]]));

    dssum_(z, &n_x, &n_y, &n_z);

}

// FEM Assembly definition
void fem_assembly()
{
    /*
     * Assembles the low-order FEM matrices from the spectral element mesh
     *
     * Returns A_fem and B_fem
     */

    // Variables
    int i, j, k, e, d, t, q;
    int row, idx;

    // Rank and prepare data to be mapped to rows of the matrix so it can 
    // be solved with Hypre (Ranking done on the Fortran side since here it fails
    // for some reason I couldn't figure out)
    long *ranking = mem_alloc_1D_long(n_xyze);

    // Ranking fix: If there are no Dirichlet pressure nodes offset is 1, otherwise is 2
    int one = 1;
    long offset = 1;

    for (idx = 0; idx < n_xyze; idx++)
    {
        if (pmask[idx] == 0.0)
        {
            offset = 2;
            break;
        }
    }

    offset = i8glmax_(&offset, &one);

    for (idx = 0; idx < n_xyze; idx++)
        ranking[idx] = glo_num[idx] - offset;

    row_start = 0;
    row_end = 0;

    for (idx = 0; idx < n_xyze; idx++)
        if (ranking[idx] >= 0) row_end = maximum(row_end, ranking[idx]);

    if (proc_id < num_procs - 1)
        MPI_Send(&row_end, 1, MPI_LONG, proc_id + 1, 0, comm);

    if (proc_id > 0)
    {
        MPI_Recv(&row_start, 1, MPI_LONG, proc_id - 1, 0, comm, MPI_STATUS_IGNORE);
        row_start += 1;
    }

    num_loc_dofs = row_end - row_start + 1;
    dof_map = mem_alloc_1D_long(num_loc_dofs);

    for (idx = 0; idx < n_xyze; idx++)
    {
        if ((row_start <= ranking[idx]) && (ranking[idx] <= row_end))
        {
            dof_map[ranking[idx] - row_start] = idx;
        }
    }

    // Assemble FE matrices with boundary conditions applied
    HYPRE_IJMatrixCreate(comm, row_start, row_end, row_start, row_end, &A_bc);
    HYPRE_IJMatrixSetObjectType(A_bc, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A_bc);

    HYPRE_IJMatrixCreate(comm, row_start, row_end, row_start, row_end, &B_bc);
    HYPRE_IJMatrixSetObjectType(B_bc, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(B_bc);

    HYPRE_IJVectorCreate(comm, row_start, row_end, &Bd_bc);
    HYPRE_IJVectorSetObjectType(Bd_bc, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(Bd_bc);

    double *Bd_sum = mem_alloc_1D_double(n_xyze);
    for (i = 0; i < n_xyze; i++) Bd_sum[i] = 0.0;

    // Set quadrature rule
    int n_quad = (n_dim == 2) ? 3 : 4;
    double **q_r;
    double *q_w;

    quadrature_rule(&q_r, &q_w, n_quad, n_dim);

    // Mesh connectivity (Can be changed to fill-out or one-per-vertex)
    int num_fem;
    int **v_coord;
    int **t_map;

    if (n_dim == 2)
        num_fem = (meshing_type == 0) ? 4 : 2;
    else
        num_fem = (meshing_type == 0) ? 8 : 6;

    mesh_connectivity(&v_coord, &t_map, num_fem, n_dim);

    // Finite element assembly
    double **A_loc = mem_alloc_2D_double(n_dim + 1, n_dim + 1);
    double **B_loc = mem_alloc_2D_double(n_dim + 1, n_dim + 1);
    double **J_xr = mem_alloc_2D_double(n_dim, n_dim);
    double **J_rx = mem_alloc_2D_double(n_dim, n_dim);
    double **x_t = mem_alloc_2D_double(n_dim, n_dim + 1);
    double *q_x = mem_alloc_1D_double(n_dim);
    double *dp = mem_alloc_1D_double(n_dim);

    Basis phi[n_dim + 1];
    DBasis dphi[n_dim + 1];

    if (n_dim == 2)
    {
        phi[0] = phi_2D_1; phi[1] = phi_2D_2; phi[2] = phi_2D_3;
        dphi[0] = dphi_2D_1; dphi[1] = dphi_2D_2; dphi[2] = dphi_2D_3;
    }
    else
    {
        phi[0] = phi_3D_1; phi[1] = phi_3D_2; phi[2] = phi_3D_3; phi[3] = phi_3D_4;
        dphi[0] = dphi_3D_1; dphi[1] = dphi_3D_2; dphi[2] = dphi_3D_3; dphi[3] = dphi_3D_4;
    }

    int s_x, s_y, s_z;
    int E_x = n_x - 1;
    int E_y = n_y - 1;
    int E_z = (n_dim == 2) ? 1 : n_z - 1;

    for (e = 0; e < n_elem; e++)
    {
        // Cycle through collocated quads/hexes
        for (s_z = 0; s_z < E_z; s_z++)
        {
            for (s_y = 0; s_y < E_y; s_y++)
            {
                for (s_x = 0; s_x < E_x; s_x++)
                {
                    // Get indices
                    int s[n_dim];

                    if (n_dim == 2)
                    {
                        s[0] = s_x; s[1] = s_y;
                    }
                    else
                    {
                        s[0] = s_x; s[1] = s_y; s[2] = s_z;
                    }

                    int idx[(int)(pow(2, n_dim))];

                    for (i = 0; i < pow(2, n_dim); i++)
                    {
                        idx[i] = 0;

                        for (d = 0; d < n_dim; d++)
                        {
                            idx[i] += (s[d] + v_coord[i][d]) * pow(n_x, d);
                        }
                    }

                    // Cycle through collocated triangles/tets
                    for (t = 0; t < num_fem; t++)
                    {
                        // Get vertices
                        for (i = 0; i < n_dim + 1; i++)
                        {
                            if (n_dim == 2)
                            {
                                x_t[0][i] = x_m[idx[t_map[t][i]] + e * n_xyz];
                                x_t[1][i] = y_m[idx[t_map[t][i]] + e * n_xyz];
                            }
                            else
                            {
                                x_t[0][i] = x_m[idx[t_map[t][i]] + e * n_xyz];
                                x_t[1][i] = y_m[idx[t_map[t][i]] + e * n_xyz];
                                x_t[2][i] = z_m[idx[t_map[t][i]] + e * n_xyz];
                            }
                        }

                        // Local FEM matrices
                        // Reset local stiffness and mass matrices
                        for (i = 0; i < n_dim + 1; i++)
                        {
                            for (j = 0; j < n_dim + 1; j++)
                            {
                                A_loc[i][j] = 0.0;
                                B_loc[i][j] = 0.0;
                            }
                        }

                        // Build local stiffness matrices by applying quadrature rules
                        for (q = 0; q < n_quad; q++)
                        {
                            // From r to x
                            x_map(&q_x, q_r[q], x_t, n_dim, phi);
                            J_xr_map(&J_xr, q_r[q], x_t, n_dim, dphi);
                            inverse(&J_rx, J_xr, n_dim);
                            double det_J_xr = determinant(J_xr, n_dim);

                            // Integrand
                            for (i = 0; i < n_dim + 1; i++)
                            {
                                for (j = 0; j < n_dim + 1; j++)
                                {
                                    int alpha, beta;
                                    double func = 0.0;

                                    for (alpha = 0; alpha < n_dim; alpha++)
                                    {
                                        double a = 0.0, b = 0.0;

                                        for (beta = 0; beta < n_dim; beta++)
                                        {
                                            dphi[i](&dp, q_r[q]);
                                            a += dp[beta] * J_rx[beta][alpha];

                                            dphi[j](&dp, q_r[q]);
                                            b += dp[beta] * J_rx[beta][alpha];
                                        }

                                        func += a * b;
                                    }

                                    A_loc[i][j] += func * det_J_xr * q_w[q];
                                    B_loc[i][j] += phi[i](q_r[q]) * phi[j](q_r[q]) * det_J_xr * q_w[q];
                                }
                            }
                        }

                        // Add to global matrix
                        for (i = 0; i < n_dim + 1; i++)
                        {
                            for (j = 0; j < n_dim + 1; j++)
                            {
                                if ((pmask[idx[t_map[t][i]] + e * n_xyz] > 0.0) && (pmask[idx[t_map[t][j]] + e * n_xyz] > 0.0))
                                {
                                    int row = ranking[idx[t_map[t][i]] + e * n_xyz];
                                    int col = ranking[idx[t_map[t][j]] + e * n_xyz];

                                    double A_val = A_loc[i][j];
                                    double B_val = B_loc[i][j];

                                    int ncols = 1;
                                    int insert_error;

                                    if (fabs(A_val) > 1.0e-14)
                                    {
                                        insert_error = HYPRE_IJMatrixAddToValues(A_bc, 1, &ncols, &row, &col, &A_val);
                                    }

                                    if (fabs(B_val) > 1.0e-14)
                                    {
                                        insert_error = HYPRE_IJMatrixAddToValues(B_bc, 1, &ncols, &row, &col, &B_val);
                                    }

                                    if (insert_error != 0)
                                    {
                                        if (proc_id == 0)
                                            printf("There was an error with entry A(%d, %d) = %f or B(%d, %d) = %f\n", row, col, A_val, row, col, B_val);

                                        exit(EXIT_FAILURE);
                                    }
                                }

                                int row = idx[t_map[t][i]] + e * n_xyz;
                                int col = idx[t_map[t][j]] + e * n_xyz;
                                Bd_sum[row] += B_loc[i][j];
                                Bd_sum[col] += B_loc[i][j];
                            }
                        }
                    }
                }
            }
        }
    }

    for (idx = 0; idx < n_xyze; idx++)
        Bd_sum[idx] /= 2.0;

    dssum_(Bd_sum, &n_x, &n_y, &n_z);

    for (row = row_start; row <= row_end; row++)
        HYPRE_IJVectorSetValues(Bd_bc, 1, &row, &(Bd_sum[dof_map[row - row_start]]));

    HYPRE_IJMatrixAssemble(A_bc);
    HYPRE_IJMatrixGetObject(A_bc, (void**) &A_fem);

    HYPRE_IJMatrixAssemble(B_bc);
    HYPRE_IJMatrixGetObject(B_bc, (void**) &B_fem);

    HYPRE_IJVectorAssemble(Bd_bc);
    HYPRE_IJVectorGetObject(Bd_bc, (void**) &Bd_fem);

    // Free memory
    mem_free_1D_double(&Bd_sum, n_xyze);
    mem_free_2D_double(&q_r, n_quad, n_dim);
    mem_free_1D_double(&q_w, n_quad);
    mem_free_2D_int(&v_coord, pow(n_dim, 2), n_dim);
    mem_free_2D_int(&t_map, num_fem, n_dim + 1);
    mem_free_2D_double(&A_loc, n_dim + 1, n_dim + 1);
    mem_free_2D_double(&B_loc, n_dim + 1, n_dim + 1);
    mem_free_2D_double(&J_xr, n_dim, n_dim);
    mem_free_2D_double(&J_rx, n_dim, n_dim);
    mem_free_2D_double(&x_t, n_dim, n_dim + 1);
    mem_free_1D_double(&q_x, n_dim);
    mem_free_1D_double(&dp, n_dim);
}

void quadrature_rule(double ***q_r, double **q_w, int n_quad, int n_dim)
{
    (*q_r) = mem_alloc_2D_double(n_quad, n_dim);
    (*q_w) = mem_alloc_1D_double(n_quad);

    if (n_dim == 2)
    {
        if (n_quad == 3)
        {
            (*q_r)[0][0] = 1.0 / 6.0; (*q_r)[0][1] = 1.0 / 6.0;
            (*q_r)[1][0] = 2.0 / 3.0; (*q_r)[1][1] = 1.0 / 6.0;
            (*q_r)[2][0] = 1.0 / 6.0; (*q_r)[2][1] = 2.0 / 3.0;

            (*q_w)[0] = 1.0 / 6.0;
            (*q_w)[1] = 1.0 / 6.0;
            (*q_w)[2] = 1.0 / 6.0;
        }
        else if (n_quad == 4)
        {
            (*q_r)[0][0] = 1.0 / 3.0; (*q_r)[0][1] = 1.0 / 3.0;
            (*q_r)[1][0] = 1.0 / 5.0; (*q_r)[1][1] = 3.0 / 5.0;
            (*q_r)[2][0] = 1.0 / 5.0; (*q_r)[2][1] = 1.0 / 5.0;
            (*q_r)[3][0] = 3.0 / 5.0; (*q_r)[3][1] = 1.0 / 5.0;

            (*q_w)[0] = - 27.0 / 96.0;
            (*q_w)[1] = 25.0 / 96.0;
            (*q_w)[2] = 25.0 / 96.0;
            (*q_w)[3] = 25.0 / 96.0;
        }
        else
        {
            printf("No quadrature rule for %d points available\n", n_quad);
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        if (n_quad == 4)
        {
            double a = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
            double b = (5.0 - sqrt(5.0)) / 20.0;

            (*q_r)[0][0] = a; (*q_r)[0][1] = b; (*q_r)[0][2] = b;
            (*q_r)[1][0] = b; (*q_r)[1][1] = a; (*q_r)[1][2] = b;
            (*q_r)[2][0] = b; (*q_r)[2][1] = b; (*q_r)[2][2] = a;
            (*q_r)[3][0] = b; (*q_r)[3][1] = b; (*q_r)[3][2] = b;

            (*q_w)[0] = 1.0 / 24.0;
            (*q_w)[1] = 1.0 / 24.0;
            (*q_w)[2] = 1.0 / 24.0;
            (*q_w)[3] = 1.0 / 24.0;
        }
        else if (n_quad == 5)
        {
            (*q_r)[0][0] = 1.0 / 2.0; (*q_r)[0][1] = 1.0 / 6.0; (*q_r)[0][2] = 1.0 / 6.0;
            (*q_r)[1][0] = 1.0 / 6.0; (*q_r)[1][1] = 1.0 / 2.0; (*q_r)[1][2] = 1.0 / 6.0;
            (*q_r)[2][0] = 1.0 / 6.0; (*q_r)[2][1] = 1.0 / 6.0; (*q_r)[2][2] = 1.0 / 2.0;
            (*q_r)[3][0] = 1.0 / 6.0; (*q_r)[3][1] = 1.0 / 6.0; (*q_r)[3][2] = 1.0 / 6.0;
            (*q_r)[4][0] = 1.0 / 4.0; (*q_r)[4][1] = 1.0 / 4.0; (*q_r)[4][2] = 1.0 / 4.0;

            (*q_w)[0] = 9.0 / 20.0;
            (*q_w)[1] = 9.0 / 20.0;
            (*q_w)[2] = 9.0 / 20.0;
            (*q_w)[3] = 9.0 / 20.0;
            (*q_w)[4] = - 4.0 / 5.0;
        }
        else
        {
            printf("No quadrature rule for %d points available\n", n_quad);
            exit(EXIT_FAILURE);
        }
    }
}

void mesh_connectivity(int ***v_coord, int ***t_map, int num_fem, int n_dim)
{
    (*v_coord) = mem_alloc_2D_int(pow(n_dim, 2), n_dim);
    (*t_map) = mem_alloc_2D_int(num_fem, n_dim + 1);

    if (n_dim == 2)
    {
        (*v_coord)[0][0] = 0; (*v_coord)[0][1] = 0;
        (*v_coord)[1][0] = 1; (*v_coord)[1][1] = 0;
        (*v_coord)[2][0] = 0; (*v_coord)[2][1] = 1;
        (*v_coord)[3][0] = 1; (*v_coord)[3][1] = 1;

        if (num_fem == 2)
        {
            (*t_map)[0][0] = 0; (*t_map)[0][1] = 1; (*t_map)[0][2] = 3;
            (*t_map)[1][0] = 0; (*t_map)[1][1] = 3; (*t_map)[1][2] = 2;
        }
        else if (num_fem == 4)
        {
            (*t_map)[0][0] = 1; (*t_map)[0][1] = 2; (*t_map)[0][2] = 0;
            (*t_map)[1][0] = 3; (*t_map)[1][1] = 0; (*t_map)[1][2] = 1;
            (*t_map)[2][0] = 0; (*t_map)[2][1] = 3; (*t_map)[2][2] = 2;
            (*t_map)[3][0] = 2; (*t_map)[3][1] = 1; (*t_map)[3][2] = 3;
        }
        else
        {
            printf("Wrong number of triangles\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        (*v_coord)[0][0] = 0; (*v_coord)[0][1] = 0; (*v_coord)[0][2] = 0;
        (*v_coord)[1][0] = 1; (*v_coord)[1][1] = 0; (*v_coord)[1][2] = 0;
        (*v_coord)[2][0] = 0; (*v_coord)[2][1] = 1; (*v_coord)[2][2] = 0;
        (*v_coord)[3][0] = 1; (*v_coord)[3][1] = 1; (*v_coord)[3][2] = 0;
        (*v_coord)[4][0] = 0; (*v_coord)[4][1] = 0; (*v_coord)[4][2] = 1;
        (*v_coord)[5][0] = 1; (*v_coord)[5][1] = 0; (*v_coord)[5][2] = 1;
        (*v_coord)[6][0] = 0; (*v_coord)[6][1] = 1; (*v_coord)[6][2] = 1;
        (*v_coord)[7][0] = 1; (*v_coord)[7][1] = 1; (*v_coord)[7][2] = 1;

        if (num_fem == 6)
        {
            (*t_map)[0][0] = 0; (*t_map)[0][1] = 2; (*t_map)[0][2] = 1; (*t_map)[0][3] = 5;
            (*t_map)[1][0] = 1; (*t_map)[1][1] = 2; (*t_map)[1][2] = 3; (*t_map)[1][3] = 5;
            (*t_map)[2][0] = 0; (*t_map)[2][1] = 4; (*t_map)[2][2] = 2; (*t_map)[2][3] = 5;
            (*t_map)[3][0] = 5; (*t_map)[3][1] = 3; (*t_map)[3][2] = 7; (*t_map)[3][3] = 2;
            (*t_map)[4][0] = 4; (*t_map)[4][1] = 5; (*t_map)[4][2] = 6; (*t_map)[4][3] = 2;
            (*t_map)[5][0] = 5; (*t_map)[5][1] = 7; (*t_map)[5][2] = 6; (*t_map)[5][3] = 2;
        }
        else if (num_fem == 8)
        {
            (*t_map)[0][0] = 0; (*t_map)[0][1] = 2; (*t_map)[0][2] = 1; (*t_map)[0][3] = 4;
            (*t_map)[1][0] = 1; (*t_map)[1][1] = 0; (*t_map)[1][2] = 3; (*t_map)[1][3] = 5;
            (*t_map)[2][0] = 2; (*t_map)[2][1] = 6; (*t_map)[2][2] = 3; (*t_map)[2][3] = 0;
            (*t_map)[3][0] = 3; (*t_map)[3][1] = 2; (*t_map)[3][2] = 7; (*t_map)[3][3] = 1;
            (*t_map)[4][0] = 4; (*t_map)[4][1] = 5; (*t_map)[4][2] = 6; (*t_map)[4][3] = 0;
            (*t_map)[5][0] = 5; (*t_map)[5][1] = 7; (*t_map)[5][2] = 4; (*t_map)[5][3] = 1;
            (*t_map)[6][0] = 6; (*t_map)[6][1] = 7; (*t_map)[6][2] = 2; (*t_map)[6][3] = 4;
            (*t_map)[7][0] = 7; (*t_map)[7][1] = 3; (*t_map)[7][2] = 6; (*t_map)[7][3] = 5;
        }
        else
        {
            printf("Wrong number of tetrahedrals\n");
            exit(EXIT_SUCCESS);
        }
    }
}

void x_map(double **x, double *r, double **x_t, int n_dim, Basis *phi)
{
    int i, d;

    for (d = 0; d < n_dim; d++)
    {
        (*x)[d] = 0.0;

        for (i = 0; i < n_dim + 1; i++)
        {
            (*x)[d] += x_t[d][i] * phi[i](r);
        }
    }
}

void J_xr_map(double ***J_xr, double *r, double **x_t, int n_dim, DBasis *dphi)
{
    int i, j, k;
    double *deriv = mem_alloc_1D_double(n_dim);

    for (i = 0; i < n_dim; i++)
    {
        for (j = 0; j < n_dim; j++)
        {
            (*J_xr)[i][j] = 0.0;

            for (k = 0; k < n_dim + 1; k++)
            {
                dphi[k](&deriv, r);

                (*J_xr)[i][j] += x_t[i][k] * deriv[j];
            }
        }
    }

    mem_free_1D_double(&deriv, n_dim);
}

// Basis functions and derivatives in 2D
double phi_2D_1(double *r) { return r[0]; }
double phi_2D_2(double *r) { return r[1]; }
double phi_2D_3(double *r) { return 1.0 - r[0] - r[1]; }
void dphi_2D_1(double **dp, double *r) { (*dp)[0] = 1.0; (*dp)[1] = 0.0; }
void dphi_2D_2(double **dp, double *r) { (*dp)[0] = 0.0; (*dp)[1] = 1.0; }
void dphi_2D_3(double **dp, double *r) { (*dp)[0] = -1.0; (*dp)[1] = -1.0; }

// Basis functions and derivatives in 3D
double phi_3D_1(double *r) { return r[0]; }
double phi_3D_2(double *r) { return r[1]; }
double phi_3D_3(double *r) { return r[2]; }
double phi_3D_4(double *r) { return 1.0 - r[0] - r[1] - r[2]; }
void dphi_3D_1(double **dp, double *r) { (*dp)[0] = 1.0; (*dp)[1] = 0.0; (*dp)[2] = 0.0; }
void dphi_3D_2(double **dp, double *r) { (*dp)[0] = 0.0; (*dp)[1] = 1.0; (*dp)[2] = 0.0; }
void dphi_3D_3(double **dp, double *r) { (*dp)[0] = 0.0; (*dp)[1] = 0.0; (*dp)[2] = 1.0; }
void dphi_3D_4(double **dp, double *r) { (*dp)[0] = -1.0; (*dp)[1] = -1.0; (*dp)[2] = -1.0; }

// Math functions
long maximum(long a, long b)
{
    return a > b ? a : b;
}

double determinant(double** A, int n)
{
    /*
     * Computes the determinant of a matrix
     */
    if (n == 2)
    {
        double d_1 = A[0][0] * A[1][1] - A[1][0] * A[0][1];

        return d_1;
    }
    else if (n == 3)
    {
        double d_1 = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]);
        double d_2 = A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]);
        double d_3 = A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);

        return d_1 - d_2 + d_3;
    }
    else if (n == 4)
    {
        double d_1 = A[0][0];
        double d_11 = A[1][1] * (A[2][2] * A[3][3] - A[3][2] * A[2][3]);
        double d_12 = A[1][2] * (A[2][1] * A[3][3] - A[3][1] * A[2][3]);
        double d_13 = A[1][3] * (A[2][1] * A[3][2] - A[3][1] * A[2][2]);

        double d_2 = A[0][1];
        double d_21 = A[1][0] * (A[2][2] * A[3][3] - A[3][2] * A[2][3]);
        double d_22 = A[1][2] * (A[2][0] * A[3][3] - A[3][0] * A[2][3]);
        double d_23 = A[1][3] * (A[2][0] * A[3][2] - A[3][0] * A[2][2]);

        double d_3 = A[0][2];
        double d_31 = A[1][0] * (A[2][1] * A[3][3] - A[3][1] * A[2][3]);
        double d_32 = A[1][1] * (A[2][0] * A[3][3] - A[3][0] * A[2][3]);
        double d_33 = A[1][3] * (A[2][0] * A[3][1] - A[3][0] * A[2][1]);

        double d_4 = A[0][3];
        double d_41 = A[1][0] * (A[2][1] * A[3][2] - A[3][1] * A[2][2]);
        double d_42 = A[1][1] * (A[2][0] * A[3][2] - A[3][0] * A[2][2]);
        double d_43 = A[1][2] * (A[2][0] * A[3][1] - A[3][0] * A[2][1]);

        return d_1 * (d_11 - d_12 + d_13) - d_2 * (d_21 - d_22 + d_23) + d_3 * (d_31 - d_32 + d_33) - d_4 * (d_41 - d_42 + d_43);
    }
    else
    {
        printf("ERROR: No determinant implementation for matrices of size > 4\n");
        exit(EXIT_FAILURE);
    }
}

void inverse(double*** inv_A, double** A, int n)
{
    /*
     * Computes the inverse of a matrix
     */
    double det_A = determinant(A, n);

    if (n == 2)
    {
        (*inv_A)[0][0] = (1.0 / det_A) * A[1][1];
        (*inv_A)[0][1] = -(1.0 / det_A) * A[0][1];
        (*inv_A)[1][0] = -(1.0 / det_A) * A[1][0];
        (*inv_A)[1][1] = (1.0 / det_A) * A[0][0];
    }
    else if (n == 3)
    {
        (*inv_A)[0][0] = (1.0 / det_A) * (A[1][1] * A[2][2] - A[2][1] * A[1][2]);
        (*inv_A)[0][1] = (1.0 / det_A) * (A[0][2] * A[2][1] - A[2][2] * A[0][1]);
        (*inv_A)[0][2] = (1.0 / det_A) * (A[0][1] * A[1][2] - A[1][1] * A[0][2]);
        (*inv_A)[1][0] = (1.0 / det_A) * (A[1][2] * A[2][0] - A[2][2] * A[1][0]);
        (*inv_A)[1][1] = (1.0 / det_A) * (A[0][0] * A[2][2] - A[2][0] * A[0][2]);
        (*inv_A)[1][2] = (1.0 / det_A) * (A[0][2] * A[1][0] - A[1][2] * A[0][0]);
        (*inv_A)[2][0] = (1.0 / det_A) * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
        (*inv_A)[2][1] = (1.0 / det_A) * (A[0][1] * A[2][0] - A[2][1] * A[0][0]);
        (*inv_A)[2][2] = (1.0 / det_A) * (A[0][0] * A[1][1] - A[1][0] * A[0][1]);
    }
    else
    {
        printf("ERROR: No inverse implementation for matrices of size > 3\n");
        exit(EXIT_FAILURE);
    }
}

// Memory management
int* mem_alloc_1D_int(int n)
{
    int *array = (int*) malloc(n * sizeof(int));

    return array;
}

long* mem_alloc_1D_long(int n)
{
    long *array = (long*) malloc(n * sizeof(long));

    return array;
}

double* mem_alloc_1D_double(int n)
{
    double *array = (double*) malloc(n * sizeof(double));

    return array;
}

int** mem_alloc_2D_int(int n, int m)
{
    int i;
    int **array = (int**) malloc(n * sizeof(int*));

    for (i = 0; i < n; i++)
        array[i] = (int*) malloc(m * sizeof(int));

    return array;
}

long** mem_alloc_2D_long(int n, int m)
{
    int i;
    long **array = (long**) malloc(n * sizeof(long*));

    for (i = 0; i < n; i++)
        array[i] = (long*) malloc(m * sizeof(long));

    return array;
}

double** mem_alloc_2D_double(int n, int m)
{
    int i;
    double **array = (double**) malloc(n * sizeof(double*));

    for (i = 0; i < n; i++)
        array[i] = (double*) malloc(m * sizeof(double));

    return array;
}

void mem_free_1D_int(int **array, int n)
{
    free((*array));
}

void mem_free_1D_long(long **array, int n)
{
    free((*array));
}

void mem_free_1D_double(double **array, int n)
{
    free((*array));
}

void mem_free_2D_int(int ***array, int n, int m)
{
    int i;

    for (i = 0; i < n; i++)
        free((*array)[i]);

    free((*array));
}

void mem_free_2D_long(long ***array, int n, int m)
{
    int i;

    for (i = 0; i < n; i++)
        free((*array)[i]);

    free((*array));
}

void mem_free_2D_double(double ***array, int n, int m)
{
    int i;

    for (i = 0; i < n; i++)
        free((*array)[i]);

    free((*array));
}
#else
void fem_amg_setup(int *precond_type_, int *meshing_type_,
                   int *n_x_, int *n_y_, int *n_z_, int *n_elem_, int *n_dim_,
                   double *x_m_, double *y_m_, double *z_m_,
                   long *glo_num_, double *pmask_, double* binv_)
{
     if (proc_id == 0)
       printf("ERROR: This requires HYPRE, please recompile!\n");
     exit(EXIT_FAILURE);
}
void fem_amg_solve(double *z, double *w)
{

     if (proc_id == 0)
       printf("ERROR: This requires HYPRE, please recompile!\n");
     exit(EXIT_FAILURE);
}
#endif
