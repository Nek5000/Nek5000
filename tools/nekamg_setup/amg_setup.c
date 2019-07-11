#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "_hypre_utilities.h"
#include "_hypre_parcsr_ls.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE.h"
#include "amg_setup.h"

/* 
    Code for performing the AMG setup for Nek5000 using the linear algebra
    library Hypre.

    Based on the original implementation of the AMG setup in Matlab by 
    James Lottes. A thorough description of the original setup can be found in
    his Ph. D. thesis "Towards Robust Algebraic Multigrid Methods for 
    Nonsymmetric Problems".
*/

int main(int argc, char *argv[])
{   
    /* Get user's input for the setup */
    int coars_strat, interp_strat, ret;
    setbuf(stdout, NULL);

    char sname[132];
    char session[132];
    printf("Enter name prefix of input file(s):\n");
    fgets(sname, sizeof sname, stdin);
    ret = sscanf(sname,"%s",&session);
    if (ret == -1) exit(1);

    /* Coarsening strategy */
    printf("Choose a coarsening method. Available options are:\n");
    //printf(" - 0: CLJP,\n");
    printf(" - 3: Ruege-Stuben,\n");
    printf(" - 6: Falgout (default),\n");
    printf(" - 8: PMIS,\n");
    printf(" - 10: HMIS,\n");
    //printf(" - 21: CGC,\n");
    //printf(" - 22: CGC-E.\n");
    printf("Choice: ");
    char scoars[30];
    fgets(scoars, sizeof scoars, stdin);
    ret = sscanf(scoars, "%d", &coars_strat);
    if (ret == -1)
    {
        coars_strat = 6; // default
    }

    /* Interpolation strategy */
    printf("Choose an interpolation method. Available options are:\n");
    printf(" - 0: classical modified interpolation,\n");
    //    printf(" - 1: LS interpolation,\n");
    //printf(" - 2: classical modified interpolation for hyperbolic PDEs,\n");
    //printf(" - 3: direct interpolation,\n");
    //printf(" - 4: multipass interpolation,\n");
    //printf(" - 5: multipass interpolation (with separation of weights),\n");
    printf(" - 6: extended + i interpolation (default),\n");
    //printf(" - 7: extended + i (if no common C neighbour) interpolation,\n");
    //printf(" - 8: standard interpolation,\n");
    //printf(" - 9: standard interpolation (with separation of weights),\n");
    //    printf(" - 10: classical block interpolation,\n");
    //    printf(" - 11: classical block interpolation with diagonalized diagonal blocks,\n");
    //printf(" - 12: FF interpolation,\n");
    //printf(" - 13: FF1 interpolation,\n");
    //printf(" - 14: extended interpolation,\n");
    printf("Choice: ");
    char sinterp[30];
    fgets(sinterp, sizeof sinterp, stdin);
    ret = sscanf(sinterp, "%d", &interp_strat);
    if (ret == -1)
    {
        interp_strat = 6; // default
    }

    /* Smoother tolerance */
    double tol;
    char stol[30];
    printf("Enter smoother tolerance [0.5]:\n");
    fgets(stol, sizeof stol, stdin);
    ret = sscanf(stol, "%lf", &tol);
    if (ret == -1)
    {
        tol = 0.5; // default
    }
    if (tol <= 0.)
    {
        printf("ERRORr: Smoother tolerance should be >0.\n");
        exit(1);
    }

     
    /* Verbose level */
    int print_level = 3;  // Print solve info + parameters

    /* Read data and convert rows and columns to integers */

    printf("Reading AMG dump files... ");
    char str1[100],str2[100],str3[100];
    sprintf(str1,"%s.amgdmp_i.dat",session);
    sprintf(str2,"%s.amgdmp_j.dat",session);
    sprintf(str3,"%s.amgdmp_p.dat",session);
    int n = filesize(str1);
    double *v   = malloc( n    * sizeof(double));
    double *Aid = malloc((n-1) * sizeof(double));
    double *Ajd = malloc((n-1) * sizeof(double));
    double *Av  = malloc((n-1) * sizeof(double));

    readfile(v,n,str1);
    memcpy(Aid, v+1, (n-1) * sizeof (double));
    readfile(v,n,str2);
    memcpy(Ajd, v+1, (n-1) * sizeof (double));
    readfile(v,n,str3);
    memcpy(Av , v+1, (n-1) * sizeof (double));
    printf("Done.\n");

    int *Ai = malloc((n-1) * sizeof (int));
    int *Aj = malloc((n-1) * sizeof (int));

    int iuppermax=0;

	/* Get matrix dimensions */
    int i;
    for (i=0;i<n-1;i++) 
    {
        Ai[i] = (int)Aid[i]-1;
        Aj[i] = (int)Ajd[i]-1;
        if (iuppermax < Ai[i]) iuppermax = Ai[i];
    }

    /* Look for zeros rows */
    int *isrownonzero = malloc((iuppermax+1) * sizeof(int));
    for (i=0;i<(iuppermax+1);i++) 
    {
        isrownonzero[i] = 0;
    }
    for (i=0;i<n-1;i++) 
    {
        if (Av[i] != 0.) isrownonzero[Ai[i]] = 1;
    }
 
    /* Data structure containing setup info */
    struct amg_setup_data *data = malloc(sizeof (struct amg_setup_data));

    /* Initialize local to global correspondence for id */
    data->id_l2g = malloc((iuppermax+1) * sizeof(int));
    int  *id_g2l = malloc((iuppermax+1) * sizeof(int));
    int c= 0;
    int ilower=0, iupper = -1;
    for (i=0;i<iuppermax+1;i++) 
    {
        if (isrownonzero[i] != 0) 
        {
            data->id_l2g[c] = i+1; // Global indices start at 1
            id_g2l[i] = c;         // Local indices start at 0
            c++;
            iupper++;
        }
    }

    /* Set local ids */
    data->idl = malloc((iupper+1) * sizeof(int));
    for (i=0;i<iupper+1;i++)
    {
        data->idl[i]=i;
    }

    /* Build Hypre IJ matrix */
    HYPRE_IJMatrix ij_matrix;
	/* Init matrix */
    HYPRE_IJMatrixCreate(0, ilower, iupper, ilower, iupper, 
                         &ij_matrix);
    HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(ij_matrix);

	/* Set matrix entries */
    for (i=0;i<n-1;i++) 
    {
		int ncols=1;
    	if (Av[i] != 0)
        {
            HYPRE_IJMatrixSetValues(ij_matrix, 1, &ncols, &(id_g2l[Ai[i]]), 
                                &(id_g2l[Aj[i]]), &Av[i]);
        }
    }

    free(v);
    free(Aid);
    free(Ajd);
    free(Ai);
    free(Aj);
    free(Av);
    free(id_g2l);

    HYPRE_IJMatrixAssemble(ij_matrix);

    /* Build Hypre CSR matrix */
    HYPRE_ParCSRMatrix A;
    HYPRE_IJMatrixGetObject(ij_matrix, (void **) &A);

    /* Create dummy rhs */
    HYPRE_ParVector b;
    HYPRE_ParVector x;

    /* AMG setup */
    HYPRE_Solver solver;
    HYPRE_BoomerAMGCreate(&solver); // Create solver
 
    /* Set parameters (See Reference Manual for more parameters) */

    HYPRE_BoomerAMGSetStrongThreshold(solver,0.5);
    if (coars_strat == 8 || coars_strat == 10) 
      HYPRE_BoomerAMGSetStrongThreshold(solver,0.25); 

    HYPRE_BoomerAMGSetPrintLevel(solver, print_level);
    HYPRE_BoomerAMGSetCoarsenType(solver, coars_strat);
    HYPRE_BoomerAMGSetInterpType(solver, interp_strat); 
    HYPRE_BoomerAMGSetMaxCoarseSize(solver, 1);
    
    /* Perform setup */
    printf("BoomerAMGSetup... ");
    HYPRE_BoomerAMGSetup(solver, A, b, x);
    printf("Done.\n");

    /* Access solver data */
    hypre_ParAMGData *amg_data = (hypre_ParAMGData*) solver; 
            // structure hypre_ParAMGData is described in parcsr_lspar_amg.h
    int numlvls = (amg_data)->num_levels; // number of levels
    hypre_ParCSRMatrix **A_array = (amg_data)->A_array; 
    hypre_ParCSRMatrix **P_array = (amg_data)->P_array;// Interpolation operator
    int **CF_marker_array        = (amg_data)->CF_marker_array;

    /* Check if last level is made of 1 element */
    if (A_array[numlvls-1]->diag->num_rows != 1)
    {
        printf("ERROR: Invalid size of last level, try different coarsening!");
	exit(1);
    }

    /* Initialize data structure */
    data->nlevels = numlvls;
    data->n     = malloc( numlvls    * sizeof (double));
    data->nnz   = malloc( numlvls    * sizeof (double));
    data->nnzf  = malloc((numlvls-1) * sizeof (double));
    data->nnzfp = malloc((numlvls-1) * sizeof (double));
    data->m     = malloc((numlvls-1) * sizeof (double));
    data->rho   = malloc((numlvls-1) * sizeof (double));
    data->idc   = malloc( numlvls    * sizeof (int*));
    data->idf   = malloc( numlvls    * sizeof (int*));
    data->D     = malloc((numlvls-1) * sizeof (double*));
    data->Af    = malloc((numlvls-1) * sizeof (hypre_CSRMatrix *));
    data->W     = malloc((numlvls-1) * sizeof (hypre_CSRMatrix *));
    data->AfP   = malloc((numlvls-1) * sizeof (hypre_CSRMatrix *));
    
    /* Exctract data and compute smoother at each level */
    int lvl;
    for (lvl=0;lvl<numlvls;lvl++)
    {
        /* Save matrix A */
        hypre_ParCSRMatrix *Alvl = A_array[lvl];
        hypre_CSRMatrix *diagA = Alvl->diag; 

        /* Update number of rows and nonzeros in data */
        int num_rows_A = diagA->num_rows;
        int num_nnzs_A = diagA->num_nonzeros;
        data->n[lvl] = num_rows_A;
        data->nnz[lvl] = num_nnzs_A;

        /* If not last level */
        if (lvl < numlvls-1)
        {
            printf("--------\nLevel %d\n--------\n", lvl+1);
            /* Extract fine and coarse variables */
            int *CFlvl = CF_marker_array[lvl]; // CFlvl[i]==1 if i is coarse

            /* Save interpolation operator W */
            hypre_ParCSRMatrix *Plvl = P_array[lvl];
            hypre_CSRMatrix *diagP = Plvl->diag;
            int nF = diagP->num_rows - diagP->num_cols;
            int nC = diagP->num_cols;

            /* W = P(F,:) */
            int *monec = malloc(nC * sizeof (int));
            int k; 
            for (k=0;k<nC;k++) monec[k] = -1; 
            hypre_CSRMatrix *diagW;
            sub_mat(&diagW, Plvl->diag, CFlvl, monec);// Extract rows "manually"
            free(monec);
            data->W[lvl] = diagW;

            /* Save idf, idc and idl */
            data->idf[lvl] = malloc(nF * sizeof (int));
            data->idc[lvl] = malloc(nC * sizeof (int));
            int *idl;
            if (lvl == 0) idl = data->idl;
            else          idl = data->idc[lvl-1];

            /* Save ids for coarse and fine variables */
            int *indicesF = malloc((nF+1) * sizeof (int));
            indicesF[0] = nF;
            int j, cc=0, cf=0;
            for (j=0;j<num_rows_A;j++) 
            {
                if (CFlvl[j] == 1) 
                {
                    (data->idc[lvl])[cc] = idl[j];
                    cc++;
                }
                else
                {
                    (data->idf[lvl])[cf] = idl[j];
                    indicesF[cf+1] = j;
                    cf++;
                }
            }

            /* Extract submatrices
               - Aff = A(F,F) 
               - Afc = A(F,C) */
            hypre_ParCSRMatrix **submatricesAf =
                    malloc(4*sizeof(hypre_ParCSRMatrix*));
            hypre_ParCSRMatrixExtractSubmatrices(Alvl, 
                                                 indicesF, &submatricesAf);
            hypre_ParCSRMatrix *Aff = (submatricesAf)[0];
            hypre_ParCSRMatrix *Afc = (submatricesAf)[1];
            hypre_CSRMatrix *diagAff = Aff->diag;
            diagAff->num_rows = nF; // Correct bug
            diagAff->num_cols = nF; // Correct bug
            hypre_CSRMatrix *diagAfc = Afc->diag;
            free(indicesF);

            /* AfP = Af*W + Afc */
            hypre_CSRMatrix *diagAfW;
            diagAfW = hypre_CSRMatrixMultiply(diagAff, diagW);
            hypre_CSRMatrix *diagAfP;
            diagAfP = hypre_CSRMatrixAdd(diagAfW, diagAfc);
            hypre_CSRMatrixDestroy(diagAfW); // destroy useless matrix

            /* Save matrices Af and AfP */
            data->Af[lvl] = diagAff;
            data->AfP[lvl] = diagAfP;

            /* Set nnzs and nnzfp in data */
            int nnzAff = diagAff->num_nonzeros;
            int nnzAfP = diagAfP->num_nonzeros;
            data->nnzf[lvl]  = (double)nnzAff;
            data->nnzfp[lvl] = (double)nnzAfP;

/* Smoother ----------------------------------------------------------------- */
            printf("Computing diagonal smoother... ");
            double *af = malloc(nnzAff * sizeof (double)); 
            memcpy(af, diagAff->data, nnzAff * sizeof (double));
            double *af2 = af; // af2 = Af.*Af
            vv_op(af2, af2, nnzAff, ewmult); // af2 = Af.*Af
            double *s = malloc(nF * sizeof (double)); // s = 1./sum(Af.*Af)
            int i;
            for (i=0; i<nF; i++)
            {
                int js = diagAff->i[i];
                int je = diagAff->i[i+1]; 
                int nsum = je-js;

                s[i] = array_op(af2, nsum, sum_op); // s = sum(af2)
                af2 += nsum;
            }
            free(af);
            array_op(s, nF, minv_op); // s = 1./s

            double *D = malloc(nF * sizeof (double)); // D = diag(Af)' .* s
            diag(D, diagAff);
            vv_op(D, s, nF, ewmult);
            free(s);
            printf("Done.\n");

            if (nF >= 2)
            {
                printf("Running Lanczos... ");
                double *Dh = malloc(nF * sizeof (double)); // Dh = sqrt(D)
                memcpy(Dh, D, nF * sizeof (double));
                array_op(Dh, nF, sqrt_op);

                hypre_CSRMatrix *diagDhAfDh;
                diagDhAfDh = hypre_CSRMatrixCreate(nF, nF, nnzAff);
                hypre_CSRMatrixInitialize(diagDhAfDh);

                hypre_CSRMatrixCopy(diagAff, diagDhAfDh, 1);

                dxm(diagDhAfDh, Dh, dmult); // DhAfDh = Dh*Af
                dxm(diagDhAfDh, Dh, multd); // DhAfDh = Dh*Af*Dh

                double *lambda; // vector of eigenvalues
                int k = lanczos(&lambda, diagDhAfDh); // k=number of eigenvalues

                printf("k = %d, ", (int)k);
                printf("[lambda[0], lambda[%d]] = [%lf, %lf],",(int)k-1,
                        lambda[0], lambda[k-1]);

                double a = lambda[0]; // First and last eigenvalues
                double b = lambda[k-1];

                ar_scal_op(D, 2./(a+b), nF, mult_op);
                
                data->D[lvl] = malloc(nF * sizeof (double)); 
                memcpy(data->D[lvl], D, nF * sizeof (double));

                double rho = (b-a)/(b+a);
                data->rho[lvl] = rho;
                printf(" rho = %lf\n", rho);
                
                double m, c; 
                double gamma2 = 1. - sqrt(1. - tol);
                chebsim(&m, &c, rho, gamma2);
                data->m[lvl] = m;
                printf("Chebyshev smoother iterations: %d\tContraction: %lf\n",
                        (int)m,c);
                if((int)m > 3){
                  printf("ERROR: Smoother iterations too large, try different coarsening!\n");
                  exit(1);
                }

                free(Dh);    
                free(lambda);
                hypre_CSRMatrixDestroy(diagDhAfDh);
            }
            else
            {
                data->D[lvl] = malloc(nF * sizeof (double)); 
                memcpy(data->D[lvl], D, nF * sizeof (double));     

                data->rho[lvl] = 0;
                data->m[lvl] = 1;
            }
/* -------------------------------------------------------------------------- */ 
        }

        if (lvl == (numlvls - 1))
        {
            data->alast = diagA->data[0];
            if (data->alast <= 1e-9) data->nullspace = 1;
            else                        data->nullspace = 0;
            printf("Nullspace = %d\n", data->nullspace);
        }
    }

    printf("Setup finished... Exporting data.\n");
    amg_export(data,session);

    /* Destroy matrix ij */
    HYPRE_IJMatrixDestroy(ij_matrix);

    /* Destroy solver */
    HYPRE_BoomerAMGDestroy(solver);

    return 1;
}

/* 
    Extract sub-matrix: 
*/
static void sub_mat(hypre_CSRMatrix **subA, const hypre_CSRMatrix *A, 
    const int* idr, const int *idc)
{
    int rn = A->num_rows, cn = A->num_cols;
    int *row_off = A->i, *col = A->j;
    double *a = A->data;
    int subrn = 0, subcn = 0, subncol = 0;
    int i, j;
    /* Compute number of rows and of non-zero elements for sub-matrix */
    for (i=0; i<rn; i++)
    {
        if (idr[i] < 0)
        {
            subrn += 1;
            int je=row_off[i+1]; 
            for(j=row_off[i]; j<je; j++)
            {
                if (idc[col[j]] < 0)
                {
                    subncol += 1;
                }
            }
        }
    }

    /* Compute cn for sub-matrix */
    int c = 0; // counter
    int *g2lcol = malloc(cn * sizeof (int)); // correspondence between local and 
                                             // global column ids
    for (i=0; i<cn; i++)
    {
        if (idc[i] < 0)
        {
            subcn += 1;
            g2lcol[i] = c;
            c++;
        }
        else
        {
            g2lcol[i] = -1; // dummy
        }
    }

    /* Initialize and build sub-matrix */
    *subA = hypre_CSRMatrixCreate(subrn, subcn, subncol);
    hypre_CSRMatrix *subA_ptr = *subA;
    hypre_CSRMatrixInitialize(subA_ptr);

    int *subrow_off = subA_ptr->i;
    int *subcol     = subA_ptr->j;
    double *suba    = subA_ptr->data;

    *subrow_off++ = 0;
    int roffset = 0;

    for (i=0; i<rn; i++)
    {
        if (idr[i] < 0)
        {
            int je=row_off[i+1]; 
            for(j=row_off[i]; j<je; j++)
            {
                if (idc[col[j]] < 0)
                {
                    roffset += 1;
                    *subcol++ = g2lcol[col[j]];
                    *suba++   = a[j];
                }
            }
            *subrow_off++ = roffset;
        }
    }

    free(g2lcol);
}

/*
    Export data from the AMG setup to correct format.
*/
static void amg_export(const struct amg_setup_data *data,char *session)
{
    int nlevels = data->nlevels;
    int n = data->n[0];
    
    /* Find what is the last level at which each variable appears */
    int *lvl = malloc(n * sizeof (int));
    int i, j;

    for (i=0;i<n;i++) lvl[i] = 1;
    for (i=0;i<nlevels-1;i++)
    {
        int nl = data->n[i+1];
        int *idl = data->idc[i];
        for (j=0;j<nl;j++)
        {
            lvl[idl[j]] += 1;
        }
    }

    /* Reorder diagonal smoother (except for last level) */
    double *dvec = malloc(n * sizeof (double));
    for (i=0;i<nlevels-1;i++)
    {
        int nl = data->n[i]-data->n[i+1];
        int *idl = data->idf[i];
        double *Dl = data->D[i];
        for (j=0;j<nl;j++)
        {
            dvec[idl[j]] = Dl[j];
        }
    }

    /* Set smoother at last level */
    int k = (data->idc[nlevels-2][0]);
    if (data->nullspace != 0) 
    {
        dvec[k] = 0.;
    }
    else 
    {
        double a = data->alast;
        dvec[k] = 1./a;
    }

    /* Save matrices */
    char str1[132],str2[132],str3[132],str4[132];
      sprintf(str1,"%s.amg_W.dat",session);
      sprintf(str2,"%s.amg_AfP.dat",session);
      sprintf(str3,"%s.amg_Aff.dat",session);
      sprintf(str4,"%s.amg.dat",session);
    int *W_len = malloc(n * sizeof (int));
    savemats(W_len, n, nlevels-1, lvl, data->idc, data->id_l2g, data->W,
             str1);        

    int *AfP_len = malloc(n * sizeof (int));
    savemats(AfP_len, n, nlevels-1, lvl, data->idc, data->id_l2g, data->AfP, 
             str2);

    int *Aff_len = malloc(n * sizeof (int));
    savemats(Aff_len, n, nlevels-1, lvl, data->idf, data->id_l2g, data->Af, 
             str3);

    /* Save vector */
    savevec(nlevels, data, n, lvl, W_len, AfP_len, Aff_len, dvec,str4);

    /* Free allocated memory */
    free(W_len);
    free(AfP_len);
    free(Aff_len);
    free(dvec);
    free(lvl);
}

/*
    Function to save matrices
*/
static void savemats(int *len, const int n, const int nl, const int *lvl, 
    int **id, const int *id_l2g, hypre_CSRMatrix **mat, const char *filename)
{
    const double magic = 3.14159;
    FILE *f = fopen(filename,"w");
    if (f == NULL)
    {
        printf("\nERROR: File %s not found!\n", filename);
        exit(1);
    }
    int max=0;
    int i;
    int *row;
    double *buf;
    fwrite(&magic,sizeof (double),1,f);
    for(i=0;i<nl;++i) {int l=max_row_nnz(mat[i]); if(l>max) max=l; }
    buf = malloc(2*max * sizeof (double));
    row = malloc(nl    * sizeof (int));
    for(i=0;i<nl;++i) row[i]=0;
    for(i=0;i<n;++i) 
    {
        int l = lvl[i]-1;
        hypre_CSRMatrix *M;
        int j,k,kb,ke;
        double *p;
        if(l>nl) { printf("level out of bounds\n"); continue; }
        if(l==nl) { len[i]=0; continue; }
        M = mat[l];
        j = row[l]++;
        if(j>=M->num_rows) { printf("row out of bounds\n"); continue; }
        kb=M->i[j],ke=M->i[j+1];
        p = buf;
        for(k=kb;k!=ke;++k) *p++ = id_l2g[id[l][M->j[k]]], *p++ = M->data[k];
        len[i] = ke-kb;
        fwrite(buf,sizeof (double),(double)(2*(ke-kb)),f);
    }
    for(i=0;i<nl;++i) 
    {
        if(row[i]!=mat[i]->num_rows) printf("matrices not exhausted\n");
    }
    free(row);
    free(buf);
    fclose(f);

    printf("%d matrices written to %s\n",(int)nl, filename);
}

/* 
    Returns max number of non zero elements on a row
*/
static int max_row_nnz(const hypre_CSRMatrix *mat)
{
    int n=mat->num_rows,max=0;     

    int i;
    for(i=0;i<n;++i) 
    {
        int l=mat->i[i+1]-mat->i[i];
        if(l>max) max=l;
    }
    return max;
}

/*
    Function to save vectors
*/
static void savevec(const int nl, const struct amg_setup_data *data, 
    const int n, const int *lvl, const int *W_len, const int *AfP_len, 
    const int *Aff_len, const double *dvec, const const char *filename)
{
    int ntot = 2+2*nl+6*n;
    double *q = malloc(ntot * sizeof (double));
    const double magic = 3.14159;
    const double stamp = 2.01;

    q[0] = magic;
    q[1] = stamp;
    q[2] = (double)nl;
    memcpy(&q[3]       , data->m  , sizeof(double)*(nl-1));
    memcpy(&q[3+(nl-1)], data->rho, sizeof(double)*(nl-1));
    q[3+2*(nl-1)] = (double)n;

    int i;
    int qi = 2+2*nl;
    for (i=0;i<n;i++)
    {
            q[qi++] = (double)data->id_l2g[data->idl[i]];
            q[qi++] = (double)lvl[i];
            q[qi++] = (double)W_len[i];
            q[qi++] = (double)AfP_len[i];
            q[qi++] = (double)Aff_len[i];
            q[qi++] = dvec[i];
    }

    FILE *f = fopen(filename,"w");
    if (f == NULL)
    {
        printf("\nERROR: File %s not found!\n", filename);
        exit(1);
    }
    fwrite(q,sizeof(double),ntot,f);
    fclose(f);
    free(q);

    printf("%d doubles written to %s\n", (int)ntot, filename);
}

/*
    Chebsim: computes number of iteration and contraction factor for the    
    Chebyshev relaxation.
*/
static void chebsim(double *m, double *c, const double rho, const double tol)
{
    double alpha = 0.25*rho*rho;
    *m = 1;
    double cp = 1;
    *c = rho;
    double gamma = 1;
    double d, cn;
    
    while (*c > tol)
    {
        *m += 1;
        d = alpha*(1+gamma);
        gamma = d/(1-d);
        cn = (1+gamma)*rho*(*c) - gamma*cp;
        cp = *c;
        *c = cn;
    }
}

/* 
    Compute eigenvalues by Lanczos algorithm
*/
static int lanczos(double **lambda, hypre_CSRMatrix *A)
{
    int rn = A->num_rows;
    int cn = A->num_cols;
    double *r = malloc(rn * sizeof (double));
    
    int i;
    for (i=0; i<rn; i++)
    {
        r[i] = (double)rand() / (double)RAND_MAX;
    }
    
    int kmax = 299; // Fixed length for max value of k
    *lambda = malloc(kmax * sizeof (double));
    double *l = *lambda;
    double *y = malloc(kmax * sizeof (double));
    double *a = malloc(kmax * sizeof (double));
    double *b = malloc(kmax * sizeof (double));
    double *d = malloc((kmax+1) * sizeof (double));
    double *v = malloc(kmax * sizeof (double));

    double beta = array_op(r, rn, norm2_op);
    double beta2 = beta*beta;

    int k = 0; 
    double change = 0.0;

    double *eye = malloc(rn * sizeof (double));
    init_array(eye, rn, 1.);

    /* Frobenius norm */
    hypre_CSRMatrix *Ae;
    dpm(&Ae, A, eye, dsub);
    /* A assumed to be real ==> Frobenius norm is 2-norm of A->a */
    double fronorm = array_op(Ae->data, Ae->i[rn], norm2_op);

    hypre_CSRMatrixDestroy(Ae);

    if (fronorm < 1e-11)
    {
        l[0] = 1;
        l[1] = 1;
        y[0] = 0;
        y[1] = 0;
        k = 2;
        change = 0.0;
    }
    else
    {
        change = 1.0;
    }

    if (rn == 1)
    {
        double A00 = A->data[0];

        l[0] = A00;
        l[1] = A00;
        y[0] = 0;
        y[1] = 0;
        k = 2;
        change = 0.0; 
    }

    double *qk = malloc(rn * sizeof (double));
    init_array(qk, cn, 0.);
    double *qkm1 = malloc(rn * sizeof (double));
    double *alphaqk = malloc(rn * sizeof (double)); // alpha * qk vector
    double *Aqk = malloc(rn * sizeof (double));
    int na = 0, nb = 0;

    hypre_Vector *Aqk_hyper, *qk_hyper;
    qk_hyper = hypre_SeqVectorCreate(rn);
    hypre_SeqVectorInitialize(qk_hyper);

    Aqk_hyper = hypre_SeqVectorCreate(rn);
    hypre_SeqVectorInitialize(Aqk_hyper);

    while (k < kmax && ( change > 1e-5 || y[0] > 1e-3 || y[k-1] > 1e-3))
    {
        k++;
        memcpy(qkm1, qk, rn*sizeof(double)); // qkm1 = qk
        memcpy(qk, r, rn*sizeof(double)); // qk = r/beta
        ar_scal_op(qk, 1./beta, rn, mult_op);
        qk_hyper->data = qk; 
        hypre_CSRMatrixMatvec( 1, A, qk_hyper, 0, Aqk_hyper); // Aqk = A*qk
        Aqk = Aqk_hyper->data;
        double alpha = vv_dot(qk, Aqk, rn);  // alpha = qk'*Aqk
        a[na++] = alpha;//a = [a; alpha];
        /* r = Aqk - alpha*qk - beta*qkm1 */
        memcpy(alphaqk, qk, rn*sizeof(double)); // alphaqk = qk
        ar_scal_op(alphaqk, alpha, rn, mult_op); // alphaqk = alpha*qk
        ar_scal_op(qkm1, beta, rn, mult_op); // qkm1 = beta*qkm1

        memcpy(r, Aqk, rn*sizeof(double)); // r = Aqk
        vv_op(r, alphaqk, rn, minus); // r = Aqk - alpha*qk
        vv_op(r, qkm1, rn, minus); // r = Aqk - alpha*qk - beta*qkm1
        
        if (k == 1)
        {
            l[0] = alpha;
            y[0] = 1;
        }
        else
        {
            double l0 = l[0];
            double lkm2 = l[k-2];
            d[0] = 0;
            for (i=1; i<k; i++) d[i] = l[i-1];
            d[k] = 0;
            v[0] = alpha;
            for (i=1; i<k; i++) v[i] = beta*y[i-1]; // y assumed to be real !!!
            tdeig(l, y, d, v, k-1);
            change = fabs(l0 - l[0]) + fabs(lkm2 - l[k-1]);
        }        

        beta = array_op(r, rn, norm2_op);
        beta2 = beta*beta;
        beta = sqrt(beta2);
        b[nb++] = beta;

        if (beta == 0) {break;}
    }

    int n = 0;

    for (i=0; i<k; i++)
    {
        if (y[i] < 0.01)
        {
            (*lambda)[n++] = l[i];
        }
    }  

    /* Free allocated memory */
    free(r);
    free(eye);
    free(qkm1);
    free(alphaqk);
    free(y);
    free(a);
    free(b);
    free(d);
    free(v);
    hypre_SeqVectorDestroy(qk_hyper);
    hypre_SeqVectorDestroy(Aqk_hyper);

    return n;
}

/* 
    Minimizes cancellation error (but not round-off ...) 
*/
static double sum_3(const double a, const double b, const double c)
{
  if     ( (a>=0 && b>=0) || (a<=0 && b<=0) ) return (a+b)+c;
  else if( (a>=0 && c>=0) || (a<=0 && c<=0) ) return (a+c)+b;
  else return a+(b+c);
}

/* 
    Solve 
             c
          - --- + b + a x == 0        with sign(x) = sign
             x
*/
static double rat_root(const double a, const double b, const double c,
                       const double sign)
{
  double bh = (fabs(b) + sqrt(b*b + 4*a*c))/2;
  return sign * (b*sign <= 0 ? bh/a : c/bh);
}

/*
    Find d[ri] <= lambda <= d[ri+1]
    such that 0 = lambda - v[0] + \sum_i^n v[i]^2 / (d[i] - lambda)
*/
#define EPS (128*DBL_EPSILON)
static double sec_root(double *y, const double *d, const double *v,
                       const int ri, const int n)
{
  double dl = d[ri], dr = d[ri+1], L = dr-dl;
  double x0l = L/2, x0r = -L/2;
  int i;
  double al, ar, bln, blp, brn, brp, cl, cr;
  double fn, fp, lambda0, lambda;
  double tol = L;
  if(fabs(dl)>tol) tol=fabs(dl);
  if(fabs(dr)>tol) tol=fabs(dr);
  tol *= EPS;
  for(;;) {
    if(fabs(x0l)==0 || x0l < 0) { *y=0; return dl; }
    if(fabs(x0r)==0 || x0r > 0) { *y=0; return dr; }
    lambda0 = fabs(x0l) < fabs(x0r) ? dl + x0l : dr + x0r;
    al = ar = cl = cr = bln = blp = brn = brp = 0;
    fn = fp = 0;
    for(i=1;i<=ri;++i) {
      double den = (d[i]-dl)-x0l;
      double fac = v[i]/den;
      double num = sum_3(d[i],-dr,-2*x0r);
      fn += v[i]*fac;
      fac *= fac;
      ar += fac;
      if(num > 0) brp += fac*num; else brn += fac*num;
      bln += fac*(d[i]-dl);
      cl  += fac*x0l*x0l;
    }
    for(i=ri+1;i<=n;++i) {
      double den = (d[i]-dr)-x0r;
      double fac = v[i]/den;
      double num = sum_3(d[i],-dl,-2*x0l);
      fp += v[i]*fac;
      fac *= fac;
      al += fac;
      if(num > 0) blp += fac*num; else bln += fac*num;
      brp += fac*(d[i]-dr);
      cr  += fac*x0r*x0r;
    }
    if(lambda0>0) fp+=lambda0; else fn+=lambda0;
    if(v[0]<0) fp-=v[0],blp-=v[0],brp-=v[0];
          else fn-=v[0],bln-=v[0],brn-=v[0];
    if(fp+fn > 0) { /* go left */
      x0l = rat_root(1+al,sum_3(dl,blp,bln),cl,1);
      lambda = dl + x0l;
      x0r = x0l - L;
    } else { /* go right */
      x0r = rat_root(1+ar,sum_3(dr,brp,brn),cr,-1);
      lambda = dr + x0r;
      x0l = x0r + L;
    }
    if( fabs(lambda-lambda0) < tol ) {
      double ty=0, fac;
      for(i=1;i<=ri;++i) fac = v[i]/((d[i]-dl)-x0l), ty += fac*fac;
      for(i=ri+1;i<=n;++i) fac = v[i]/((d[i]-dr)-x0r), ty += fac*fac;
      *y = 1/sqrt(1+ty);
      return lambda;
    }
  }
}

/*
    Find eigenvalues
*/
static void tdeig(double *lambda, double *y, double *d, const double *v,
                  const int n)
{
  int i;
  double v1norm = 0, min=v[0], max=v[0];
  for(i=1;i<=n;++i) {
    double vi = fabs(v[i]), a=d[i]-vi, b=d[i]+vi;
    v1norm += vi;
    if(a<min) min=a;
    if(b>max) max=b;
  }
  d[0]   = v[0] - v1norm < min ? v[0] - v1norm : min;
  d[n+1] = v[0] + v1norm > max ? v[0] + v1norm : max;
  for(i=0;i<=n;++i) lambda[i] = sec_root(&y[i],d,v,i,n);
}

/*
    Dot product between two vectors
*/
static double vv_dot(const double *a, const double *b, const int n)
{
    double r = 0;
    int i;
    for (i=0;i<n;i++)
    {
        r += (*a++)*(*b++);
    }

    return r;
}

/* 
    Array initialization
*/
static void init_array(double *a, const int n, const double v)
{
    int i;
    for (i=0;i<n;i++)
    {
        *a++ = v;
    }
}

/*
    Operations between an array and a scalar
*/
static void ar_scal_op(double *a, const double scal, const int n, 
    enum ar_scal_ops op)
{
    int i;
    switch(op)
    {
        case(mult_op)  :  for (i=0;i<n;i++) *a = (*a)*scal, a++; break;
        case(add_op)   :  for (i=0;i<n;i++) *a = (*a)+scal, a++; break;
    }
}

/* 
    Extract diagonal from square sparse matrix 
*/
static void diag(double *D, const hypre_CSRMatrix *A)
{
    int i; 
    int rn=A->num_rows;
    int *row_off = A->i, *col = A->j;
    double *a = A->data;
    
    for(i=0;i<rn;++i) 
    {
        int j; 
        const int jb=row_off[i], je=row_off[i+1];
        if (jb == je) *D++ = 0.; // A(i,:) is empty -> A(i,i)=0
        for(j=jb; j<je; ++j) 
        {  
            if (col[j] == i) 
            {
                *D++ = a[j]; 
                break;
            }
            if (j == je-1) *D++ = 0.; // A(i,i) is not present -> A(i,i)=0
        }
    }
}

/*
    Addition/subtraction between csr and diagonal matrices
*/
static void dpm(hypre_CSRMatrix **AD, hypre_CSRMatrix *A, const double *D,
    enum dpm_ops op)   
{
    int rn = A->num_rows;

    /* Build Hypre IJ diagonal matrix */
    HYPRE_IJMatrix D_ij;
	/* Init matrix */
    HYPRE_IJMatrixCreate(0, 0, rn-1, 0, rn-1, &D_ij);
    HYPRE_IJMatrixSetObjectType(D_ij, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(D_ij);

    /* Set matrix entries */
    int i;
    switch(op)
    {
        case(dadd): for (i=0;i<rn;i++) 
                    {
		                int ncols=1;
    	                if (D[i] != 0)
                        {
                            HYPRE_IJMatrixSetValues(D_ij, 1, &ncols, &i, &i, 
                                                    &D[i]);
                        }
                    }
                    break;
        case(dsub): for (i=0;i<rn;i++) 
                    {
		                int ncols=1;
    	                if (D[i] != 0)
                        {
                            double d = -D[i];
                            HYPRE_IJMatrixSetValues(D_ij, 1, &ncols, &i, &i, 
                                                    &d);
                        }
                    }
                    break;
    }
    HYPRE_IJMatrixAssemble(D_ij);

    /* Build diagonal Hypre CSR matrix */
    HYPRE_ParCSRMatrix ParD;
    HYPRE_IJMatrixGetObject(D_ij, (void **) &ParD);
    hypre_ParCSRMatrix *T = ParD;
    hypre_CSRMatrix *Dcsr = T->diag;
    /* Add A and D */
    *AD = hypre_CSRMatrixAdd(A, Dcsr);
}

/*
    Addition/subtraction between csr and diagonal matrices
*/
static void dxm(hypre_CSRMatrix *A, const double *D, enum dxm_ops op)
{
    int i;
    int rn = A->num_rows;
    int *row_off = A->i, *col = A->j;
    double *a = A->data;

    switch(op)
    {
        case(dmult):  for(i=0;i<rn;++i) 
                      {
                          int j; 
                          const int jb=row_off[i], je=row_off[i+1];
                          for(j=jb; j<je; ++j) 
                          {       
                              a[j] = a[j]*D[i];
                          }
                      }
                      break;
        case(multd):  for(i=0;i<rn;++i) 
                      {
                          int j; 
                          const int jb=row_off[i], je=row_off[i+1];
                          for(j=jb; j<je; ++j) 
                          {       
                              a[j] = a[j]*D[col[j]];
                          }
                      }
                      break;
    }
}

/*
    Array operations
*/
static double array_op(double *a, const int n, enum array_ops op)
{
    double r = 0;
    int i;
    switch(op)
    {
        case(abs_op):  for (i=0;i<n;i++) *a = fabs(*a), a++; break;
        case(sqrt_op): for (i=0;i<n;i++) *a = sqrt(*a), a++; break;
        case(minv_op): for (i=0;i<n;i++) *a = 1./(*a), a++; break;
        case(sqr_op):  for (i=0;i<n;i++) *a = (*a)*(*a), a++; break;
        case(sum_op)  :  for (i=0;i<n;i++) r += *a++; break;
        case(norm2_op):  for (i=0;i<n;i++)
                         { 
                             r += (*a)*(*a);
                             a++; 
                         }
                         r = sqrt(r); 
                         break;
    }

    return r;
}

/*
    Vector-vector operations
    OUTPUT:
    - a: array on which to perform operation (a[i] = a[i] (op) b[i])
    INPUT:
    - a: first array
    - b: second array
    - n: length of the array
    - op: operation to be performed
*/
static void vv_op(double *a, const double *b, const int n, enum vv_ops op)
{
    int i;
    switch(op)
    {
        case(plus):   for (i=0;i<n;i++) *a = *a + *b, a++, b++; break;
        case(minus):  for (i=0;i<n;i++) *a = *a - *b, a++, b++; break;
        case(ewmult): for (i=0;i<n;i++) *a = *a * (*b), a++, b++; break;
        case(ewdiv):  for (i=0;i<n;i++) *a = *a / (*b), a++, b++; break;
    }
}

/*
    Swap bytes for big endian / small endian consistency
*/
static double byteswap(double x)
{
#define N sizeof(double)
  char buf[N]; char t;
  memcpy(buf,&x,N);
#define SWAP(i) if(N>2*(i)+1) t=buf[i],buf[i]=buf[N-1-(i)],buf[N-1-(i)]=t
#define SWAP2(i) SWAP(i); SWAP((i)+1)
#define SWAP4(i) SWAP2(i); SWAP2((i)+2)
#define SWAP8(i) SWAP4(i); SWAP4((i)+4)
#define SWAP16(i) SWAP8(i); SWAP8((i)+8)
  SWAP16(0);
#undef SWAP
#undef SWAP2
#undef SWAP4
#undef SWAP8
#undef SWAP16
  memcpy(&x,buf,N);
  return x;
#undef N
}

/*
    Get size file in number of doubles
*/
static long filesize(const char *name) 
{
  long n;
  FILE *f = fopen(name,"r");
  if (f == NULL)
  {
      printf("\nERROR: File %s not found!\n", name);
      exit(1);
  }
  fseek(f,0,SEEK_END);
  n = ftell(f)/sizeof(double);
  fclose(f);
  return n;
}

/*
    Read input files
*/
static long readfile(double *data, long max, const char *name)
{
  const double magic = 3.14159;
  long n;
  FILE *f = fopen(name,"r");
  if (f == NULL)
  {
      printf("\nERROR: File %s not found!\n", name);
      exit(1);
  }  
  
  fseek(f,0,SEEK_END);
  n = ftell(f)/sizeof(double);
  if(n>max) printf("file longer than expected"),n=max;
  fseek(f,0,SEEK_SET);
  fread(data,sizeof(double),n,f);
  fclose(f);
  if(n>0 && fabs(data[0]-magic)>0.000001) {
    long i;
    /*printf("swapping byte order"); */
    if(fabs(byteswap(data[0])-magic)>0.000001) {
      printf("magic number for endian test not found");
    } else
      for(i=0;i<n;++i) data[i]=byteswap(data[i]);
  }
  return n;
}

/* 
    Print CSR matrix
    For debugging purposes only
*/
void print_csrmat(hypre_CSRMatrix *A)
{
    int num_rows = A->num_rows;
    int num_cols = A->num_cols;
    int num_nonzeros = A->num_nonzeros;

    printf("PRINT MATRIX\n"); 
    printf(" Dimensions:\n");
    printf("  num_rows = %d, num_cols = %d, num_nonzeros = %d\n", num_rows,
           num_cols, num_nonzeros);

    int *row_off = A->i;
    int *col = A->j;
    double *val = A->data;

    printf(" Data:\n"); 
    int ip, jp, jpe, jps;
    for (ip=0;ip<num_rows;ip++)
    {
        jps = row_off[ip];
        jpe = row_off[ip+1];   
        printf("   js = %d, je = %d\n", jps, jpe);      
        for (jp=jps;jp<jpe;jp++)
        {
            printf("  A(%d,%d) = %lf\n", ip+1, col[jp]+1, val[jp]); 
        }
    }
}
