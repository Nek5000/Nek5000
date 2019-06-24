#ifndef AMG_HYPRE
#define AMG_HYPRE

/*******************************************************************************
* Matrix in compressed sparse row (CSR) format
*******************************************************************************/
struct csr_mat 
{
    int rn, cn, *row_off, *col; // number of rows, number of columns,
        // offset for the rows, column number
    double *a; // data
};

/*******************************************************************************
* Main data structure
*******************************************************************************/
struct amg_setup_data 
{
    double *n; // number of variables at each level
    double *nnz; // number of non zeros of matrix A at each level
    double *nnzf; // number of non zeros of matrix Af at each level
    double *nnzfp; // number of non zeros of matrix AfP at each level
    double *m; // number of Chebyshev iterations at each level
    double *rho; // contraction for one Chebyshev iteration at each level
    int *idl; // local id (ranges from 0 to # of nonzero lines-1)
    int *id_l2g; // local to global id (ranges from 1 to max of global id)
    int **idc; // coarse subset of idl
    int **idf; // fine subset of idl
    double **D; // diagonal smoother at each level
    double alast; // A at the last level (one variable only)
    hypre_CSRMatrix **Af; // Af (fine rows, fine columns) at each level
    hypre_CSRMatrix **W; // interpolation operator at each level
    hypre_CSRMatrix **AfP; // AfP at each level
    int nlevels; // number of levels
    int nullspace; // 0 if no nullspace / 1 otherwise
};

/*******************************************************************************
* Algebraic functions
*******************************************************************************/
/* 
    Compute eigenvalues by Lanczos algorithm
    OUTPUT:
    - n = number of eigenvalues computed
    - lambda: array of eigenvalues
    INPUT:
    - A: matrix
*/
static int lanczos(double **lambda, hypre_CSRMatrix *A);

/* 
    Extract sub-matrix: 
    OUTPUT:    
    - subA = A(idr, idc)
    INPUT:
    - A: original matrix of dim. rn x cn
    - idr: if idr[i] < 0, i=0...rn-1 -> keep row i in subA
    - idc: if idc[j] < 0, j=0...cn-1 -> keep column j in subA
*/
static void sub_mat(hypre_CSRMatrix **subA, const hypre_CSRMatrix *A, 
    const int* idr, const int *idc);

/* 
    Extract diagonal from square sparse matrix 
    OUTPUT:
    - D: array of double containing diagonal of A
    INPUT:
    - A: matrix from which to extract diagonal
*/
static void diag(double *D, const hypre_CSRMatrix *A);

/*
    Multiplication between csr and diagonal matrices
    OUTPUT:
    - A: A=A*D or A=D*A
    INPUT:
    - A: matrix
    - D: array of doubles representing the diagonal
    - op: operation to be performed (A*D or D*A)
*/
enum dxm_ops {dmult, multd}; // A = {D*A, A*D}
static void dxm(hypre_CSRMatrix *A, const double *D, enum dxm_ops op);

/*
    Addition/subtraction between csr and diagonal matrices
    OUTPUT:
    - AD: AD = A +/- D
    INPUT:
    - A: matrix
    - D: array of doubles representing the diagonal
    - op: operation to be performed (addition of subtraction)
*/
enum dpm_ops {dadd, dsub}; // A = {A+D, A-D}
static void dpm(hypre_CSRMatrix **AD, hypre_CSRMatrix *A, const double *D,
    enum dpm_ops op);

/*
    tdeig: find the eigenvalues of
  
      d[1]           v[1]
             d[2]      v[2]
              d[n] v[n]
      v[1] v[2] v[n] v[0]
  
      sets d[0], d[n+1] to Gershgorin bounds
  
      also gives (n+1)th component of each orthonormal eigenvector in y
*/
static void tdeig(double *lambda, double *y, double *d, const double *v,
                  const int n);

/*
    Find d[ri] <= lambda <= d[ri+1]
    such that 0 = lambda - v[0] + \sum_i^n v[i]^2 / (d[i] - lambda)
*/
static double sec_root(double *y, const double *d, const double *v,
                       const int ri, const int n);

/* 
    Solve 
             c
          - --- + b + a x == 0        with sign(x) = sign
             x
*/
static double rat_root(const double a, const double b, const double c,
                       const double sign);

/* 
    Minimizes cancellation error (but not round-off ...) 
    OUTPUT: a + b + c 
*/
static double sum_3(const double a, const double b, const double c);

/*
    Chebsim: computes number of iteration and contraction factor for the    
    Chebyshev relaxation.
    OUTPUT
    - m: number of iteration
    - c: contraction factor
    INPUT:
    - rho 
    - tol: tolerance
*/
static void chebsim(double *m, double *c, const double rho, const double tol);

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
enum vv_ops {plus, minus, ewmult, ewdiv}; 
    //+, -, element-wise multiplication, element-wise division
static void vv_op(double *a, const double *b, const int n, enum vv_ops op);

/*
    Dot product between two vectors
    OUTPUT:
    - r = a (dot) b
    INPUT:
    - a: array of n doubles
    - b: array of n doubles
    - n: size of arrays
*/
static double vv_dot(const double *a, const double *b, const int n);

/*
    Array operations
    OUTPUT:
    - a: array on which to perform operation (a[i] = op(a[i]))
    INPUT:
    - a: initial array
    - n: length of the array
    - op: operation to be performed

    If op = absolute value, sqrt, multiplicative inverse, square:
        function returns 0
    If op = sum, 2-norm: 
        function returns corresponding value
*/
enum array_ops {abs_op, sqrt_op, minv_op, sqr_op, sum_op, norm2_op}; 
    // absolute value, sqrt, multiplicative inverse, square, sum, 2-norm
static double array_op(double *a, const int n, enum array_ops op);

/* 
    Array initialization
    OUTPUT:
    - a: array to be initialized (a = v*ones(n, 1))
    INPUT:
    - a: initial array
    - n: length of the array
    - v: value to initialize a
*/
static void init_array(double *a, const int n, const double v);

/*
    Operations between an array and a scalar
    OUTPUT:
    - a: array on which to perform operation (a[i] = a[i] op scal)
    INPUT:
    - a: initial array
    - scal: scalar value used for operation
    - n: length of the array
    - op: operation to be performed (addition or multiplication)
*/
enum ar_scal_ops {mult_op, add_op}; // mult, add
static void ar_scal_op(double *a, const double scal, const int n, 
    enum ar_scal_ops op);

/*******************************************************************************
* AMG Import/Export
*******************************************************************************/
/*
    Export data from the AMG setup to correct format.
    OUTPUT: 4 files
    - amg.dat
    - amg_Aff.dat
    - amg_W.dat
    - amg_AfP.dat
    INPUT:
    - structure data
*/
static void amg_export(const struct amg_setup_data *data,char *session);

/*
    Save matrices
    OUTPUT:
    - file 'filename' with matrix data
    - len(n): length of each row
    INPUT:
    - nl: number of levels
    - n: number of points
    - lvl(n): last level at which each point appears
    - id(nl, n(l)): local id for each row
    - id_l2g(n): local to global id converter
    - mat(nl): array pointing to matrices of the different levels
    - filename: name of the file
*/
static void savemats(int *len, const int n, const int nl, const int *lvl, 
    int **id, const int *id_l2g, hypre_CSRMatrix **mat, const char *filename);

/* 
    Returns max number of non zero elements on a row
*/
static int max_row_nnz(const hypre_CSRMatrix *mat);

/*
    Save vectors
    OUTPUT:
    - file 'filename' with vector data
    INPUT:
    - nl: number of levels
    - data: data with all info about setup
    - n: number of points
    - lvl(n): last level at which each point appears
    - W_len(n), Aff_len(n), AfP_len(n): length of each row of corresponding matrix
    - dvec(n): diagonal smoother for each point
    - filename: name of the file
*/
static void savevec(const int nl, const struct amg_setup_data *data, 
    const int n, const int *lvl, const int *W_len, const int *AfP_len, 
    const int *Aff_len, const double *dvec, const const char *filename);

/*
    Swap bytes for big endian / small endian consistency
*/
static double byteswap(double x);

/*
    Get size file in number of doubles
*/
static long filesize(const char *name);

/*
    Read input files
*/
static long readfile(double *data, long max, const char *name);

/*******************************************************************************
* For debugging purposes only
*******************************************************************************/
void print_csrmat(hypre_CSRMatrix *A);

#endif
