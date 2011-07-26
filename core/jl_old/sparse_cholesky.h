#ifndef SPARSE_CHOLESKY_H
#define SPARSE_CHOLESKY_H

/* requires "types.h" and "errmem.h" */
#if !defined(TYPES_H) || !defined(ERRMEM_H)
#warning "sparse_cholesky.h" requires "types.h" and "errmem.h"
#endif

typedef struct {
  uint n, *Lrp, *Lj;
  real *L, *D;
} sparse_cholesky_data;

/* input data is the usual CSR
   matrix is n by n
   Arp has n+1 elements
   elements of row i are A [Arp[i]], ..., A [Arp[i+1]-1]
              in columns Aj[Arp[i]], ..., Aj[Arp[i+1]-1]
*/
void sparse_cholesky_factor(uint n, const uint *Arp, const uint *Aj,
                            const real *A,
                            sparse_cholesky_data *out, buffer *buf);
                            
/* x = A^(-1) b;  works when x and b alias */
void sparse_cholesky_solve(real *x, const sparse_cholesky_data *fac, real *b);

void sparse_cholesky_free(sparse_cholesky_data *fac);

#endif

