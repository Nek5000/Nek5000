#ifndef SPARSE_CHOLESKY_H
#define SPARSE_CHOLESKY_H

#if !defined(TYPES_H) || !defined(ERRMEM_H)
#warning "sparse_cholesky.h" requires "types.h" and "errmem.h"
#endif

#ifdef PREFIX
#  define sparse_cholesky_factor TOKEN_PASTE(PREFIX,sparse_cholesky_factor)
#  define sparse_cholesky_solve  TOKEN_PASTE(PREFIX,sparse_cholesky_solve )
#  define sparse_cholesky_free   TOKEN_PASTE(PREFIX,sparse_cholesky_free  )
#endif

typedef struct {
  uint n, *Lrp, *Lj;
  double *L, *D;
} sparse_cholesky_data;

/* input data is the usual CSR
   matrix is n by n
   Arp has n+1 elements
   elements of row i are A [Arp[i]], ..., A [Arp[i+1]-1]
              in columns Aj[Arp[i]], ..., Aj[Arp[i+1]-1]
*/
void sparse_cholesky_factor(uint n, const uint *Arp, const uint *Aj,
                            const double *A,
                            sparse_cholesky_data *out, buffer *buf);
                            
/* x = A^(-1) b;  works when x and b alias */
void sparse_cholesky_solve(
  double *x, const sparse_cholesky_data *fac, double *b);

void sparse_cholesky_free(sparse_cholesky_data *fac);

#endif

