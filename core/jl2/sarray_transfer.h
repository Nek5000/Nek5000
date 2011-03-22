#ifndef SARRAY_TRANSFER_H
#define SARRAY_TRANSFER_H

#if !defined(CRYSTAL_H)
#warning "sarray_transfer.h" requires "crystal.h"
#endif

/*
  High-level interface for the crystal router.
  Given an array of structs, transfers each to the process indicated
  by a field of the struct, which gets set to the source process on output.
  
  For the dynamic "array" type, see "mem.h".
  
  Requires a "crystal router" object:
  
    struct comm c;
    crystal_data cr;
    
    comm_init(&c, MPI_COMM_WORLD);
    crystal_init(&cr, &c);
    
  Example sarray_transfer usage:
  
    struct T { ...; uint proc; ...; };
    array A = null_array;
    struct T *p, *e;
    
    // resize A to 100 struct T's, fill up with data
    p = array_reserve(struct T, &A, 100), A.n=100;
    for(e=p+A.n;p!=e;++p) {
      ...
      p->proc = ...;
      ...
    }
    
    // array A represents the array
    //   struct T ar[A.n]    where &ar[0] == A.ptr
    // transfer ar[i] to processor ar[i].proc  for each i=0,...,A.n-1:
    
    sarray_transfer(struct T, A, proc, &cr);
    
    // now array A represents a different array with a different size
    //   struct T ar[A.n]    where &ar[0] == A.ptr
    // where now ar[i] came from processor ar[i].proc
    // the ordering is arbitrary
    
    sarray_transfer(struct T, A, proc, &cr);
    // note: this second call should return A to its original state,
    //       up to ordering
 
  Cleanup:
    array_free(&A);
    crystal_free(&cr);
    comm_free(&c);

  Example sarray_transfer_ext usage:
  
    struct T { ... };
    array A;
    uint proc[A.n];
    
    // array A represents the array
    //   struct T ar[A.n]    where &ar[0] == A.ptr
    // transfer ar[i] to processor proc[i]  for each i=0,...,A.n-1:
    sarray_transfer_ext(struct T, &A, proc, &cr);
    
    // no information is available now on where each struct came from

*/

#define sarray_transfer_     PREFIXED_NAME(sarray_transfer_    )
#define sarray_transfer_ext_ PREFIXED_NAME(sarray_transfer_ext_)

void sarray_transfer_(array *A, size_t size, size_t off,
                      int fixed, crystal_data *cr);
void sarray_transfer_ext_(array *A, size_t size, const uint *proc,
                          crystal_data *cr);

#define sarray_transfer(T,A,proc_field,cr) \
  sarray_transfer_(A,sizeof(T),offsetof(T,proc_field),0,cr)

#define sarray_transfer_ext(T,A,proc,cr) \
  sarray_transfer_ext_(A,sizeof(T),proc,cr)

#endif
