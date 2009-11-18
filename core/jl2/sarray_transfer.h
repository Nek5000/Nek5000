#ifndef SARRAY_TRANSFER_H
#define SARRAY_TRANSFER_H

#if !defined(CRYSTAL_H)
#warning "sarray_transfer.h" requires "crystal.h"
#endif

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
