#ifndef SARRAY_TRANSFER_H
#define SARRAY_TRANSFER_H

#if !defined(CRYSTAL_H)
#warning "sarray_transfer.h" requires "crystal.h"
#endif

#ifdef PREFIX
#  define sarray_transfer_ TOKEN_PASTE(PREFIX,sarray_transfer_)
#endif

void sarray_transfer_(array *A, crystal_data *cr, size_t off, size_t size);

#define sarray_transfer(T,A,proc_field,cr) \
  sarray_transfer_(A,cr,offsetof(T,proc_field),sizeof(T))

#endif
