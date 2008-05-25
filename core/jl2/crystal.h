#ifndef CRYSTAL_H
#define CRYSTAL_H

#if !defined(ERRMEM_H) || !defined(TYPES_H) || !defined(COMM_H)
#warning "crystal.h" requires "errmem.h", "types.h", and "comm.h"
#endif

#ifdef PREFIX
#  define crystal_init   TOKEN_PASTE(PREFIX,crystal_init)
#  define crystal_free   TOKEN_PASTE(PREFIX,crystal_free)
#  define crystal_router TOKEN_PASTE(PREFIX,crystal_router)
#endif

typedef struct {
  comm_t comm;
  buffer data, work;
  uint n;
} crystal_data;

void crystal_init(crystal_data *p, const comm_t *comm);
void crystal_free(crystal_data *p);
void crystal_router(crystal_data *p);

#endif
