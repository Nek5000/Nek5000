#ifndef CRYSTAL_H
#define CRYSTAL_H

#if !defined(COMM_H) || !defined(MEM_H)
#warning "crystal.h" requires "comm.h" and "mem.h"
#endif

#define crystal_init   PREFIXED_NAME(crystal_init  )
#define crystal_free   PREFIXED_NAME(crystal_free  )
#define crystal_router PREFIXED_NAME(crystal_router)

typedef struct {
  struct comm comm;
  buffer data, work;
  uint n;
} crystal_data;

void crystal_init(crystal_data *p, const struct comm *comm);
void crystal_free(crystal_data *p);
void crystal_router(crystal_data *p);

#endif
