/* I can't believe I need a function to do this... */
#ifdef MPI
#include "mpi.h"
#else
#define MPI_Comm int
#endif

void moab_comm_f2c_(int *fcomm, MPI_Comm *ccomm) 
{
#ifdef MPI  
  *ccomm = MPI_Comm_f2c(*fcomm);
#endif
}
