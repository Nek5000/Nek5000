/* I can't believe I need a function to do this... */
#include "mpi.h"

void moab_comm_f2c_(int *fcomm, MPI_Comm *ccomm) 
{
  *ccomm = MPI_Comm_f2c(*fcomm);
}
