#include <stdio.h>
#include <mpi.h>

#define NTIMERS 2

int SYNC=1;
int TIMING=1;
double TIMERS[NTIMERS]={0.0};


#pragma weak NEK_COMM_SETTINGS   = nek_comm_settings
#pragma weak nek_comm_settings_  = nek_comm_settings
#pragma weak nek_comm_settings__ = nek_comm_settings
void nek_comm_settings(int *sync,int *timing)
{
     SYNC   = *sync;
     TIMING = *timing; 
}

#pragma weak NEK_COMM_GETTIMERS   = nek_comm_gettimers
#pragma weak nek_comm_gettimers_  = nek_comm_gettimers
#pragma weak nek_comm_gettimers__ = nek_comm_gettimers
void nek_comm_gettimers(double *timers)
{
     int i;
     for(i=0; i<NTIMERS; i++) timers[i]=TIMERS[i];
}

#pragma weak NEK_COMM_RESETTIMERS   = nek_comm_resettimers
#pragma weak nek_comm_resettimers_  = nek_comm_resettimers
#pragma weak nek_comm_resettimers__ = nek_comm_resettimers
void nek_comm_resettimers(void)
{
     int i;
     for(i=0; i<NTIMERS; i++) TIMERS[i]=0.0;
}

#pragma weak MPI_ALLREDUCE   = mpi_allreduce_f
#pragma weak mpi_allreduce   = mpi_allreduce_f
#pragma weak mpi_allreduce_  = mpi_allreduce_f
#pragma weak mpi_allreduce__ = mpi_allreduce_f
void mpi_allreduce_f(char *sendbuf, char *recvbuf, MPI_Fint *count,
                    MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *comm,
                    MPI_Fint *ierr)
{
    MPI_Comm c_comm;
    MPI_Datatype c_type;
    MPI_Op c_op;
    double t0,t1;

    c_comm = MPI_Comm_f2c(*comm);
    c_type = MPI_Type_f2c(*datatype);
    c_op   = MPI_Op_f2c(*op);

    if(TIMING) t0 = MPI_Wtime(); 
    if(SYNC) *ierr = MPI_Barrier(c_comm);
    if(TIMING) t1 = MPI_Wtime(); 
    *ierr = MPI_Allreduce(sendbuf, recvbuf, *count, c_type, c_op, c_comm);
    if(TIMING) {TIMERS[0] += MPI_Wtime()-t1; TIMERS[1] += t1-t0;}

}

#pragma weak MPI_Allreduce  = mpi_allreduce_c
#pragma weak MPI_Allreduce_ = mpi_allreduce_c
int mpi_allreduce_c(void *sendbuf, void *recvbuf, int count,
                    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    int ierr;
    double t0,t1;

    if(TIMING) t0 = MPI_Wtime(); 
    if(SYNC) ierr = MPI_Barrier(comm);
    if(TIMING) t1 = MPI_Wtime(); 
    ierr = MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    if(TIMING) {TIMERS[0] += MPI_Wtime()-t1; TIMERS[1] += t1-t0;}
    return ierr;
}
