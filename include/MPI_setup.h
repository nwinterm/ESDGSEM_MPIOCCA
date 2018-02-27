#ifndef MPI_setup_H
#define MPI_setup_H

#include "mpi.h"
#include <stdio.h>










class MPI_setup
{
public:
//        MPI_setup(int argc, char *argv[]);
    MPI_setup();
    int  numtasks, rank, len, rc;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Request *reqs;   // required variable for non-blocking calls
    MPI_Request *Send_b_reqs;   // required variable for non-blocking calls
    MPI_Request *Send_q_reqs;   // required variable for non-blocking calls
    MPI_Request *Send_qX_reqs;   // required variable for non-blocking calls
    MPI_Request *Send_qY_reqs;   // required variable for non-blocking calls
    MPI_Request *Send_ViscPar_reqs;   // required variable for non-blocking calls

    MPI_Request *Recv_b_reqs;   // required variable for non-blocking calls
    MPI_Request *Recv_q_reqs;   // required variable for non-blocking calls
    MPI_Request *Recv_qX_reqs;   // required variable for non-blocking calls
    MPI_Request *Recv_qY_reqs;   // required variable for non-blocking calls
    MPI_Request *Recv_ViscPar_reqs;   // required variable for non-blocking calls

    MPI_Status *stats;   // required variable for Waitall routine

private:


};

#endif // SW1D_H

