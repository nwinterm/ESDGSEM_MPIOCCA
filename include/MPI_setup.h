//==================================================================================================================================
// Copyright (c) 2019 Niklas Wintermeyer
// Copyright (c) 2019 Gregor Gassner
// Copyright (c) 2019 Andrew Winters
//
// This file is part of ESDGSEM_MPIOCCA (github.com/ESDGSEM_MPIOCCA). ESDGSEM_MPIOCCA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version.
//
// ESDGSEM_MPIOCCA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
//
// You should have received a copy of the GNU General Public License along with ESDGSEM_MPIOCCA. If not, see <http://www.gnu.org/licenses/>.


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

