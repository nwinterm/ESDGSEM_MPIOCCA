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


#include "MPI_setup.h"



//MPI_setup::MPI_setup(int argc, char *argv[]){

MPI_setup::MPI_setup()
{


//    MPI_Init(&argc,&argv);
    MPI_Init(NULL,NULL);

    // get number of tasks
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    // get my rank
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // this one is obvious
    MPI_Get_processor_name(hostname, &len);
    printf ("Number of tasks= %d My rank= %d Running on %s\n", numtasks,rank,hostname);


    reqs = new MPI_Request[numtasks];
    Send_b_reqs = new MPI_Request[numtasks];
    Send_q_reqs= new MPI_Request[numtasks];
    Send_qX_reqs= new MPI_Request[numtasks];
    Send_qY_reqs= new MPI_Request[numtasks];
    Send_ViscPar_reqs= new MPI_Request[numtasks];
    Recv_b_reqs= new MPI_Request[numtasks];
    Recv_q_reqs= new MPI_Request[numtasks];
    Recv_qX_reqs= new MPI_Request[numtasks];
    Recv_qY_reqs= new MPI_Request[numtasks];
    Recv_ViscPar_reqs= new MPI_Request[numtasks];



    for (int i=0; i<numtasks; i++)
    {

        reqs[i]   = MPI_REQUEST_NULL;
        Send_b_reqs[i] = MPI_REQUEST_NULL;
        Send_q_reqs[i]= MPI_REQUEST_NULL;
        Send_qX_reqs[i]= MPI_REQUEST_NULL;
        Send_qY_reqs[i]= MPI_REQUEST_NULL;
        Send_ViscPar_reqs[i]= MPI_REQUEST_NULL;
        Recv_b_reqs[i]= MPI_REQUEST_NULL;
        Recv_q_reqs[i]= MPI_REQUEST_NULL;
        Recv_qX_reqs[i]= MPI_REQUEST_NULL;
        Recv_qY_reqs[i]= MPI_REQUEST_NULL;
        Recv_ViscPar_reqs[i]= MPI_REQUEST_NULL;
    }

    stats = new MPI_Status[numtasks];   // required variable for Waitall routine
}


