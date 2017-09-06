#include "MPI_setup.h"



MPI_setup::MPI_setup(int argc, char *argv[]){

    MPI_Init(&argc,&argv);

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



    for (int i=0; i<numtasks;i++){

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


