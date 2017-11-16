
#include "MPI_Communication.h"

//#include "Constants.hpp"
  void ShareInputData(MPI_setup MPI,
                      int *N,
                      dfloat *CFL,
                      dfloat *DFL,
                      dfloat *T,
                      dfloat *g_const,
                      int *ArtificialViscosity,
                      int *PositivityPreserving,
                      dfloat *epsilon_0,
                      dfloat *sigma_min,
                      dfloat *sigma_max,
                      int *PlotVar,
                      int *NumPlots,
                      int *NumTimeChecks,
                      int *Testcase,
                      int *ES,
                      int *NumFlux,
                      int *FluxDifferencing,
                      int *rkorder,
                      int *rkSSP,
                      int *NEpad,
                      int *NEsurfpad,
                      int *Nedgepad,
                      int* KernelVersion,
                      int *KernelVersionSTD){


	 MPI_Bcast (&*N,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*CFL,1,MPI_DFLOAT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*DFL,1,MPI_DFLOAT,0,MPI_COMM_WORLD);
	 MPI_Bcast (&*T,1,MPI_DFLOAT,0,MPI_COMM_WORLD);
	 MPI_Bcast (&*g_const,1,MPI_DFLOAT,0,MPI_COMM_WORLD);
	 MPI_Bcast (&*ArtificialViscosity,1,MPI_INT,0,MPI_COMM_WORLD);
	 MPI_Bcast (&*PositivityPreserving,1,MPI_INT,0,MPI_COMM_WORLD);
	 MPI_Bcast (&*epsilon_0,1,MPI_DFLOAT,0,MPI_COMM_WORLD);
	 MPI_Bcast (&*sigma_min,1,MPI_DFLOAT,0,MPI_COMM_WORLD);
	 MPI_Bcast (&*sigma_max,1,MPI_DFLOAT,0,MPI_COMM_WORLD);
	 MPI_Bcast (&*PlotVar,1,MPI_INT,0,MPI_COMM_WORLD);
	 MPI_Bcast (&*NumPlots,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*NumTimeChecks,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*Testcase,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*ES,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*NumFlux,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*FluxDifferencing,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*rkorder,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*rkSSP,1,MPI_INT,0,MPI_COMM_WORLD);

     MPI_Bcast (&*NEpad,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*NEsurfpad,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*Nedgepad,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*KernelVersion,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast (&*KernelVersionSTD,1,MPI_INT,0,MPI_COMM_WORLD);
//int *PlotVar,int *NumPlots,int *NumTimeChecks,int *Testcase, bool *ES,int *NumFlux, bool *FluxDifferencing, bool *Cartesian,int *rkorder, bool *rkSSP
  }

  void CollectEdgeDataMPI(MPI_setup MPI, const MeshPartitioning MeshSplit, dfloat qL[], dfloat qR[]){

    int ngl = MeshSplit.ngl;
    int ngl2 = ngl*ngl;

//
//    cout << "Num Elements: " << MeshSplit.ElementsPerProc(1,1) << "\n";

  for (int i=1;i<=MeshSplit.NumProcessors;i++){


    if (MPI.rank +1 != i ) {
        int startIndex = MeshSplit.ProcIndex[i-1];
        int endIndex = MeshSplit.ProcIndex[MeshSplit.NumProcessors+i-1];
        int EdgesToSend = endIndex-startIndex+1;
        int cpuR = i-1;

        if (EdgesToSend > 0){
            int id = (startIndex) * ngl*Neq;
            int tagSend = MeshSplit.CommTags[MPI.rank*MeshSplit.NumProcessors + cpuR];//(MPI.rank+1,cpuR+1);
            int tagRecv = MeshSplit.CommTags[cpuR*MeshSplit.NumProcessors + MPI.rank];//(cpuR+1,MPI.rank+1);


            MPI_Isend(&qL[id],EdgesToSend*ngl*Neq,MPI_DFLOAT,cpuR,tagSend,MPI_COMM_WORLD,&MPI.Send_q_reqs[cpuR]);

            MPI_Irecv(&qR[id],EdgesToSend*ngl*Neq,MPI_DFLOAT,cpuR,tagRecv,MPI_COMM_WORLD,&MPI.Recv_q_reqs[cpuR]);


        }

    }

  }



  }

  void CollectEdgeDataMPI_Bonly(MPI_setup MPI, const MeshPartitioning MeshSplit, dfloat bL[], dfloat bR[]){

    int ngl = MeshSplit.ngl;
    int ngl2 = ngl*ngl;

//
//    cout << "Num Elements: " << MeshSplit.ElementsPerProc(1,1) << "\n";

  for (int i=1;i<=MeshSplit.NumProcessors;i++){


    if (MPI.rank +1 != i ) {
        int startIndex = MeshSplit.ProcIndex[i-1];
        int endIndex = MeshSplit.ProcIndex[MeshSplit.NumProcessors+i-1];
        int EdgesToSend = endIndex-startIndex+1;
        int cpuR = i-1;

        if (EdgesToSend > 0){

            int idx = (startIndex) * ngl;
            int tagSend = MeshSplit.CommTags[MPI.rank*MeshSplit.NumProcessors + cpuR];//(MPI.rank+1,cpuR+1);
            int tagRecv = MeshSplit.CommTags[cpuR*MeshSplit.NumProcessors + MPI.rank];//(cpuR+1,MPI.rank+1);


            MPI_Isend(&bL[idx],EdgesToSend*ngl,MPI_DFLOAT,cpuR,MeshSplit.NumProcessors*MeshSplit.NumProcessors+tagSend,MPI_COMM_WORLD,&MPI.Send_b_reqs[cpuR]);

            MPI_Irecv(&bR[idx],EdgesToSend*ngl,MPI_DFLOAT,cpuR,MeshSplit.NumProcessors*MeshSplit.NumProcessors+tagRecv,MPI_COMM_WORLD,&MPI.Recv_b_reqs[cpuR]);

        }

    }

  }



  }



  void CollectViscoseEdgeDataMPI(MPI_setup MPI, const MeshPartitioning MeshSplit, dfloat ViscParaL[], dfloat ViscParaR[], dfloat qGradXL[], dfloat qGradXR[], dfloat qGradYL[], dfloat qGradYR[]){

    int ngl = MeshSplit.ngl;
    int ngl2 = ngl*ngl;

//
//    cout << "Num Elements: " << MeshSplit.ElementsPerProc(1,1) << "\n";



  for (int i=1;i<=MeshSplit.NumProcessors;i++){


    if (MPI.rank +1 != i ) {
        int startIndex = MeshSplit.ProcIndex[i-1];
        int endIndex = MeshSplit.ProcIndex[MeshSplit.NumProcessors+i-1];
        int EdgesToSend = endIndex-startIndex+1;
        int cpuR = i-1;

        if (EdgesToSend > 0){
            int id = (startIndex) * ngl*Neq;
            int idx = (startIndex) * ngl;
            int tagSend = MeshSplit.CommTags[MPI.rank*MeshSplit.NumProcessors + cpuR];//(MPI.rank+1,cpuR+1);
            int tagRecv = MeshSplit.CommTags[cpuR*MeshSplit.NumProcessors + MPI.rank];//(cpuR+1,MPI.rank+1);
//            int tagSend = MeshSplit.CommTags(MPI.rank+1,cpuR+1);
//            int tagRecv = MeshSplit.CommTags(cpuR+1,MPI.rank+1);


            MPI_Isend(&qGradXL[id],EdgesToSend*ngl*Neq,MPI_DFLOAT,cpuR,  tagSend   ,MPI_COMM_WORLD,&MPI.Send_qX_reqs[cpuR]);
            MPI_Isend(&qGradYL[id],EdgesToSend*ngl*Neq,MPI_DFLOAT,cpuR,MeshSplit.NumProcessors*MeshSplit.NumProcessors  +tagSend   ,MPI_COMM_WORLD,&MPI.Send_qY_reqs[cpuR]);
            MPI_Isend(&ViscParaL[startIndex],EdgesToSend,MPI_DFLOAT,cpuR,2*MeshSplit.NumProcessors*MeshSplit.NumProcessors  +tagSend      ,MPI_COMM_WORLD,&MPI.Send_ViscPar_reqs[cpuR]);

            MPI_Irecv(&qGradXR[id],EdgesToSend*ngl*Neq,MPI_DFLOAT,cpuR,tagRecv,MPI_COMM_WORLD,&MPI.Recv_qX_reqs[cpuR]);
            MPI_Irecv(&qGradYR[id],EdgesToSend*ngl*Neq,MPI_DFLOAT,cpuR,MeshSplit.NumProcessors*MeshSplit.NumProcessors+tagRecv,MPI_COMM_WORLD,&MPI.Recv_qY_reqs[cpuR]);
            MPI_Irecv(&ViscParaR[startIndex],EdgesToSend,MPI_DFLOAT,cpuR,2*MeshSplit.NumProcessors*MeshSplit.NumProcessors  +tagRecv,MPI_COMM_WORLD,&MPI.Recv_ViscPar_reqs[cpuR]);





        }

    }

  }
  }




  void CollectSolution(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat Q[], dfloat Q_global[]){

    int ngl = MeshSplit.ngl;
    int ngl2 = ngl*ngl;

//
//    cout << "Num Elements: " << MeshSplit.ElementsPerProc(1,1) << "\n";

  for (int ie = 0; ie<MeshSplit.ElementsPerProc[0]; ie++){
//    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
    int eleID = MeshSplit.ElementLocalToGlobal[ie]-1;
    for (int i = 0; i<ngl;i++){
        for (int j=0;j<ngl;j++){
            int Gid = eleID*ngl2*Neq   +j*ngl+i;
            int Lid = ie*ngl2*Neq + j*ngl+i;
            Q_global[Gid] = Q[Lid];Gid+=ngl2;Lid+=ngl2;
            Q_global[Gid] = Q[Lid];Gid+=ngl2;Lid+=ngl2;
            Q_global[Gid] = Q[Lid];
        }

    }
  }



  for (int iproc = 1; iproc <MeshSplit.NumProcessors;iproc++){

    int varDim = ngl2*Neq * MeshSplit.ElementsPerProc[iproc];

    dfloat * q_tmp = (dfloat*) calloc(varDim,sizeof(dfloat));
//    dfloat q_tmp[varDim];

    MPI_Recv(&q_tmp[0],varDim,MPI_DFLOAT,iproc,iproc,MPI_COMM_WORLD, MPI.stats);




      for (int ie = 0; ie<MeshSplit.ElementsPerProc[iproc]; ie++){

        int eleID = MeshSplit.ElementLocalToGlobal[iproc*MeshSplit.ElementsPerProc[0]+ie]-1;
        for (int i = 0; i<ngl;i++){
            for (int j=0;j<ngl;j++){
                int Gid = eleID*ngl2*Neq   +j*ngl+i;
                int Lid = ie*ngl2*Neq + j*ngl+i;
                Q_global[Gid] = q_tmp[Lid];Gid+=ngl2;Lid+=ngl2;
                Q_global[Gid] = q_tmp[Lid];Gid+=ngl2;Lid+=ngl2;
                Q_global[Gid] = q_tmp[Lid];
            }

        }
      }

      free(q_tmp);

  }

  }



  void SendSolution(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat Q[]){

    int ngl = MeshSplit.ngl;
    int ngl2 = ngl*ngl;


    int varDim = ngl2*Neq * MeshSplit.NumElements;

    MPI_Send(&Q[0],varDim,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD);




  }











    void CollectViscPara(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat ViscPara[], dfloat ViscPara_global[]){

//
//    cout << "Num Elements: " << MeshSplit.ElementsPerProc(1,1) << "\n";

  for (int ie = 0; ie<MeshSplit.ElementsPerProc[0]; ie++){
//    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
    int eleID = MeshSplit.ElementLocalToGlobal[ie]-1;
            ViscPara_global[eleID] = ViscPara[ie];
  }



  for (int iproc = 1; iproc <MeshSplit.NumProcessors;iproc++){

    int varDim = MeshSplit.ElementsPerProc[iproc];

    dfloat * q_tmp = (dfloat*) calloc(varDim,sizeof(dfloat));
//    dfloat q_tmp[varDim];

    MPI_Recv(&q_tmp[0],varDim,MPI_DFLOAT,iproc,iproc,MPI_COMM_WORLD, MPI.stats);




      for (int ie = 0; ie<MeshSplit.ElementsPerProc[iproc]; ie++){

        int eleID = MeshSplit.ElementLocalToGlobal[iproc*MeshSplit.ElementsPerProc[0]+ie]-1;

                ViscPara_global[eleID] = q_tmp[ie];

      }

      free(q_tmp);

  }

  }



  void SendViscPara(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat ViscPara[]){




    int varDim = MeshSplit.NumElements;

    MPI_Send(&ViscPara[0],varDim,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD);



  }








    void CollectViscosity(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat Qx[],const dfloat Qy[], dfloat Qx_global[], dfloat Qy_global[]){

    int ngl = MeshSplit.ngl;
    int ngl2 = ngl*ngl;

//
//    cout << "Num Elements: " << MeshSplit.ElementsPerProc(1,1) << "\n";

  for (int ie = 0; ie<MeshSplit.ElementsPerProc[0]; ie++){
//    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
    int eleID = MeshSplit.ElementLocalToGlobal[ie]-1;
    for (int i = 0; i<ngl;i++){
        for (int j=0;j<ngl;j++){
            int Gid = eleID*ngl2*Neq   +j*ngl+i;
            int Lid = ie*ngl2*Neq + j*ngl+i;
            Qx_global[Gid] = Qx[Lid];
            Qy_global[Gid] = Qy[Lid];
            Gid+=ngl2;Lid+=ngl2;
            Qx_global[Gid] = Qx[Lid];
            Qy_global[Gid] = Qy[Lid];
            Gid+=ngl2;Lid+=ngl2;
            Qx_global[Gid] = Qx[Lid];
            Qy_global[Gid] = Qy[Lid];
        }

    }
  }



  for (int iproc = 1; iproc <MeshSplit.NumProcessors;iproc++){

    int varDim = ngl2*Neq * MeshSplit.ElementsPerProc[iproc];

    dfloat * qx_tmp = (dfloat*) calloc(varDim,sizeof(dfloat));
    dfloat * qy_tmp = (dfloat*) calloc(varDim,sizeof(dfloat));
//    dfloat q_tmp[varDim];

    MPI_Recv(&qx_tmp[0],varDim,MPI_DFLOAT,iproc,iproc,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[iproc], MPI.stats);
    MPI_Recv(&qy_tmp[0],varDim,MPI_DFLOAT,iproc,iproc,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[iproc], MPI.stats);


      for (int ie = 0; ie<MeshSplit.ElementsPerProc[iproc]; ie++){

        int eleID = MeshSplit.ElementLocalToGlobal[iproc*MeshSplit.ElementsPerProc[0]+ie]-1;
        for (int i = 0; i<ngl;i++){
            for (int j=0;j<ngl;j++){
                int Gid = eleID*ngl2*Neq   +j*ngl+i;
                int Lid = ie*ngl2*Neq + j*ngl+i;
                Qx_global[Gid] = qx_tmp[Lid];
                Qy_global[Gid] = qy_tmp[Lid];
                Gid+=ngl2;Lid+=ngl2;
                Qx_global[Gid] = qx_tmp[Lid];
                Qy_global[Gid] = qy_tmp[Lid];
                Gid+=ngl2;Lid+=ngl2;
                Qx_global[Gid] = qx_tmp[Lid];
                Qy_global[Gid] = qy_tmp[Lid];
            }

        }
      }

      free(qx_tmp);
      free(qy_tmp);

  }

  }



  void SendViscosity(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat Qx[], const dfloat Qy[]){

    int ngl = MeshSplit.ngl;
    int ngl2 = ngl*ngl;


    int varDim = ngl2*Neq * MeshSplit.NumElements;

    MPI_Send(&Qx[0],varDim,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    MPI_Send(&Qy[0],varDim,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


  }







  void GetGlobalLambdaMax(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat LocalLambdas[], dfloat * LambdaMax){



    dfloat LocalLambdaMax=0.0;
    for (int ie=0;ie<MeshSplit.NumElements;ie++){
        LocalLambdaMax = fmax(LocalLambdaMax,LocalLambdas[ie]);

    }
	//cout << "Local Lambda max: " << LocalLambdaMax <<"\n";
    MPI_Allreduce(&LocalLambdaMax, &*LambdaMax, 1, MPI_DFLOAT, MPI_MAX,MPI_COMM_WORLD);
	//cout << "Local Lambda max: " << LocalLambdaMax <<"\n";




  }



  void GetGlobalMinEleSize(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat minEleSize, dfloat * globalMinEleSize){



 dfloat localMinEleSize = minEleSize;
    MPI_Allreduce(&localMinEleSize, &*globalMinEleSize, 1, MPI_DFLOAT, MPI_MIN,MPI_COMM_WORLD);

  }


  void GetGlobalViscParaMax(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat LocalViscPara[], dfloat * ViscParaMax){



    dfloat LocalViscParaMax=0.0;
    for (int ie=0;ie<MeshSplit.NumElements;ie++){
        LocalViscParaMax = fmax(LocalViscParaMax,LocalViscPara[ie]);

    }

    MPI_Allreduce(&LocalViscParaMax, &*ViscParaMax, 1, MPI_DFLOAT, MPI_MAX,MPI_COMM_WORLD);





  }
