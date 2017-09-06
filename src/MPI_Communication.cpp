
#include "MPI_Communication.h"

//#include "Constants.hpp"


  void CollectEdgeDataMPI(MPI_setup MPI, const MeshPartitioning MeshSplit, dfloat qL[], dfloat qR[]){

    int ngl = MeshSplit.ngl;
    int ngl2 = ngl*ngl;

//
//    cout << "Num Elements: " << MeshSplit.ElementsPerProc(1,1) << "\n";

  for (int i=1;i<=MeshSplit.NumProcessors;i++){


    if (MPI.rank +1 != i ) {
        int startIndex = MeshSplit.ProcIndex(1,i);
        int endIndex = MeshSplit.ProcIndex(2,i);
        int EdgesToSend = endIndex-startIndex+1;
        int cpuR = i-1;

        if (EdgesToSend > 0){
            int id = (startIndex) * ngl*Neq;
            int tagSend = MeshSplit.CommTags(MPI.rank+1,cpuR+1);
            int tagRecv = MeshSplit.CommTags(cpuR+1,MPI.rank+1);
//            int PackSizeQ = EdgesToSend*ngl*Neq;
//            int PackSizeB = EdgesToSend*ngl;
//            dfloat * q_Recv = (dfloat*) calloc(PackSizeQ,sizeof(dfloat));
//            dfloat * q_Send = (dfloat*) calloc(PackSizeQ,sizeof(dfloat));
//            dfloat * b_Recv = (dfloat*) calloc(PackSizeB,sizeof(dfloat));
//            dfloat * b_Send = (dfloat*) calloc(PackSizeB,sizeof(dfloat));
//
//            for (int is = 1;is<=EdgesToSend;is++){
//                int qid = (is-1) * ngl*Neq;
//                int bid = (is-1) * ngl;
//                if ( MPI.rank = MeshSplit.MyEdgeInfo(8,is)){
//                    // we are left to the edge
//                } else{
//
//                }
//            }

//            cout << "MPI CASE! cpuL = " << MPI.rank << " cpuR= " << cpuR <<"\n";
//            cout <<"i am: "<< MPI.rank << " Sending " << EdgesToSend <<" Edges to " << cpuR << "with send Tag: " << tagSend << " and receive tag: "<< tagRecv <<"\n";
//            cout <<"startindex: " << startIndex << " endIndex: "<< endIndex << " id: "<< id  <<"\n";
//            cout <<"idx: "<< idx  <<"\n";
            MPI_Isend(&qL[id],EdgesToSend*ngl*Neq,MPI_DOUBLE,cpuR,tagSend,MPI_COMM_WORLD,&MPI.Send_q_reqs[cpuR]);

            MPI_Irecv(&qR[id],EdgesToSend*ngl*Neq,MPI_DOUBLE,cpuR,tagRecv,MPI_COMM_WORLD,&MPI.Recv_q_reqs[cpuR]);


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
        int startIndex = MeshSplit.ProcIndex(1,i);
        int endIndex = MeshSplit.ProcIndex(2,i);
        int EdgesToSend = endIndex-startIndex+1;
        int cpuR = i-1;

        if (EdgesToSend > 0){

            int idx = (startIndex) * ngl;
            int tagSend = MeshSplit.CommTags(MPI.rank+1,cpuR+1);
            int tagRecv = MeshSplit.CommTags(cpuR+1,MPI.rank+1);
//            int PackSizeQ = EdgesToSend*ngl*Neq;
//            int PackSizeB = EdgesToSend*ngl;
//            dfloat * q_Recv = (dfloat*) calloc(PackSizeQ,sizeof(dfloat));
//            dfloat * q_Send = (dfloat*) calloc(PackSizeQ,sizeof(dfloat));
//            dfloat * b_Recv = (dfloat*) calloc(PackSizeB,sizeof(dfloat));
//            dfloat * b_Send = (dfloat*) calloc(PackSizeB,sizeof(dfloat));
//
//            for (int is = 1;is<=EdgesToSend;is++){
//                int qid = (is-1) * ngl*Neq;
//                int bid = (is-1) * ngl;
//                if ( MPI.rank = MeshSplit.MyEdgeInfo(8,is)){
//                    // we are left to the edge
//                } else{
//
//                }
//            }

//            cout << "MPI CASE! cpuL = " << MPI.rank << " cpuR= " << cpuR <<"\n";
//            cout <<"i am: "<< MPI.rank << " Sending " << EdgesToSend <<" Edges to " << cpuR << "with send Tag: " << tagSend << " and receive tag: "<< tagRecv <<"\n";
//            cout <<"startindex: " << startIndex << " endIndex: "<< endIndex << " id: "<< id  <<"\n";
//            cout <<"idx: "<< idx  <<"\n";
            MPI_Isend(&bL[idx],EdgesToSend*ngl,MPI_DOUBLE,cpuR,MeshSplit.NumProcessors*MeshSplit.NumProcessors+tagSend,MPI_COMM_WORLD,&MPI.Send_b_reqs[cpuR]);

            MPI_Irecv(&bR[idx],EdgesToSend*ngl,MPI_DOUBLE,cpuR,MeshSplit.NumProcessors*MeshSplit.NumProcessors+tagRecv,MPI_COMM_WORLD,&MPI.Recv_b_reqs[cpuR]);

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
        int startIndex = MeshSplit.ProcIndex(1,i);
        int endIndex = MeshSplit.ProcIndex(2,i);
        int EdgesToSend = endIndex-startIndex+1;
        int cpuR = i-1;

        if (EdgesToSend > 0){
            int id = (startIndex) * ngl*Neq;
            int idx = (startIndex) * ngl;
            int tagSend = MeshSplit.CommTags(MPI.rank+1,cpuR+1);
            int tagRecv = MeshSplit.CommTags(cpuR+1,MPI.rank+1);


            MPI_Isend(&qGradXL[id],EdgesToSend*ngl*Neq,MPI_DOUBLE,cpuR,  tagSend   ,MPI_COMM_WORLD,&MPI.Send_qX_reqs[cpuR]);
            MPI_Isend(&qGradYL[id],EdgesToSend*ngl*Neq,MPI_DOUBLE,cpuR,MeshSplit.NumProcessors*MeshSplit.NumProcessors  +tagSend   ,MPI_COMM_WORLD,&MPI.Send_qY_reqs[cpuR]);
            MPI_Isend(&ViscParaL[startIndex],EdgesToSend,MPI_DOUBLE,cpuR,2*MeshSplit.NumProcessors*MeshSplit.NumProcessors  +tagSend      ,MPI_COMM_WORLD,&MPI.Send_ViscPar_reqs[cpuR]);

            MPI_Irecv(&qGradXR[id],EdgesToSend*ngl*Neq,MPI_DOUBLE,cpuR,tagRecv,MPI_COMM_WORLD,&MPI.Recv_qX_reqs[cpuR]);
            MPI_Irecv(&qGradYR[id],EdgesToSend*ngl*Neq,MPI_DOUBLE,cpuR,MeshSplit.NumProcessors*MeshSplit.NumProcessors+tagRecv,MPI_COMM_WORLD,&MPI.Recv_qY_reqs[cpuR]);
            MPI_Irecv(&ViscParaR[startIndex],EdgesToSend,MPI_DOUBLE,cpuR,2*MeshSplit.NumProcessors*MeshSplit.NumProcessors  +tagRecv,MPI_COMM_WORLD,&MPI.Recv_ViscPar_reqs[cpuR]);





        }

    }

  }
  }




  void CollectSolution(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat Q[], dfloat Q_global[]){

    int ngl = MeshSplit.ngl;
    int ngl2 = ngl*ngl;

//
//    cout << "Num Elements: " << MeshSplit.ElementsPerProc(1,1) << "\n";

  for (int ie = 0; ie<MeshSplit.ElementsPerProc(1,1); ie++){
//    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
    int eleID = MeshSplit.ElementLocalToGlobal(ie+1,1)-1;
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

    int varDim = ngl2*Neq * MeshSplit.ElementsPerProc(iproc+1,1);

    dfloat * q_tmp = (dfloat*) calloc(varDim,sizeof(dfloat));
//    dfloat q_tmp[varDim];

    MPI_Irecv(&q_tmp[0],varDim,MPI_DOUBLE,iproc,iproc,MPI_COMM_WORLD,&MPI.reqs[iproc]);
    MPI_Wait(&MPI.reqs[iproc], MPI.stats);



      for (int ie = 0; ie<MeshSplit.ElementsPerProc(iproc+1,1); ie++){

        int eleID = MeshSplit.ElementLocalToGlobal(ie+1,iproc+1)-1;
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

    MPI_Isend(&Q[0],varDim,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);




  }









  void GetGlobalLambdaMax(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat LocalLambdas[], dfloat * LambdaMax){



    dfloat * AllLambdas = (dfloat*) calloc(MeshSplit.NumProcessors,sizeof(dfloat));




    dfloat LocalLambdaMax=0.0;
    for (int ie=0;ie<MeshSplit.NumElements;ie++){
        LocalLambdaMax = max(LocalLambdaMax,LocalLambdas[ie]);

    }




    AllLambdas[MPI.rank] = LocalLambdaMax;


    if (MPI.rank == 0){
          for (int i=1;i<MeshSplit.NumProcessors;i++){
                MPI_Irecv(&AllLambdas[i],1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&MPI.reqs[i]);
            }

    }else{

        MPI_Isend(&LocalLambdaMax,1,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);

    }



    MPI_Waitall(MeshSplit.NumProcessors,MPI.reqs, MPI.stats);

    if (MPI.rank == 0){
        dfloat globalLambdaMax=0.0;
        for (int i=0;i<MeshSplit.NumProcessors;i++){
            globalLambdaMax = max(globalLambdaMax,AllLambdas[i]);

        }
        for (int i=1;i<MeshSplit.NumProcessors;i++){
            int cpuR = i;
            MPI_Isend(&globalLambdaMax,1,MPI_DOUBLE,cpuR,cpuR,MPI_COMM_WORLD,&MPI.reqs[cpuR]);
        }
        *LambdaMax=globalLambdaMax;

    }else{
        MPI_Irecv(&*LambdaMax,1,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);

    }


    MPI_Waitall(MeshSplit.NumProcessors,MPI.reqs, MPI.stats);


    free(AllLambdas);


  }



  void GetGlobalMinEleSize(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat minEleSize, dfloat * globalMinEleSize){



    dfloat * AllEleSizes;


//    cout << "RANK: " << MPI.rank << " now in GetGlobalMinEleSize \n" ;

    if (MPI.rank == 0){
        AllEleSizes = (dfloat*) calloc(MeshSplit.NumProcessors,sizeof(dfloat));
        AllEleSizes[0] = minEleSize;
          for (int i=1;i<MeshSplit.NumProcessors;i++){
                MPI_Irecv(&AllEleSizes[i],1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&MPI.reqs[i]);
            }

    }else{
        dfloat LocalMinEleSize = minEleSize;
        MPI_Isend(&LocalMinEleSize,1,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);

    }


//    cout << "RANK: " << MPI.rank << " now Waiting for first SendReceive \n" ;
    MPI_Waitall(MeshSplit.NumProcessors,MPI.reqs, MPI.stats);
//    cout << "RANK: " << MPI.rank << " now Resuming after first SendReceive \n" ;
    if (MPI.rank == 0){
        dfloat glbMinEleSize=AllEleSizes[0];
        for (int i=1;i<MeshSplit.NumProcessors;i++){
            glbMinEleSize = min(glbMinEleSize,AllEleSizes[i]);

        }
        for (int i=1;i<MeshSplit.NumProcessors;i++){
            int cpuR = i;
            MPI_Isend(&glbMinEleSize,1,MPI_DOUBLE,cpuR,cpuR,MPI_COMM_WORLD,&MPI.reqs[cpuR]);
        }
        *globalMinEleSize=glbMinEleSize;

    }else{
        MPI_Irecv(&*globalMinEleSize,1,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);

    }

//    cout << "RANK: " << MPI.rank << " now Waiting for 2nd SendReceive \n" ;
    MPI_Waitall(MeshSplit.NumProcessors,MPI.reqs, MPI.stats);
//    cout << "RANK: " << MPI.rank << " now Resuming after 2nd SendReceive \n" ;
     if (MPI.rank == 0){
        free(AllEleSizes);
     }

  }



    void CollectViscPara(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat ViscPara[], dfloat ViscPara_global[]){

//
//    cout << "Num Elements: " << MeshSplit.ElementsPerProc(1,1) << "\n";

  for (int ie = 0; ie<MeshSplit.ElementsPerProc(1,1); ie++){
//    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
    int eleID = MeshSplit.ElementLocalToGlobal(ie+1,1)-1;
            ViscPara_global[eleID] = ViscPara[ie];
  }



  for (int iproc = 1; iproc <MeshSplit.NumProcessors;iproc++){

    int varDim = MeshSplit.ElementsPerProc(iproc+1,1);

    dfloat * q_tmp = (dfloat*) calloc(varDim,sizeof(dfloat));
//    dfloat q_tmp[varDim];

    MPI_Irecv(&q_tmp[0],varDim,MPI_DOUBLE,iproc,iproc,MPI_COMM_WORLD,&MPI.reqs[iproc]);
    MPI_Wait(&MPI.reqs[iproc], MPI.stats);



      for (int ie = 0; ie<MeshSplit.ElementsPerProc(iproc+1,1); ie++){

        int eleID = MeshSplit.ElementLocalToGlobal(ie+1,iproc+1)-1;

                ViscPara_global[eleID] = q_tmp[ie];

      }

      free(q_tmp);

  }

  }



  void SendViscPara(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat ViscPara[]){




    int varDim = MeshSplit.NumElements;

    MPI_Isend(&ViscPara[0],varDim,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);




  }








    void CollectViscosity(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat Qx[],const dfloat Qy[], dfloat Qx_global[], dfloat Qy_global[]){

    int ngl = MeshSplit.ngl;
    int ngl2 = ngl*ngl;

//
//    cout << "Num Elements: " << MeshSplit.ElementsPerProc(1,1) << "\n";

  for (int ie = 0; ie<MeshSplit.ElementsPerProc(1,1); ie++){
//    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
    int eleID = MeshSplit.ElementLocalToGlobal(ie+1,1)-1;
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

    int varDim = ngl2*Neq * MeshSplit.ElementsPerProc(iproc+1,1);

    dfloat * qx_tmp = (dfloat*) calloc(varDim,sizeof(dfloat));
    dfloat * qy_tmp = (dfloat*) calloc(varDim,sizeof(dfloat));
//    dfloat q_tmp[varDim];

    MPI_Irecv(&qx_tmp[0],varDim,MPI_DOUBLE,iproc,iproc,MPI_COMM_WORLD,&MPI.reqs[iproc]);
    MPI_Wait(&MPI.reqs[iproc], MPI.stats);
    MPI_Irecv(&qy_tmp[0],varDim,MPI_DOUBLE,iproc,iproc,MPI_COMM_WORLD,&MPI.reqs[iproc]);
    MPI_Wait(&MPI.reqs[iproc], MPI.stats);


      for (int ie = 0; ie<MeshSplit.ElementsPerProc(iproc+1,1); ie++){

        int eleID = MeshSplit.ElementLocalToGlobal(ie+1,iproc+1)-1;
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

    MPI_Isend(&Qx[0],varDim,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    MPI_Isend(&Qy[0],varDim,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


  }








  void GetGlobalViscParaMax(MPI_setup MPI, const MeshPartitioning MeshSplit, const dfloat LocalViscPara[], dfloat * ViscParaMax){



    dfloat * AllViscParas;




    dfloat LocalViscParaMax=0.0;
    for (int ie=0;ie<MeshSplit.NumElements;ie++){
        LocalViscParaMax = max(LocalViscParaMax,LocalViscPara[ie]);

    }




    if (MPI.rank == 0){

         AllViscParas = (dfloat*) calloc(MeshSplit.NumProcessors,sizeof(dfloat));
         AllViscParas[0] = LocalViscParaMax;
          for (int i=1;i<MeshSplit.NumProcessors;i++){
                MPI_Irecv(&AllViscParas[i],1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&MPI.reqs[i]);
            }

    }else{

        MPI_Isend(&LocalViscParaMax,1,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);

    }



    MPI_Waitall(MeshSplit.NumProcessors,MPI.reqs, MPI.stats);

    if (MPI.rank == 0){
        dfloat globalViscParaMax=0.0;
        for (int i=0;i<MeshSplit.NumProcessors;i++){
            globalViscParaMax = max(globalViscParaMax,AllViscParas[i]);

        }
        for (int i=1;i<MeshSplit.NumProcessors;i++){
            int cpuR = i;
            MPI_Isend(&globalViscParaMax,1,MPI_DOUBLE,cpuR,cpuR,MPI_COMM_WORLD,&MPI.reqs[cpuR]);
        }
        *ViscParaMax=globalViscParaMax;

    }else{
        MPI_Irecv(&*ViscParaMax,1,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);

    }


    MPI_Waitall(MeshSplit.NumProcessors,MPI.reqs, MPI.stats);

    if (MPI.rank==0){
        free(AllViscParas);
    }



  }
