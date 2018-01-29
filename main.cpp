#include <iostream>
#include <fstream>
#include "occa.hpp"
#include "Constants.hpp"
#include "MPI_setup.h"
#include "basis.h"
#include "SW2D.h"
#include "RungeKutta.h"
#include "Mesh.h"
#include "ReadInput.h"
#include "plots.h"
#include "MeshPartitioning.h"
#include "MPI_Communication.h"
#include <ctime>
#include <sys/time.h>


using namespace std;
//using namespace Constants;


int main(int argc, char *argv[])
{

    //   int  numtasks, rank, len, rc;
    //   char hostname[MPI_MAX_PROCESSOR_NAME];


    //MPI_setup MPI(argc, argv);
    MPI_setup MPI;

    if(argc<2)
    {
        printf("usage:   ./main [model] \n");
        printf("example: ./main  CUDA\n");
        printf("example: ./main  OpenCL\n");
        printf("example: ./main  OpenMP\n");
        exit(-1);
    }


    if (MPI.rank == 0)
    {
        occa::printAvailableDevices();
    }

    if(MPI.rank==0)
    {
        cout << "Initializing device...";
    }

    occa::device device;
    if(!strcmp(argv[1], "Serial"))
        device.setup("mode = Serial");
    else if(!strcmp(argv[1], "OpenMP"))
        device.setup("mode = OpenMP  , schedule = compact, chunk = 10");
    else if(!strcmp(argv[1], "OpenCL"))
        if(MPI.rank==0)  //HACK WAY FOR 2 GPUs
        {
            device.setup("mode = OpenCL  , platformID = 0, deviceID = 0");
        }
        else
        {
            device.setup("mode = OpenCL  , platformID = 0, deviceID = 1");
        }
    else if(!strcmp(argv[1], "CUDA"))
        if(MPI.rank==0)     //HACK WAY FOR 2 GPUs
        {
            device.setup("mode = CUDA    , deviceID = 0");
        }
        else
        {
            device.setup("mode = CUDA    , deviceID = 1");
        }

    occa::kernelInfo info;

    if(!strcmp(argv[1], "OpenMP"))
    {
        info.addCompilerFlag("-O3");
        //    info.addCompilerFlag("-ftree-vectorizer-verbose=7");
    }

    if(!strcmp(argv[1], "CUDA"))
    {
        info.addCompilerFlag("--ftz=true");
        info.addCompilerFlag("--prec-div=false");
        info.addCompilerFlag("--prec-sqrt=false");
        info.addCompilerFlag("--use_fast_math");
        info.addCompilerFlag("--fmad=true"); // compiler option for cuda
    }

    if(!strcmp(argv[1], "OpenCL"))
    {
        info.addCompilerFlag("-cl-strict-aliasing");
        info.addCompilerFlag("-cl-mad-enable");
        info.addCompilerFlag("-cl-no-signed-zeros");
        info.addCompilerFlag("-cl-unsafe-math-optimizations");
        info.addCompilerFlag("-cl-finite-math-only");
        info.addCompilerFlag("-cl-fast-relaxed-math");
    }

    if(MPI.rank==0)
    {
        cout <<" ... device initialized.\n";
        cout <<" reading in possible other input parameters ... \n.";
    }








	cout.precision(17);
    //initialize polynomial basis

    int ngl, ngl2;;
    int NoDofs,NoSpaceDofs;
    int Nfaces, Nelem;
    int Nfaces_global, Nelem_global,NoDofs_global,NoSpaceDofs_global;
    int PlotVar;
    int EntropyPlot;
	int ArtificialViscosity;
    int PositivityPreserving;
    int rkorder;
    int rkSSP;
    int ES,Cartesian,Fluxdifferencing;
    int NumFlux;
    string meshFile;
    dfloat dt,T,CFL,DFL;
    dfloat g_const;
    dfloat epsilon_0,sigma_min,sigma_max;
    dfloat t=0.0;
    dfloat PosPresTOL=0.0;
    int NumPlots,NumTimeChecks,Testcase;
    int Nedgepad;
    int NEsurfpad;
    int KernelVersion=-1;
    int KernelVersionSTD=-1;
	

    int N=0;
    int NelemX=0;
    int NEpad=0;

    for (int i = 2; i < argc; i+=2)
    {
        /* We will iterate over argv[] to get the parameters stored inside.
        							  * Note that we're starting on 1 because we don't need to know the
        							  * path of the program, which is stored in argv[0] */
        if (i + 1 != argc) // Check that we haven't finished parsing already
            if(!strcmp(argv[i], "Nelem"))
            {
                NelemX = atoi(argv[i + 1]);
            }
            else if(!strcmp(argv[i], "N"))
            {
                N =atoi(argv[i + 1]);
            }
            else if(!strcmp(argv[i], "Nepad"))
            {
                NEpad =atoi(argv[i + 1]);
            }
            else if(!strcmp(argv[i], "STDKernelV"))
            {
                KernelVersionSTD =atoi(argv[i + 1]);
            }
            else if(!strcmp(argv[i], "KernelV"))
            {
                KernelVersion =atoi(argv[i + 1]);
            }
            else
            {
                std::cout << "Not enough or invalid arguments, please try again.\n";
                std::cout <<"invalid argument: " << argv[i]<<"\n";
                MPI_Finalize();
                exit(0);
            }
    }

    if(MPI.rank==0)
    {
        cout << "rank 0 reading input file \n";
        ReadInputFile(&N,
                      &meshFile,
                      &CFL,
                      &DFL,
                      &T,
                      &g_const,
                      &ArtificialViscosity,
                      &PositivityPreserving,
                      &PosPresTOL,
		      &epsilon_0,
                      &sigma_min,
                      &sigma_max,
                      &PlotVar,
					  &EntropyPlot,
                      &NumPlots,
                      &NumTimeChecks,
                      &Testcase,
                      &ES,
                      &NumFlux,
                      &Fluxdifferencing,
                      &Cartesian,
                      &rkorder,
                      &rkSSP,
                      &NEpad,
                      &NEsurfpad,
                      &Nedgepad,
                      &KernelVersion,
                      &KernelVersionSTD);
    }

    ShareInputData(MPI,
                   &N,
                   &CFL,
                   &DFL,
                   &T,
                   &g_const,
                   &ArtificialViscosity,
                   &PositivityPreserving,
                   &PosPresTOL,
		   &epsilon_0,
                   &sigma_min,
                   &sigma_max,
                   &PlotVar,
				   &EntropyPlot,
                   &NumPlots,
                   &NumTimeChecks,
                   &Testcase,
                   &ES,
                   &NumFlux,
                   &Fluxdifferencing,
                   &rkorder,
                   &rkSSP,
                   &NEpad,
                   &NEsurfpad,
                   &Nedgepad,
                   &KernelVersion,
                   &KernelVersionSTD);

    if (Testcase == 31)
    {
        dfloat h0 = 0.1;
        dfloat a=1.0;
        dfloat omega = sqrt(2*g_const*h0)/a;
        T = 2*PI / omega;
        T=2*T;
    }
    ngl=N+1;
    ngl2=ngl*ngl;


    if(MPI.rank==0)
    {
        cout <<"Polynomial Order: " <<N <<"\n";
        cout <<"Testcase: " <<Testcase <<"\n";
        cout <<"CFL: " <<CFL <<"\n";
        cout <<"T: " <<T <<"\n";
        cout <<"g_const: " <<g_const <<"\n";
        cout <<"ArtificialViscosity: " <<ArtificialViscosity <<"\n";
        cout <<"PositivityPreserving: " <<PositivityPreserving <<"\n";
        cout <<"PosPresTOL: " <<PosPresTOL <<"\n";
	cout <<"epsilon_0: " <<epsilon_0 <<"\n";
        cout <<"sigma_min: " <<sigma_min <<"\n";
        cout <<"sigma_max: " <<sigma_max <<"\n";
        cout <<"PlotVar: " <<PlotVar <<"\n";
        cout <<"NumPlots: " <<NumPlots <<"\n";
        cout <<"Fluxdifferencing: " <<Fluxdifferencing <<"\n";
        cout <<"NumFlux: " <<NumFlux <<"\n";
        cout <<"Cartesian: " <<Cartesian <<"\n";
        cout <<"rkorder: " <<rkorder <<"\n";
        cout <<"rkSSP: " <<rkSSP <<"\n";
        cout <<"Setting up Basis for N=" << N <<".\n";
    }


    basis DGBasis(N,Fluxdifferencing);



    if(MPI.rank==0)
    {
        if (ES)
        {
            cout <<"Entropy Stability activated.\n";
        }
        else
        {
            cout <<"Entropy Stability deactivated.\n";
        }
    }


    cout <<"Creating Mesh Partitioning .\n";
    MeshPartitioning DGMeshPartition(MPI.numtasks, ngl);

    Mesh DGMesh(DGBasis.x_GL,ngl,NelemX);

    if(MPI.rank==0)
    {
        if (Cartesian)
        {
            cout <<"Cartesian Mesh will be generated.\n";
        }
        else
        {
            cout <<"Reading Global Mesh from Inputfile.\n";
        }
        DGMesh.InitMesh(meshFile,Cartesian,Testcase);
        Nelem_global=DGMesh.m_num_elements;
        Nfaces_global=DGMesh.m_num_edges;
        NoDofs_global=ngl2*Nelem_global*Neq;
        NoSpaceDofs_global=ngl2*Nelem_global;





        //cout << " We solve " << DGMesh.m_num_edges << " Faces with "<< ngl << " nodes each, so we should have " << DGMesh.m_num_edges*ngl << " IDs! \n" ;
        cout <<"Global Mesh created.\n";
        DGMeshPartition.DivideMesh(DGMesh,MPI);
        cout <<"Mesh Partitions created! \n";
    }
    else
    {

        DGMeshPartition.ReceiveMesh(MPI);

    }



    //if(MPI.rank==3){cout << "I did get my Mesh Partition! \n...";}
    //cout << " We solve " << DGMeshPartition.NumEdges << " Faces with "<< ngl << " nodes each, so we should have " << DGMeshPartition.NumEdges*ngl << " IDs! \n" ;


    Nelem=DGMeshPartition.NumElements;
    Nfaces=DGMeshPartition.NumEdges;
    DGBasis.setNelemLocal(Nelem);
    if(MPI.rank==0)
    {
        DGBasis.setNelemGlobal(Nelem_global);
    }

    NoDofs=ngl2*Nelem*Neq;
    NoSpaceDofs=ngl2*Nelem;



    cout <<"Sorting Edges for MPI! \n";
    DGMeshPartition.SortMPIEdges(MPI);
    cout <<"done! \n";

    dfloat LambdaMax;

    dfloat * ViscPara = (dfloat*) calloc(Nelem,sizeof(dfloat));

    dfloat * ViscPara_Global;


    dfloat * Q_global ;
	dfloat * EntropyOverTime;
	dfloat * MassOverTime;
	dfloat * EntropyTimes;
    dfloat * QtVisc_global;
    dfloat * Qx_global;
    dfloat * Qy_global ;

    dfloat * J_global;
    dfloat * b_global;

    dfloat * x_phy_global;
    dfloat * y_phy_global;

    if (MPI.rank ==0 )
    {
        ViscPara_Global = (dfloat*) calloc(Nelem_global,sizeof(dfloat));


        Q_global = (dfloat*) calloc(NoDofs_global,sizeof(dfloat));
        if (ArtificialViscosity==1)
        {
            Qx_global = (dfloat*) calloc(NoDofs_global,sizeof(dfloat));
            Qy_global = (dfloat*) calloc(NoDofs_global,sizeof(dfloat));
            QtVisc_global = (dfloat*) calloc(NoDofs_global,sizeof(dfloat));
        }

        J_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));
        b_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));

        x_phy_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));
        y_phy_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));
		
		if (EntropyPlot){
			EntropyOverTime = (dfloat*) calloc(NumPlots,sizeof(dfloat));
			MassOverTime = (dfloat*) calloc(NumPlots,sizeof(dfloat));
			EntropyTimes = (dfloat*) calloc(NumPlots,sizeof(dfloat));
		}




        for(int ie=0; ie<Nelem_global; ++ie)
        {

            for(int j=0; j<ngl; ++j)
            {
                for(int i=0; i<ngl; ++i)
                {
                    //erst alle Punkte die zu j=0 gehoeren, dann j=1, usw.
                    int id = ie*ngl2   +j*ngl+i;
                    int Qid =ie*ngl2*Neq + j*ngl+i;
                    //        int id = (ieY*NelemX+ieX)*ngl2   +j*ngl+i+1;
                    J_global[id] = 1.0/DGMesh.J_global[id];
                    x_phy_global[id] = DGMesh.x_global[id];
                    y_phy_global[id] = DGMesh.y_global[id];


                }
            }
        }
        InitB(0,DGMeshPartition,Testcase,Nelem_global,ngl,ngl2,x_phy_global,y_phy_global,b_global);



    }






    if(MPI.rank==0)
    {
        cout <<"Initializing Host Memory... ";
    }
    // THESE ARE TEMPORARY VARIABLES AND FREED BEFORE THE TIME LOOP
    dfloat * Qt = (dfloat*) calloc(NoDofs,sizeof(dfloat));
    dfloat * gRK = (dfloat*) calloc(NoDofs,sizeof(dfloat));
    dfloat * ElementSizes= (dfloat*) calloc(Nelem,sizeof(dfloat));
    dfloat minEleSize, globalMinEleSize;
    dfloat * x_xi = (dfloat*) calloc(NoSpaceDofs,sizeof(dfloat));
    dfloat * y_xi = (dfloat*) calloc(NoSpaceDofs,sizeof(dfloat));
    dfloat * x_eta = (dfloat*) calloc(NoSpaceDofs,sizeof(dfloat));
    dfloat * y_eta = (dfloat*) calloc(NoSpaceDofs,sizeof(dfloat));
    dfloat * Bx = (dfloat*) calloc(NoSpaceDofs,sizeof(dfloat));
    dfloat * By = (dfloat*) calloc(NoSpaceDofs,sizeof(dfloat));
    int * EdgeData = (int*) calloc(8*(Nfaces),sizeof(int));
    int * ElemToEdge = (int*) calloc(4*Nelem,sizeof(int));
    int * ElemEdgeMasterSlave = (int*) calloc(4*Nelem,sizeof(int));
    int * ElemEdgeOrientation = (int*) calloc(4*Nelem,sizeof(int));
    dfloat * Scal = (dfloat*) calloc(ngl*(Nfaces),sizeof(dfloat));
    dfloat * normalsX = (dfloat*) calloc(ngl*(Nfaces),sizeof(dfloat));
    dfloat * normalsY = (dfloat*) calloc(ngl*(Nfaces),sizeof(dfloat));
    dfloat * x_phy = (dfloat*) calloc(NoSpaceDofs,sizeof(dfloat));
    dfloat * y_phy = (dfloat*) calloc(NoSpaceDofs,sizeof(dfloat));
    dfloat * b = (dfloat*) calloc(NoSpaceDofs,sizeof(dfloat));
    dfloat * J = (dfloat*) calloc(NoSpaceDofs,sizeof(dfloat));




    // THESE ARE NEEDED THROUGHOUT RUNTIME FOR MPI COMMUNICATION
    dfloat * q = (dfloat*) calloc(NoDofs,sizeof(dfloat));
    dfloat * qL = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));
    dfloat * qR = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));
    dfloat * bL = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));
    dfloat * bR = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));
    dfloat * ViscParaL = (dfloat*) calloc(Nfaces,sizeof(dfloat));
    dfloat * ViscParaR = (dfloat*) calloc(Nfaces,sizeof(dfloat));
    dfloat * qGradientXL = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));
    dfloat * qGradientXR = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));
    dfloat * qGradientYL = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));
    dfloat * qGradientYR = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));

    dfloat * Qx = (dfloat*) calloc(NoDofs,sizeof(dfloat));
    dfloat * Qy = (dfloat*) calloc(NoDofs,sizeof(dfloat));

    dfloat * QtVisc = (dfloat*) calloc(NoDofs,sizeof(dfloat));


    dfloat * SurfaceParts = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));
    if(MPI.rank==0)
    {
        cout <<"... finished.\n";
    }




    if(MPI.rank==0)
    {
        cout <<"Loading Mesh Data into global arrays on HOST... \n";
    }
    for(int is=0; is<Nfaces; ++is)
    {
        for (int info=0; info<8; info++)
        {
            int id = is*8 + info;
            int id2 = is*10 + info;
            EdgeData[id] = DGMeshPartition.MyEdgeInfo[id2+2];  //left element

        }

    }
    for(int ie=1; ie<=Nelem; ++ie)
    {
        int id = (ie-1)*4;
        int idFace;
        for (int is=0; is<4; is++)
        {
            int ifa  = DGMeshPartition.MyElementToEdge[id+is];//(is+1,ie);
            idFace = (ifa-1)*8;
            if ((EdgeData[idFace+5] == MPI.rank) && (EdgeData[idFace]==ie-1))
            {
                //this is the left element to this edge!
                ElemEdgeMasterSlave[id+is]=+1;
                ElemEdgeOrientation[id+is]=1;   //order is never reversed for left element! as it is master
            }
            else
            {
                ElemEdgeMasterSlave[id+is]=-1;
                ElemEdgeOrientation[id+is]=EdgeData[idFace+4];
            }
            ElemToEdge[id+is]   = ifa;
        }
    }
    for(int ie=0; ie<Nelem; ++ie)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2   +j*ngl+i;
                x_phy[id] = DGMeshPartition.x_global[id];
                y_phy[id] = DGMeshPartition.y_global[id];
                x_xi[id] = DGMeshPartition.xXi_global[id];
                y_xi[id] = DGMeshPartition.yXi_global[id];
                x_eta[id] = DGMeshPartition.xEta_global[id];
                y_eta[id] = DGMeshPartition.yEta_global[id];
                J[id] = 1.0/DGMeshPartition.J_global[id];
            }
        }
    }
    for(int is=0; is<Nfaces; ++is)
    {
        for (int i=0; i<ngl; i++)
        {
            int id = is*ngl+i;
            normalsX[id]=DGMeshPartition.nx_global[id];
            normalsY[id]=DGMeshPartition.ny_global[id];
            Scal[id] = DGMeshPartition.scal_global[id];
        }
    }
//        for(int ie=0;ie<Nfaces;++ie){
//        cout <<"\n face: " << ie <<"\n";
//            for(int i=0;i<ngl;++i){
//                int id = ie*ngl  +i ;
//                cout << normalsX[id] << " ";
//          }
//           cout  <<"\n";
//        }


    if(MPI.rank==0)
    {
        cout <<"          ...finished.\n";
    }


    DGBasis.calcElementSizes(J,ElementSizes,&minEleSize);
    if(MPI.rank==0)
    {
        cout <<" Element Sizes calculated... \n";
    }
    //    cout << "RANK: " << MPI.rank << " now entering GetGlobalMinEleSize \n" ;
    GetGlobalMinEleSize(MPI,  DGMeshPartition,minEleSize, &globalMinEleSize);
    if(MPI.rank==0)
    {
        cout <<" ... and distributed \n";
    }
    //    cout << "RANK: " << MPI.rank << " local min ele Size: "<< minEleSize << " global min Ele Size: " << globalMinEleSize<< "\n";




    dfloat Dmat[ngl2];
    dfloat Dmat0[ngl2];
    dfloat Dhat[ngl2];
    dfloat VdmInv[ngl2];
	dfloat DCentralFD[ngl2];
	dfloat DupwindFD[ngl2];
	dfloat DdownwindFD[ngl2];
    //dfloat SubCellMat[ngl2];
    dfloat GLw[ngl];
    for (int i=0; i<ngl; i++)
    {
        for (int l=0; l<ngl; l++)
        {
            int Did =   i*ngl + l;
            Dmat[Did] =DGBasis.D[Did];
            Dmat0[Did] =DGBasis.D0[Did];;
            Dhat[Did] = DGBasis.Dhat[Did];
            VdmInv[Did] = DGBasis.VdmInv[Did];
			DCentralFD[Did] = DGBasis.DCentralFD[Did];
			DupwindFD[Did] = DGBasis.DupwindFD[Did];
			DdownwindFD[Did] = DGBasis.DdownwindFD[Did];
            //            SubCellMat[Did] = DGBasis.SubCellMat(i+1,l+1);
        }
        GLw[i] = DGBasis.w_GL[i];
    }

           cout <<"\n D central: \n";
        for(int j=0;j<ngl;++j){
            for(int i=0;i<ngl;++i){
                int id =   j*ngl+i;
                cout <<DCentralFD[id]<<"  ";
          }
            cout <<"\n";
        }

	           cout <<"\n D upwind: \n";
        for(int j=0;j<ngl;++j){
            for(int i=0;i<ngl;++i){
                int id =   j*ngl+i;
                cout <<DupwindFD[id]<<"  ";
          }
            cout <<"\n";
        }

	           cout <<"\n D downwind: \n";
        for(int j=0;j<ngl;++j){
            for(int i=0;i<ngl;++i){
                int id =   j*ngl+i;
                cout <<DdownwindFD[id]<<"  ";
          }
            cout <<"\n";
        }
    //initialise time integrator
    RungeKutta RK(rkorder,rkSSP);

    //INITIALIZE SOLUTION
    InitB(1,DGMeshPartition,Testcase,Nelem,ngl,ngl2,x_phy,y_phy,b);

    CalcBDerivatives(Nelem,ngl,ngl2,g_const,x_phy,y_phy,b,Dmat0,y_eta,y_xi,x_eta,x_xi,Bx,By,J);
//	cout << "i am rank: " << MPI.rank << " and my Nelem_global is " << DGMeshPartition.global_NumElements << "\n";

    InitQ(1,DGMeshPartition,Testcase,Nelem,ngl,ngl2,x_phy,y_phy,q,0.0,b, g_const);


//    dfloat * q_modal = (dfloat*) malloc(NoDofs*sizeof(dfloat));
//    dfloat * q_exakt = (dfloat*) calloc(NoDofs_global,sizeof(dfloat));
//    DGBasis.ConvertToModal(q, q_modal);
//           cout <<"\n Modal Coefficients: \n";
//    for (int ie=0;ie<Nelem;ie++){
//            cout <<"Ele: " << ie <<"\n";
//        for(int j=0;j<ngl;++j){
//            for(int i=0;i<ngl;++i){
//                int id = ie*ngl2*Neq +  j*ngl+i;
//                cout <<q_modal[id]<<"  ";
//          }
//            cout <<"\n";
//        }
//    }

//    DGBasis.EvaluteModalPolynomial(q_modal, q_exakt);
//           cout <<"\n q_init : \n";
//    for (int ie=0;ie<1;ie++){
//            cout <<"Ele: " << ie <<"\n";
//        for(int j=0;j<ngl;++j){
//            for(int i=0;i<ngl;++i){
//                int id = ie*ngl2*Neq +  j*ngl+i;
//               cout <<q[id]<<"  ";
//          }
//            cout <<"\n";
//        }
//    } 
//           cout <<"\n q_exakt : \n";
//    for (int ie=0;ie<1;ie++){
//            cout <<"Ele: " << ie <<"\n";
//        for(int j=0;j<ngl;++j){
//            for(int i=0;i<ngl;++i){
//                int id = ie*ngl2*Neq +  j*ngl+i;
//               cout <<q_exakt[id]<<"  ";
//          }
//            cout <<"\n";
//        }
//    }


    //
    //    dfloat L2modalVnodal[Neq];
    //    dfloat LinfmodalVnodal[Neq];
    //    for(int ik=0;ik<Neq;ik++){
    //        L2modalVnodal[ik]=0.0;
    //        LinfmodalVnodal[ik]=0.0;
    //    }
    //    DGBasis.L2Norm(q,q_exakt,J,L2modalVnodal);
    //    DGBasis.LinfNorm(q,q_exakt,LinfmodalVnodal);
    //    for(int ik=0;ik<Neq;ik++){
    //        L2modalVnodal[ik]= sqrt(L2modalVnodal[ik]);
    //    }
    //
    //    for(int ik=0;ik<Neq;ik++){
    //        cout <<"L2modalVnodal["<<ik<<"] is: " <<  L2modalVnodal[ik] <<"\n";
    //
    //    }
    //    for(int ik=0;ik<Neq;ik++){
    //        cout <<"LinfmodalVnodal["<<ik<<"] is: " <<  LinfmodalVnodal[ik] <<"\n";
    //
    //    }
    //
    //
    //free(q_modal);



    if(MPI.rank==0)
    {
        cout <<"declaring occa Kernels and Variables... ";
    }


    //    occa::kernel FluxKernel;
    //    occa::memory o_Ftilde,o_Gtilde;

    occa::kernel VolumeKernel;
    occa::kernel VolumeKernelSTD;
    occa::kernel calcNumFluxes;
    occa::kernel SurfaceKernel;
    occa::kernel UpdateKernel;
    occa::kernel calcRK;
    occa::kernel addS;
    occa::kernel CollectEdgeData;
    occa::kernel CollectEdgeData_Bottom;

    occa::kernel CollectEdgeDataGradient;
    occa::kernel calcGradient;
    occa::kernel calcNumFluxesGradient;
    occa::kernel SurfaceKernelGradient;
    occa::kernel VolumeKernelViscose;
    occa::kernel calcNumFluxesViscose;
    occa::kernel ShockCapturing;
    occa::kernel calcDiscBottomSurf;
//    occa::kernel calcEdgeValues;
    occa::kernel preservePosivitity;
    occa::kernel calcAvg;
    occa::kernel FindLambdaMax;
    occa::kernel scaleGradient;
    occa::kernel SurfaceKernelVisc;
    occa::kernel UpdateQt;
	
	occa::kernel VolumeKernelFD;

//    occa::kernel MemCopyKernel;


    occa::memory o_Qtmp; // for SSP RK
    occa::memory o_D,o_Dstrong,o_Dhat,o_Qt,o_gRK,o_q;//,o_Neq,o_ngl,o_Jac;
    occa::memory o_VdmInv;//,o_SubCellMat,;
    occa::memory o_Jac,o_Yxi,o_Yeta,o_Xxi,o_Xeta;

    occa::memory o_SurfaceParts, o_ElemEdgeOrientation, o_ElemToEdge,o_ElemEdgeMasterSlave;
    occa::memory o_nx,o_ny,o_scal;
    occa::memory o_Bx,o_By,o_B;
    occa::memory o_x,o_y;
    occa::memory o_qL,o_qR,o_bL,o_bR;
    occa::memory o_EdgeData;
    occa::memory o_qGradientX, o_qGradientY,o_SurfGradientX,o_SurfGradientY ;
    occa::memory o_qGradientXL,o_qGradientXR,o_qGradientYL,o_qGradientYR;
    occa::memory o_SurfacePartsVisc;
    occa::memory o_EleSizes, o_ViscPara;
    occa::memory o_ViscParaL,o_ViscParaR;

    occa::memory o_DBSurf1,o_DBSurf2;
    occa::memory o_LambdaMax;

    occa::memory o_QtVisc;
    occa::memory o_Qavg;

    occa::memory o_GLw;
	
	occa::memory o_ViscForPlot;
	
	occa::memory o_DcentralFD, o_DupwindFD, o_DdownwindFD;
	occa::memory o_isPartlyDry;

//   occa::memory o_PackSend, o_PackReceive;

    if(MPI.rank==0)
    {
        cout <<"... finished.\n";
        cout <<"Allocating Memory on the device...      ";
    }

    // NEEDED PARAMETERS: (Nelem,Neq,ngl,DGMetrics.Jac,o_F,o_D,o_Qt);
    //o_Ftilde  = device.malloc(NoDofs*sizeof(dfloat));
    //o_Gtilde  = device.malloc(NoDofs*sizeof(dfloat));
    //cout <<"rank: " << MPI.rank << " I AM AT: Device Memory Allocation" << "\n";
    if (rkSSP)
    {
        o_Qtmp = device.malloc(NoDofs*sizeof(dfloat));
    }

	o_DcentralFD  = device.malloc(ngl2*sizeof(dfloat));
	o_DupwindFD  = device.malloc(ngl2*sizeof(dfloat));
	o_DdownwindFD  = device.malloc(ngl2*sizeof(dfloat));

    o_D  = device.malloc(ngl2*sizeof(dfloat));
    o_Dstrong  = device.malloc(ngl2*sizeof(dfloat));
    o_Dhat  = device.malloc(ngl2*sizeof(dfloat));
    o_VdmInv  = device.malloc(ngl2*sizeof(dfloat));
    //o_SubCellMat = device.malloc(ngl2*sizeof(dfloat));
    o_x  = device.malloc(NoSpaceDofs*sizeof(dfloat));
    o_y  = device.malloc(NoSpaceDofs*sizeof(dfloat));
    o_Bx  = device.malloc(NoSpaceDofs*sizeof(dfloat));
    o_By  = device.malloc(NoSpaceDofs*sizeof(dfloat));
    o_B  = device.malloc(NoSpaceDofs*sizeof(dfloat));
    o_Qt = device.malloc(NoDofs*sizeof(dfloat));
    o_q = device.malloc(NoDofs*sizeof(dfloat));
    o_gRK = device.malloc(NoDofs*sizeof(dfloat));


    o_ElemEdgeMasterSlave = device.malloc(4*Nelem*sizeof(int));
    o_ElemEdgeOrientation = device.malloc(4*Nelem*sizeof(int));
    o_ElemToEdge = device.malloc(4*Nelem*sizeof(int));
    o_qL = device.malloc(Neq*ngl*Nfaces*sizeof(dfloat));
    o_qR = device.malloc(Neq*ngl*Nfaces*sizeof(dfloat));
    o_bL = device.malloc(ngl*Nfaces*sizeof(dfloat));
    o_bR = device.malloc(ngl*Nfaces*sizeof(dfloat));
    o_SurfaceParts = device.malloc(Neq*ngl*Nfaces*sizeof(dfloat));


    o_Jac  = device.malloc(NoSpaceDofs*sizeof(dfloat));
    o_Yxi  = device.malloc(NoSpaceDofs*sizeof(dfloat));
    o_Yeta  = device.malloc(NoSpaceDofs*sizeof(dfloat));
    o_Xxi  = device.malloc(NoSpaceDofs*sizeof(dfloat));
    o_Xeta  = device.malloc(NoSpaceDofs*sizeof(dfloat));

    o_nx = device.malloc(ngl*Nfaces*sizeof(dfloat));
    o_ny = device.malloc(ngl*Nfaces*sizeof(dfloat));
    o_scal = device.malloc(ngl*Nfaces*sizeof(dfloat));
    o_EdgeData = device.malloc(8*Nfaces*sizeof(int));


    o_DBSurf1 = device.malloc(Nfaces*ngl*sizeof(dfloat));
    o_DBSurf2 = device.malloc(Nfaces*ngl*sizeof(dfloat));
    //viscose term

    o_LambdaMax = device.malloc(Nelem*sizeof(dfloat));
    //pospres
    o_EleSizes = device.malloc(Nelem*sizeof(dfloat));
	
	o_isPartlyDry= device.malloc(Nelem*sizeof(int));


    if (PositivityPreserving == 1)
    {
        o_GLw  = device.malloc(ngl*sizeof(dfloat));
        o_Qavg = device.malloc(Nelem*4*sizeof(dfloat));

    }
    if (ArtificialViscosity == 1)
    {

        o_qGradientX = device.malloc(NoDofs*sizeof(dfloat));
        o_qGradientY = device.malloc(NoDofs*sizeof(dfloat));
        o_qGradientXL = device.malloc(Neq*ngl*Nfaces*sizeof(dfloat));
        o_qGradientXR = device.malloc(Neq*ngl*Nfaces*sizeof(dfloat));
        o_qGradientYL = device.malloc(Neq*ngl*Nfaces*sizeof(dfloat));
        o_qGradientYR = device.malloc(Neq*ngl*Nfaces*sizeof(dfloat));
        o_SurfGradientX = device.malloc(Neq*ngl*Nfaces*sizeof(dfloat));
        o_SurfGradientY = device.malloc(Neq*ngl*Nfaces*sizeof(dfloat));
        o_SurfacePartsVisc = device.malloc(Neq*ngl*Nfaces*sizeof(dfloat));
        o_QtVisc = device.malloc(NoDofs*sizeof(dfloat));


        o_ViscPara = device.malloc(Nelem*sizeof(dfloat));
		o_ViscForPlot = device.malloc(Nelem*sizeof(dfloat));
		
        o_ViscParaL = device.malloc(Nfaces*sizeof(dfloat));
        o_ViscParaR = device.malloc(Nfaces*sizeof(dfloat));

    }


    if(MPI.rank==0)
    {
        cout <<" ...Memory on Device allocated.\n";
    }




    info.addDefine("procID",MPI.rank);
    info.addDefine("NEpad",NEpad);
    info.addDefine("Nedgepad",Nedgepad);
    info.addDefine("NEsurfpad",NEsurfpad);
    info.addDefine("ngl",ngl);
    info.addDefine("ngl2",ngl2);
    info.addDefine("Neq",Neq);

    info.addDefine("PI",PI);
    //info.addDefine("Nelem",Nelem);
    info.addDefine("dfloat", dfloatString);
    info.addDefine("dfloat4", dfloat4String);
    info.addDefine("ES", ES);
    info.addDefine("Testcase", Testcase);

    info.addDefine("eps0",epsilon_0);
    info.addDefine("sigmaMax",sigma_max);
    info.addDefine("sigmaMin",sigma_min);
    int nglPad=0;
    if (ngl%4==0)
    {
        nglPad=1;
    }
    info.addDefine("nglPad",nglPad );
    dfloat TOL_PosPres = PosPresTOL;//pow(10.0,-4);
    dfloat ZeroTOL = pow(10.0,-5);
    dfloat geomface = 1.0/DGBasis.w_GL[0];
    dfloat zero = 0.0;
    dfloat half = 0.5;
    dfloat one = 1.0;
    dfloat fourth = 0.25;
	dfloat eight = 8.0;
	dfloat two = 2.0;
	dfloat onepointfive = 1.5;
    dfloat halfg = half*g_const;
    dfloat fourthg = fourth*g_const;
    info.addDefine("zero",zero );
    info.addDefine("half",half );
    info.addDefine("fourth",fourth );
    info.addDefine("one",one );
	info.addDefine("two",two );
	info.addDefine("eight",eight );
	info.addDefine("onepointfive",onepointfive );
	
    info.addDefine("half_g",halfg);
    info.addDefine("fourth_g",fourthg);
    info.addDefine("g_const",g_const);
    info.addDefine("geomFace",geomface );
    info.addDefine("PosPresTOL",TOL_PosPres);
    info.addDefine("ZeroTOL",ZeroTOL);
    //copy all permanent data onto the device
    if(MPI.rank==0)
    {
        cout <<"Copy Necessary Data onto Device...";
    }
    o_nx.copyFrom(normalsX);
    o_ny.copyFrom(normalsY);
    o_scal.copyFrom(Scal);
    o_Yxi.copyFrom(y_xi);
    o_Yeta.copyFrom(y_eta);
    o_Xxi.copyFrom(x_xi);
    o_Xeta.copyFrom(x_eta);
    o_B.copyFrom(b);
    o_Bx.copyFrom(Bx);
    o_By.copyFrom(By);
    o_D.copyFrom(Dmat);
	
    o_Dstrong.copyFrom(DGBasis.Dstrong);
    o_Dhat.copyFrom(Dhat);
    o_Jac.copyFrom(J);
    o_ElemEdgeMasterSlave.copyFrom(ElemEdgeMasterSlave);
    o_ElemEdgeOrientation.copyFrom(ElemEdgeOrientation);
    o_ElemToEdge.copyFrom(ElemToEdge);
    o_EdgeData.copyFrom(EdgeData);
    o_x.copyFrom(x_phy);
    o_y.copyFrom(y_phy);
    o_q.copyFrom(q);
    o_EleSizes.copyFrom(ElementSizes);
    o_gRK.copyFrom(gRK);
    o_Qt.copyFrom(Qt);
    o_VdmInv.copyFrom(VdmInv);
    //o_SubCellMat.copyFrom(SubCellMat);
	o_DcentralFD.copyFrom(DCentralFD);
	o_DupwindFD.copyFrom(DupwindFD);
	o_DdownwindFD.copyFrom(DdownwindFD);

    dfloat * qavgtmp = (dfloat*) calloc(Nelem*4,sizeof(dfloat));
    if(PositivityPreserving == 1)
    {
        o_GLw.copyFrom(GLw);

        o_Qavg.copyFrom(qavgtmp);

    }
    free(qavgtmp);
    if(MPI.rank==0)
    {
        cout << "         done.\n";
    }





    free(gRK);
    free(x_xi);
    free(y_xi);
    free(x_eta);
    free(y_eta);
    free(Bx);
    free(By);
    free(EdgeData);
    free(ElemToEdge);
    free(ElemEdgeMasterSlave);
    free(ElemEdgeOrientation);
    free(Scal);
    free(normalsX);
    free(normalsY);

    free(x_phy);
    free(y_phy);




    //cout <<"rank: " << MPI.rank << " I AM AT: Building Kernels" << "\n";

    if(MPI.rank==0)
    {
        cout <<"Build Kernels...\n";
    }

    switch(Testcase)
    {
    case 1:  // THIS INCLUDES DIRICHLET BOUNDARIES FOR PERIODIC CONVERGENCE TEST
    {
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/Dirichlet_ConvTest.okl","CollectEdgeData",info);
        addS = device.buildKernelFromSource("okl/ManufacturedSolutions/S_ConvTest.okl","addS",info);
        break;
    }
    case 32:  // Inflow Boundaries for 3 Mound PP test case
    {
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/3MoundInflow.okl","CollectEdgeData",info);
        break;
    }
    default:
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/SolidWalls.okl","CollectEdgeData",info);
        break;
    }
    //cout <<"rank: " << MPI.rank << " I AM AT: Collecting Kernels Post TESTCASE" << "\n";
    CollectEdgeData_Bottom=device.buildKernelFromSource("okl/GatherEdgeData/CollectEdgeData_Bottom.okl","CollectEdgeData_Bottom",info);



    if(Fluxdifferencing)
    {

        std::ostringstream oss;
        cout << "Kernel Version: V " << KernelVersion << ".\n";
        oss << "okl/DG/VolumeKernelFluxDiffV" << KernelVersion << ".okl";
        std::string var = oss.str();
        VolumeKernel=device.buildKernelFromSource(var,"VolumeKernelFluxDiff",info);

    }
    else
    {
        std::ostringstream oss;
        cout << "Kernel Version: V " << KernelVersionSTD << ".\n";
        oss << "okl/DG/VolumeKernelV" << KernelVersionSTD << ".okl";
        std::string var = oss.str();
        VolumeKernel=device.buildKernelFromSource(var,"VolumeKernel",info);
    }

	// FD = FINITE DIFFERENCE HERE
	VolumeKernelFD=device.buildKernelFromSource("okl/NEW_OPERATOR/VolumeKernelFD.okl","VolumeKernelFD",info);

    switch(NumFlux)
    {
    // ENTROPY STABLE FLUX FROM PAPERS
    case 0:
    {
        calcNumFluxes=device.buildKernelFromSource("okl/RiemannSolvers/calcESFlux.okl","calcNumFluxes",info);
        break;
    }
    //Lax Friedrich Flux (MAXIMUM EIGENVALUES ON EDGE ARE EVALUATED POINT WISE
    case 1:
    {
        calcNumFluxes=device.buildKernelFromSource("okl/RiemannSolvers/calcLaxFriedrich.okl","calcNumFluxes",info);
        break;
    }
    //Lax Friedrich Type Entropy Stable Flux   // NOT WORKING PROPERLY ATM
    case 2:
    {
        calcNumFluxes=device.buildKernelFromSource("okl/RiemannSolvers/calcLFTypeESFlux.okl","calcNumFluxes",info);
        break;
    }
	case 3:
    {
        calcNumFluxes=device.buildKernelFromSource("okl/RiemannSolvers/calcESFluxRotated.okl","calcNumFluxes",info);
        break;
    }
    }

    // calc specific value on edges like jump in b and average h
//    calcEdgeValues          =   device.buildKernelFromSource("okl/DiscontinuousBathimetry/calcEdgeValues.okl","calcEdgeValues",info);
    // adds additional surface terms due to a possibly discontinuous bottom topography
    calcDiscBottomSurf      =   device.buildKernelFromSource("okl/DiscontinuousBathimetry/calcDiscBottomSurf.okl","calcDiscBottomSurf",info);
    // standard dg kernel for the surface parts
   SurfaceKernel=device.buildKernelFromSource("okl/DG/SurfaceKernel.okl","SurfaceKernel",info);
    //kernel to compute eigenvalues
    FindLambdaMax           =   device.buildKernelFromSource("okl/DG/FindLambdaMax.okl","FindLambdaMax",info);

    if (rkSSP)
    {
        UpdateKernel=device.buildKernelFromSource("okl/RungeKutta/UpdateKernelSSP_V1.okl","UpdateKernel",info);
    }
    else
    {
        // no SSP runge kutta but Low Storage verison
        calcRK=device.buildKernelFromSource("okl/RungeKutta/rkLS.okl","calcRK",info);
        UpdateKernel=device.buildKernelFromSource("okl/RungeKutta/UpdateKernelLS.okl","UpdateKernel",info);
    }




    if (ArtificialViscosity)
    {
        UpdateQt                =   device.buildKernelFromSource("okl/ViscoseParts/UpdateQt.okl","UpdateQt",info);
        CollectEdgeDataGradient =   device.buildKernelFromSource("okl/BR1_Gradient/CollectEdgeDataGradient.okl","CollectEdgeDataGradient",info);
        SurfaceKernelGradient   =   device.buildKernelFromSource("okl/BR1_Gradient/SurfaceKernelGradient.okl","SurfaceKernelGradient",info);
        calcGradient            =   device.buildKernelFromSource("okl/BR1_Gradient/calcGradient.okl","calcGradient",info);
        calcNumFluxesGradient   =   device.buildKernelFromSource("okl/BR1_Gradient/calcNumFluxesGradient.okl","calcNumFluxesGradient",info);
        scaleGradient          =   device.buildKernelFromSource("okl/BR1_Gradient/scaleGradient.okl","scaleGradient",info);
        VolumeKernelViscose     =   device.buildKernelFromSource("okl/ViscoseParts/VolumeKernelViscose.okl","VolumeKernelViscose",info);
        calcNumFluxesViscose    =   device.buildKernelFromSource("okl/ViscoseParts/calcNumFluxesViscose.okl","calcNumFluxesViscose",info);
        SurfaceKernelVisc       =   device.buildKernelFromSource("okl/ViscoseParts/SurfaceKernelVisc.okl","SurfaceKernelVisc",info);
        ShockCapturing          =   device.buildKernelFromSource("okl/ViscoseParts/ShockCapturing.okl","ShockCapturing",info);
    }
    if (PositivityPreserving)
    {
        calcAvg      =   device.buildKernelFromSource("okl/Positivity/calcAvg.okl","calcAvg",info);
        preservePosivitity      =   device.buildKernelFromSource("okl/Positivity/PosPres.okl","PosPres",info);
    }


//    MemCopyKernel           =   device.buildKernelFromSource("okl/DG/MemCopyComparison.okl","MemCopyComparison",info);



    cout <<"rank: " << MPI.rank << " GB allocated on device: " << device.bytesAllocated()/(1024.f*1024*1024) << "\n";





    dfloat * mCheckpoints;

    if (NumPlots>0)
    {
	if (Testcase == 32){
		NumPlots=5;      
    	}
	if (Testcase == 31){
		NumPlots=5;
	}
        mCheckpoints = (dfloat*) calloc(NumPlots,sizeof(dfloat));
        mCheckpoints[0]=0.0;
        if (NumPlots>1)
        {
            dfloat mTimestep = T/(NumPlots-1);
            for(int i=1; i<NumPlots; i++)
            {
                mCheckpoints[i]=i*mTimestep;

            }
        }

	if (Testcase == 32){
		mCheckpoints[1]=8.0;
		mCheckpoints[2]=30.0;
		mCheckpoints[3]=300.0;
		mCheckpoints[4]=900.0;    
    	}
	if (Testcase == 31){
		mCheckpoints[1]=T/12.0;
		mCheckpoints[2]=T/6.0;
		mCheckpoints[3]=T/4.0;
		mCheckpoints[4]=T;    
    	}
    }
    

    dfloat * tCheckpoints;
    int timeCount=0;
    if (NumTimeChecks>0)
    {
        tCheckpoints = (dfloat*) calloc(NumTimeChecks,sizeof(dfloat));
        tCheckpoints[0]=0.0;
        dfloat tTimestep = T/(NumTimeChecks-1);
        for(int i=1; i<NumTimeChecks; i++)
        {
            tCheckpoints[i]=i*tTimestep;

        }
    }





    // Exchange information about the bottom topography only ONCE!
    CollectEdgeData_Bottom(Nfaces,o_EdgeData,o_B,o_bL,o_bR);
    o_bL.copyTo(bL);
    o_bR.copyTo(bR);
    CollectEdgeDataMPI_Bonly(MPI, DGMeshPartition, bL,  bR);
    MPI_Waitall(DGMeshPartition.NumProcessors,MPI.Recv_b_reqs, MPI.stats);
    MPI_Waitall(DGMeshPartition.NumProcessors,MPI.Send_b_reqs, MPI.stats);

    o_bL.copyFrom(bL);
    o_bR.copyFrom(bR);




    //printing initial conditions
    int plotCount=0;
    if(plotCount<NumPlots)
    {
        if(t>=mCheckpoints[plotCount])
        {
            o_q.copyTo(q);

            if(MPI.rank==0)
            {

                CollectSolution( MPI, DGMeshPartition, q, Q_global);

		if (Testcase == 31){
			dfloat * q_exakt = (dfloat*) calloc(NoDofs_global,sizeof(dfloat));
			InitQ(0,DGMeshPartition,Testcase,Nelem_global,ngl,ngl2,x_phy_global,y_phy_global,q_exakt,t,b_global,g_const);
			PlotSolutionWithExact(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,Q_global,b_global,plotCount,q_exakt);
			free(q_exakt);
		}else{
			PlotSolution(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,Q_global,b_global,plotCount);
		}

				if (EntropyPlot){
					dfloat TotalEntropy=0.0;
					DGBasis.calcTotalEntropy(g_const,Q_global,b_global,J_global,&TotalEntropy);
					EntropyOverTime[plotCount] = TotalEntropy;
					dfloat TotalMass=0.0;
					DGBasis.calcTotalMass(Q_global,J_global,&TotalMass);
					MassOverTime[plotCount] = TotalMass;
					EntropyTimes[plotCount] = t;
				}
            }
            else
            {
                SendSolution(MPI, DGMeshPartition,  q);
            }

            plotCount=plotCount+1;

        }
    }



    if(MPI.rank==0)
    {
        cout << "Entering Time Loop!\n";
    }
    //time run
    struct timeval startTV, endTV;
    if(MPI.rank==0)
    {

        gettimeofday(&startTV, NULL);
    }



    // cout << "rank: " << MPI.rank << " Entering Time Loop!\n";

    //TIME LOOP!!
    dfloat globalLambdaMax=0.0;
    dfloat maxViscPara=0.0;
    dfloat dt_i=0.0;
    dfloat dt_v=0.0;
    dfloat LocalLambdas[DGMeshPartition.NumElements];
    dfloat rkD=0.0;

    cout << "Minimum Ele Size: " << globalMinEleSize << ".\n";
    while (t<T)
    {
        FindLambdaMax(Nelem, o_q, o_LambdaMax);
        globalLambdaMax=0.0;
        o_LambdaMax.copyTo(LocalLambdas);
        GetGlobalLambdaMax(MPI,  DGMeshPartition,LocalLambdas, &globalLambdaMax);
	dt_i = globalMinEleSize/(ngl) * CFL /globalLambdaMax;
	if (Testcase==32){
		if (t==0.0){
			dt_i = 0.0001;		
		}
	}
	if ( ArtificialViscosity==1)
        {
            ShockCapturing(Nelem, o_q,o_VdmInv,o_EleSizes,o_ViscPara,o_ViscForPlot,o_isPartlyDry);
            o_ViscPara.copyTo(ViscPara);
            GetGlobalViscParaMax(MPI,  DGMeshPartition,ViscPara, &maxViscPara);
//	   cout << "Visc para max: " << maxViscPara <<"\n";
            dt_v = DFL/(pow(ngl,2)) * pow(globalMinEleSize,2) / maxViscPara;
	    dt=fmin(T-t,fmin(dt_i,dt_v));
        }else{
	    dt=fmin(T-t,dt_i);
	}


//cout << "timestep " << dt << "\n" ;
        if (rkSSP)
        {
            o_Qtmp.copyFrom(o_q);
        }

        for (int rkstage=0; rkstage<RK.rkStages; rkstage++)
        {

            dfloat rkA=RK.CoeffsA[rkstage];
            dfloat rkB=RK.CoeffsB[rkstage];
            dfloat rkC=RK.CoeffsC[rkstage];
            dfloat intermediatetime=t;
            if (rkSSP)
            {
                switch(rkstage)
                {
                case 0:
                    rkD=0.0;
                    break;
                case 1:
                    rkD=1.0;
                    break;
                case 2:
                    rkD=1.0/2.0;
                    break;
                }
                intermediatetime = t +rkD*dt;
            }
            else
            {
                intermediatetime = t +rkB*dt;
            }

            CollectEdgeData(Nfaces,o_EdgeData,o_q,o_x,o_y,o_nx,o_ny, o_bL,o_bR, o_qL, o_qR, intermediatetime);
            o_qL.copyTo(qL);
            o_qR.copyTo(qR);
            CollectEdgeDataMPI(MPI, DGMeshPartition, qL, qR);







			VolumeKernel(Nelem, o_Jac,o_Yxi,o_Yeta,o_Xxi,o_Xeta,o_q,o_D,o_Bx,o_By,o_Qt);

			o_Qt.copyTo(Qt);
           cout <<"\n q_t ESDGSEM : \n";
			for (int ie=0;ie<Nelem;ie++){
					cout <<"Ele: " << ie <<"\n";
				for(int j=0;j<ngl;++j){
					for(int i=0;i<ngl;++i){
						int id = ie*ngl2*Neq +  j*ngl+i;
					   cout <<Qt[id]<<"  ";
				  }
					cout <<"\n";
				}
			}

			VolumeKernelFD(Nelem, o_Jac,o_Yxi,o_Yeta,o_Xxi,o_Xeta,o_q,o_isPartlyDry,o_DcentralFD,o_DupwindFD,o_DdownwindFD,o_Bx,o_By,o_Qt);
			o_Qt.copyTo(Qt);
		        cout <<"\n q_t NEW : \n";
			for (int ie=0;ie<Nelem;ie++){
					cout <<"Ele: " << ie <<"\n";
				for(int j=0;j<ngl;++j){
					for(int i=0;i<ngl;++i){
						int id = ie*ngl2*Neq +  j*ngl+i;
					   cout <<Qt[id]<<"  ";
				  }
					cout <<"\n";
				}
			}


            MPI_Waitall(DGMeshPartition.NumProcessors,MPI.Recv_q_reqs, MPI.stats);
            MPI_Waitall(DGMeshPartition.NumProcessors,MPI.Send_q_reqs, MPI.stats);
            o_qL.copyFrom(qL);
            o_qR.copyFrom(qR);

            calcNumFluxes(Nfaces,o_nx,o_ny,o_scal,o_qL,o_qR,o_bL,o_bR,o_SurfaceParts);
            calcDiscBottomSurf(Nfaces,o_qL,o_qR, o_bL,o_bR,o_nx,o_ny,o_scal,o_DBSurf1,o_DBSurf2);

            SurfaceKernel(Nelem,o_Jac,o_ElemEdgeMasterSlave,o_ElemEdgeOrientation,o_ElemToEdge, o_SurfaceParts,o_DBSurf1,o_DBSurf2,o_Qt);

            if ( ArtificialViscosity==1)
            {


                // at first RK step we already now o_ViscPara from time step computation!
                if (rkstage>0)
                {
                    ShockCapturing(Nelem, o_q,o_VdmInv,o_EleSizes,o_ViscPara,o_ViscForPlot,o_isPartlyDry);
                }

                calcNumFluxesGradient(Nfaces,o_EdgeData,o_nx,o_ny,o_scal, o_qL, o_qR, o_bL,o_bR, o_SurfGradientX,o_SurfGradientY);

                calcGradient(Nelem, o_Jac,o_Yxi,o_Yeta,o_Xxi,o_Xeta,o_q,o_B,o_Dhat,o_qGradientX,o_qGradientY);


                SurfaceKernelGradient(Nelem,o_Jac,o_ElemEdgeMasterSlave,o_ElemEdgeOrientation,o_ElemToEdge, o_SurfGradientX, o_SurfGradientY,o_qGradientX,o_qGradientY);

                scaleGradient(Nelem,o_q,o_qGradientX,o_qGradientY);

                CollectEdgeDataGradient(Nfaces,o_EdgeData,o_qGradientX,o_qGradientY,o_ViscPara,o_ViscParaL,o_ViscParaR, o_qGradientXL, o_qGradientXR,o_qGradientYL, o_qGradientYR);

                o_ViscParaL.copyTo(ViscParaL);
                o_ViscParaR.copyTo(ViscParaR);
                o_qGradientXL.copyTo(qGradientXL);
                o_qGradientXR.copyTo(qGradientXR);
                o_qGradientYL.copyTo(qGradientYL);
                o_qGradientYR.copyTo(qGradientYR);
                    
                // not an OCCA kernel, this is MPI communication for exchanging gradient data
				CollectViscoseEdgeDataMPI(MPI, DGMeshPartition, ViscParaL, ViscParaL, qGradientXL,  qGradientXR,qGradientYL,qGradientYR);
                VolumeKernelViscose(Nelem, o_Jac,o_Yxi,o_Yeta,o_Xxi,o_Xeta,o_qGradientX,o_qGradientY,o_Dstrong,o_ViscPara,o_QtVisc);



                MPI_Waitall(DGMeshPartition.NumProcessors,MPI.Recv_qX_reqs, MPI.stats);
                MPI_Waitall(DGMeshPartition.NumProcessors,MPI.Recv_qY_reqs, MPI.stats);
                MPI_Waitall(DGMeshPartition.NumProcessors,MPI.Recv_ViscPar_reqs, MPI.stats);
                MPI_Waitall(DGMeshPartition.NumProcessors,MPI.Send_qX_reqs, MPI.stats);
                MPI_Waitall(DGMeshPartition.NumProcessors,MPI.Send_qY_reqs, MPI.stats);
                MPI_Waitall(DGMeshPartition.NumProcessors,MPI.Send_ViscPar_reqs, MPI.stats);
                o_ViscParaL.copyFrom(ViscParaL);
                o_ViscParaR.copyFrom(ViscParaR);
                o_qGradientXL.copyFrom(qGradientXL);
                o_qGradientXR.copyFrom(qGradientXR);
                o_qGradientYL.copyFrom(qGradientYL);
                o_qGradientYR.copyFrom(qGradientYR);

                calcNumFluxesViscose(Nfaces,o_EdgeData,o_nx,o_ny,o_scal,o_ViscParaL,o_ViscParaR,o_qGradientXL,o_qGradientXR,o_qGradientYL,o_qGradientYR,o_SurfacePartsVisc);

                SurfaceKernelVisc(Nelem,o_Jac,o_ElemEdgeMasterSlave,o_ElemEdgeOrientation,o_ElemToEdge, o_SurfacePartsVisc,o_QtVisc);



                UpdateQt(Nelem,o_QtVisc,o_Qt);

            }

            // add manufactured source term for convergence test
            if (Testcase==1)
            {
                addS(Nelem,o_Bx,o_By,o_B,o_x,o_y,intermediatetime,o_Qt);
            }

            if (rkSSP)
            {
                UpdateKernel(Nelem,rkA,rkB,rkC,dt,o_Qt,o_Qtmp,o_q);
            }
            else
            {
                calcRK(Nelem,o_Qt,rkA,rkB,o_gRK);
                UpdateKernel(Nelem,rkC,dt,o_gRK,o_q);
            }

            if(PositivityPreserving==1)
            {
//                calcAvg(Nelem,o_EleSizes,o_GLw, o_Jac,o_q,o_Qavg);

//                dfloat * qavgtmp = (dfloat*) calloc(Nelem*4,sizeof(dfloat));
////                o_Qavg.copyTo(qavgtmp);
////                for (int i=0; i<1;i++){
////
////                    cout << "Ele: " << i << " ";
////                    cout << "Avg H: " << qavgtmp[i*4] << " ";
////                    cout << "Avg Hu: " << qavgtmp[i*4+1] << " ";
////                    cout << "Avg Hv: " << qavgtmp[i*4+2] << " ";
////                    cout << "Min H: " << qavgtmp[i*4+3] << " \n";
////                }
////                free(qavgtmp);




//                preservePosivitity(Nelem,o_Qavg,o_q);

//                // old kernel header
               preservePosivitity(Nelem,o_EleSizes,o_GLw, o_Jac,o_q);
//
            }




        }



        t += dt;

        //PRINT SOLUTION
        if(plotCount<NumPlots)
        {
            if(t>=mCheckpoints[plotCount])
            {

                o_q.copyTo(q);
                if(ArtificialViscosity==1)
                {
//                    o_qGradientX.copyTo(Qx);
//                    o_QtVisc.copyTo(QtVisc);
//                    o_qGradientY.copyTo(Qy);
//                    o_ViscPara.copyTo(ViscPara);
                }
                if(MPI.rank==0)
                {
                    cout << "we have time: "<<t<< " and we do plot number " <<  plotCount <<"!\n";
                    CollectSolution( MPI, DGMeshPartition, q, Q_global);
                    
					
					if (Testcase == 31){
						dfloat * q_exakt = (dfloat*) calloc(NoDofs_global,sizeof(dfloat));
						InitQ(0,DGMeshPartition,Testcase,Nelem_global,ngl,ngl2,x_phy_global,y_phy_global,q_exakt,t,b_global,g_const);
						PlotSolutionWithExact(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,Q_global,b_global,plotCount,q_exakt);
						
						free(q_exakt);
					}else{
						PlotSolution(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,Q_global,b_global,plotCount);
					}
                    if(ArtificialViscosity==1)
                    {
					
					  o_ViscForPlot.copyTo(ViscPara);
                      CollectViscPara(MPI,   DGMeshPartition, ViscPara, ViscPara_Global);
                      PlotViscoseParameter(Nelem_global, ngl, x_phy_global,y_phy_global, ViscPara_Global, plotCount);
//                        CollectViscosity( MPI, DGMeshPartition, Qx,Qy, Qx_global, Qy_global);
//                        PlotViscosity(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,Qx_global,Qy_global,plotCount);
//                        CollectViscosity( MPI, DGMeshPartition, QtVisc,Qy, QtVisc_global, Qy_global);
//                        PlotViscosity(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,QtVisc_global,Qy_global,plotCount);

                    }

				if (EntropyPlot){
					dfloat TotalEntropy=0.0;
					DGBasis.calcTotalEntropy(g_const,Q_global,b_global,J_global,&TotalEntropy);
					EntropyOverTime[plotCount] = TotalEntropy;
					dfloat TotalMass=0.0;
					DGBasis.calcTotalMass(Q_global,J_global,&TotalMass);
					MassOverTime[plotCount] = TotalMass;
					EntropyTimes[plotCount] = t;
				}
                }
                else
                {
                    SendSolution(MPI, DGMeshPartition,  q);
                    if(ArtificialViscosity==1)
                    {
						o_ViscForPlot.copyTo(ViscPara);
                        SendViscPara(MPI, DGMeshPartition,  ViscPara);
                        //            SendViscosity(MPI, DGMeshPartition,  Qx,Qy);
//                        SendViscosity(MPI, DGMeshPartition,  QtVisc,Qy);
                    }
                }
//
                plotCount=plotCount+1;
            }
        }

        if (MPI.rank==0)
        {
            if (timeCount < NumTimeChecks)
                if(t>=tCheckpoints[timeCount])
                {
                    cout << "time: " << t << " last timestep "<< dt<< " dt_i = "<< dt_i << " dt_v " << dt_v << "\n";
                    timeCount=timeCount+1;
                }
        }

    }
    o_q.copyTo(q);



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////ANALYSIS OF SOLUTION ///////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    if(MPI.rank==0)
    {
        gettimeofday(&endTV, NULL);

        double TimeDifferenceMeasured = ((endTV.tv_sec  - startTV.tv_sec) * 1000000u + endTV.tv_usec - startTV.tv_usec) / 1.e6;


        //endTV.tv_sec  - startTV.tv_sec +  (endTV.tv_usec - startTV.tv_usec) / 1000000 ;

        cout <<"time taken in seconds = "<< TimeDifferenceMeasured<< "\n";
    }




    if(MPI.rank==0)
    {

        cout << "I am rank: " << MPI.rank << " and i collect solution now!\n";
        CollectSolution( MPI, DGMeshPartition, q, Q_global);
    }
    else
    {
        cout << "I am rank: " << MPI.rank << " and i send solution now!\n";
        SendSolution(MPI, DGMeshPartition,  q);

    }


    if(MPI.rank==0)
    {
        dfloat * q_exakt = (dfloat*) calloc(NoDofs_global,sizeof(dfloat));

		if (EntropyPlot){
			PlotEntropy(NumPlots, EntropyTimes,EntropyOverTime);
			PlotMass(NumPlots, EntropyTimes,MassOverTime);
		}

        InitQ(0,DGMeshPartition,Testcase,Nelem_global,ngl,ngl2,x_phy_global,y_phy_global,q_exakt,0,b_global,g_const);


        dfloat EntropyDelta=0.0;
        dfloat relEntropyDelta=0.0;
        DGBasis.calcEntropyDelta(g_const,Q_global,q_exakt,b_global,J_global,&EntropyDelta,&relEntropyDelta);
        cout <<"Entropydifference is: " <<  EntropyDelta <<"\n";
        cout <<"relative Entropydifference is: " <<  relEntropyDelta <<"\n";

		dfloat MassDelta=0.0;
		dfloat relMassError=0.0;
		DGBasis.checkConservation(Q_global,q_exakt,J_global,&MassDelta,&relMassError);
		cout <<"Mass difference is: " <<  MassDelta <<"\n";
		cout <<"relative mass difference is: " <<  relMassError <<"\n";

        InitQ(0,DGMeshPartition,Testcase,Nelem_global,ngl,ngl2,x_phy_global,y_phy_global,q_exakt,T,b_global,g_const);



        dfloat L2Error[Neq];
        dfloat LinfError[Neq];
        for(int ik=0; ik<Neq; ik++)
        {
            L2Error[ik]=0.0;
            LinfError[ik]=0.0;
        }
        DGBasis.L2Norm(Q_global,q_exakt,J_global,L2Error);
        DGBasis.LinfNorm(Q_global,q_exakt,LinfError);
        for(int ik=0; ik<Neq; ik++)
        {
            L2Error[ik]= sqrt(L2Error[ik]);
        }

        for(int ik=0; ik<Neq; ik++)
        {
            cout <<"L2Error["<<ik<<"] is: " <<  L2Error[ik] <<"\n";

        }
        for(int ik=0; ik<Neq; ik++)
        {
            cout <<"LinfError["<<ik<<"] is: " <<  LinfError[ik] <<"\n";

        }

        free(Q_global);
        free(q_exakt);
        free(J_global);
        free(b_global);
        free(x_phy_global);
        free(y_phy_global);
    }




    free(q);
    free(Qt);
    free(b);


    free(qL);
    free(qR);
    free(bL);
    free(bR);
    free(ViscParaL);
    free(ViscParaR);
    free(qGradientXL);
    free(qGradientXR);
    free(qGradientYL);
    free(qGradientYR);

    if (NumPlots>0)
    {
        free(mCheckpoints);
    }
    if (NumTimeChecks>0)
    {
        free(tCheckpoints);
    }




    free(SurfaceParts);


    if (rkSSP)
    {
        o_Qtmp.free();
    }


    o_D.free();
    o_Dhat.free();
//o_VdmInv  = device.malloc(ngl2*sizeof(dfloat));
//o_SubCellMat = device.malloc(ngl2*sizeof(dfloat));
    o_x.free();
    o_y.free();
    o_Bx.free();
    o_By.free();
    o_B.free();
    o_Qt.free();
    o_q.free();
    o_gRK.free();


    o_ElemEdgeMasterSlave.free();
    o_ElemEdgeOrientation.free();
    o_ElemToEdge.free();
    o_qL.free();
    o_qR.free();
    o_bL.free();
    o_bR.free();
    o_SurfaceParts.free();


    o_Jac.free();
    o_Yxi.free();
    o_Yeta.free();
    o_Xxi.free();
    o_Xeta.free();

    o_nx.free();
    o_ny.free();
    o_scal.free();
    o_EdgeData.free();


    o_DBSurf1.free();
    o_DBSurf2.free();
	
	
	o_isPartlyDry.free();
	o_DcentralFD.free();
	o_DupwindFD.free();
	o_DdownwindFD.free();

//viscose term

    o_LambdaMax.free();
//pospres
    o_EleSizes.free();


    if (PositivityPreserving == 1)
    {
        o_GLw.free();

    }
    if (ArtificialViscosity == 1)
    {

        o_qGradientX.free();
        o_qGradientY.free();
        o_qGradientXL.free();
        o_qGradientXR.free();
        o_qGradientYL.free();
        o_qGradientYR.free();
        o_SurfGradientX.free();
        o_SurfGradientY.free();
        o_SurfacePartsVisc.free();

        o_ViscPara.free();
		o_ViscForPlot.free();
        o_ViscParaL.free();
        o_ViscParaR.free();
        if (PositivityPreserving == 1)
        {
            o_Qavg.free();

        }
    }
//    o_PackSend.free();
//    o_PackReceive.free();

// done with MPI
    MPI_Finalize();

    return 0;
}






