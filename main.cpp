#include <iostream>
#include <fstream>
#include "occa.hpp"
#include "Constants.hpp"
#include "MPI_setup.h"
#include "basis.h"
#include "SW2D.h"
#include "okada.h"
#include "RungeKutta.h"
#include "Mesh.h"
#include "ReadInput.h"
#include "plots.h"
#include "MeshPartitioning.h"
#include "MPI_Communication.h"
#include <ctime>
#include <sys/time.h>
#include <deviceclass.h>


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
        MPI_Finalize();
        exit(-1);
    }
    deviceclass occa_device(MPI.rank,argc,argv);








    //occa::device device;










    cout.precision(17);
    //initialize polynomial basis

    int ngl, ngl2;;
    int NoDofs,NoSpaceDofs;
    int Nfaces, Nelem;
    int Nfaces_global, Nelem_global,NoDofs_global, NoGradDofs_global,NoSpaceDofs_global;
    int PlotVar;
    int EntropyPlot;
    int ArtificialViscosity;
    int PositivityPreserving;
    int DiscBottom;
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
    int NavgPad;
    int KernelVersion=-1;
    int KernelVersionSTD=-1;
    int ReadInBottom;
    int PartialDryTreatment;
    int FrictionTerms;


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
                      &NavgPad,
                      &KernelVersion,
                      &KernelVersionSTD,
                      &DiscBottom,
                      &ReadInBottom,
                      &PartialDryTreatment,
                      &FrictionTerms);
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
                   &NavgPad,
                   &KernelVersion,
                   &KernelVersionSTD,
                   &DiscBottom,
                   &PartialDryTreatment,
                   &FrictionTerms);

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
        cout <<"DiscBottom: " <<DiscBottom<<"\n";
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
        cout <<"PartialDryTreatment: " <<PartialDryTreatment <<"\n";
        cout <<"FrictionTerms: " <<FrictionTerms <<"\n";
        cout <<"Setting up Basis for N=" << N <<".\n";
    }


    basis DGBasis(N,Fluxdifferencing);

    SW2D SW_Problem(Testcase,PosPresTOL);



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

    Mesh DGMesh(DGBasis.x_GL,ngl,NelemX,ReadInBottom);



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
        NoGradDofs_global=ngl2*Nelem_global*(Neq-1);
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



    dfloat * ViscPara_Global;


    dfloat * Q_global ;

    dfloat * QtVisc_global;
    dfloat * Qx_global;
    dfloat * Qy_global ;

    dfloat * J_global;
    dfloat * b_global;

    dfloat * x_phy_global;
    dfloat * y_phy_global;
    dfloat h_0;


    if (MPI.rank ==0 )
    {
        ViscPara_Global = (dfloat*) calloc(Nelem_global,sizeof(dfloat));


        Q_global = (dfloat*) calloc(NoDofs_global,sizeof(dfloat));
        if (ArtificialViscosity==1)
        {
            Qx_global = (dfloat*) calloc(NoGradDofs_global,sizeof(dfloat));
            Qy_global = (dfloat*) calloc(NoGradDofs_global,sizeof(dfloat));
            QtVisc_global = (dfloat*) calloc(NoGradDofs_global,sizeof(dfloat));
        }

        J_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));
        b_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));

        x_phy_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));


        y_phy_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));


//            cout << "Writing Out Full Mesh ...\n";
//            WriteFullMesh(NoSpaceDofs_global, DGMesh.x_global,DGMesh.y_global);
//            cout << "... DONE\n";


        if (Cartesian)
        {
            if (ReadInBottom)
            {
                ReadFullMesh(NoSpaceDofs_global, b_global, &h_0);
                cout << " Water displacement: " << h_0 << "\n";
            }
        }
        for(int ie=0; ie<Nelem_global; ++ie)
        {

            for(int j=0; j<ngl; ++j)
            {
                for(int i=0; i<ngl; ++i)
                {
                    int id = ie*ngl2   +j*ngl+i;
                    int Qid =ie*ngl2*Neq + j*ngl+i;

                    J_global[id] = 1.0/DGMesh.J_global[id];
                    x_phy_global[id] = DGMesh.x_global[id];
                    y_phy_global[id] = DGMesh.y_global[id];

                    if (ReadInBottom)
                    {
                        if (!Cartesian)
                        {
                            dfloat b_min = -7.3100;// in kilometers now
                            h_0 = -b_min;
                            b_global[id]    =   h_0 + DGMesh.b_global[id];
                        }
                        else
                        {
                           // b_global[id]    =   h_0 + b_global[id];
                        }

                    }



                }
            }
        }





        if (!ReadInBottom)
        {
            cout << "Not reading bottom from file, taking initial data as prescribed for test case!\n";
            SW_Problem.InitB(0,DGMeshPartition,Nelem_global,ngl,ngl2,x_phy_global,y_phy_global,b_global);
        }




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
    int * EdgeReversed = (int*) calloc(Nfaces,sizeof(int));
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

    int * isDryElement = (int*) calloc(Nelem,sizeof(int));

    dfloat * q = (dfloat*) calloc(NoDofs,sizeof(dfloat));


    if(MPI.rank==0)
    {
        cout <<"... finished.\n";
    }



    if(MPI.rank==0)
    {
        cout <<"Distributing bottom topography data to MPI-ranks \n";
        DGMeshPartition.DivideBottom(MPI,NoSpaceDofs,b_global,b);
        //std::memcpy(&b, &b_global, sizeof b_global);
        //std::memcpy(b, b_global, NoSpaceDofs*sizeof(dfloat));


        //cout <<" b at random node " << b[50] << " b_global: " << b_global[50]  << "\n";
        //cout <<" b at random node " << b[550]<< " b_global: " << b_global[550]  <<"\n";
        //cout <<" b at random node " << b[1230]<< " b_global: " << b_global[1230]  <<"\n";
    }
    else
    {

        DGMeshPartition.ReceiveBottom(MPI,b);

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
        EdgeReversed[is] = EdgeData[is*8+4]; 		// save just the orientation
    }
    for(int ie=1; ie<=Nelem; ++ie)
    {
        int id = (ie-1)*4;
        int idFace;
        for (int is=0; is<4; is++)
        {
            int ifa  = DGMeshPartition.MyElementToEdge[id+is];
            idFace = (ifa-1)*8;


            ElemEdgeMasterSlave[id+is]=DGMeshPartition.MyElemEdgeMasterSlave[id+is];
            if (EdgeData[idFace]==ie-1)
            {
                ElemEdgeOrientation[id+is]=1;
            }
            else
            {
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

//    for(int ie=0; ie<Nfaces; ++ie)
//    {
//        if (EdgeData[ie*8 + 4] == 0)
//        {
//            cout << "\n face "<< ie << " is rotated for right element !\n  the normal is: \n" ;
//        }
//
//        for(int i=0; i<ngl; ++i)
//        {
//            int id = ie*ngl  +i ;
//            cout << normalsX[id] << " ";
//            cout << normalsY[id] << " ";
//            dfloat  normalNorm = pow(normalsX[id],2) + pow(normalsY[id],2);
//            cout  << " NORMALNORM: " << sqrt(normalNorm) <<" \n";
//
//        }
//    }


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
    cout << "RANK: " << MPI.rank << " local min ele Size: "<< minEleSize << " global min Ele Size: " << globalMinEleSize<< "\n";




    dfloat Dmat[ngl2];
    dfloat Dmat0[ngl2];
    dfloat Dhat[ngl2];
    dfloat VdmInv[ngl2];

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

            //            SubCellMat[Did] = DGBasis.SubCellMat(i+1,l+1);
        }
        GLw[i] = DGBasis.w_GL[i];
    }

//    cout <<"\n D central: \n";
//    for(int j=0; j<ngl; ++j)
//    {
//        for(int i=0; i<ngl; ++i)
//        {
//            int id =   j*ngl+i;
//            cout <<DGBasis.DCentralFD[id]<<"  ";
//        }
//       cout <<"\n";
//    }
//
//    cout <<"\n D forward: \n";
//    for(int j=0; j<ngl; ++j)
//    {
//        for(int i=0; i<ngl; ++i)
//        {
//            int id =   j*ngl+i;
//            cout <<DGBasis.DforwardFD[id]<<"  ";
//        }
//       cout <<"\n";
//    }
//
//    cout <<"\n D backward: \n";
//    for(int j=0; j<ngl; ++j)
//    {
//        for(int i=0; i<ngl; ++i)
//        {
//            int id =   j*ngl+i;
//            cout <<DGBasis.DbackwardFD[id]<<"  ";
//        }
//        cout <<"\n";
//    }



//MetricIdentities
    dfloat * MetricIdentities1 = (dfloat*) calloc(Nelem*ngl2,sizeof(dfloat));
    dfloat * MetricIdentities2 = (dfloat*) calloc(Nelem*ngl2,sizeof(dfloat));
    dfloat maxMetricIDError=0.0;
    for (int ie=0; ie<Nelem; ie++)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int ele_ij = ie*ngl2 +  j*ngl+i;
                int loc_ij = j*ngl+i;
                for (int l=0; l<ngl; l++)
                {

                    int loc_il = i*ngl+l;
                    int loc_jl = j*ngl+l;
                    int ele_il = ie*ngl2 + i*ngl+l;
                    int ele_lj = ie*ngl2 + l*ngl+j;
                    MetricIdentities1[ele_ij] += Dmat0[loc_jl] *y_eta[ele_il]- Dmat0[loc_il] *y_xi[ele_lj];
                    MetricIdentities2[ele_ij] += -Dmat0[loc_jl] *x_eta[ele_il]+ Dmat0[loc_il] *x_xi[ele_lj];
                }
            }
        }
    }

    for (int ie=0; ie<Nelem; ie++)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2+  j*ngl+i;
//                if(abs(MetricIdentities1[id]) >0.000000001)
//                    cout <<"Metric Identities1 for Ele: " << ie <<" "<<MetricIdentities1[id]<<"\n ";
//                if(abs(MetricIdentities2[id]) >0.000000001)
//                    cout <<"Metric Identities2 for Ele: " << ie <<" "<<MetricIdentities2[id]<<"\n  ";
                maxMetricIDError = max(maxMetricIDError,abs(MetricIdentities1[id]));
                maxMetricIDError = max(maxMetricIDError,abs(MetricIdentities2[id]));
            }

        }
    }
    cout << "Maximum Error in Approximation of Metric Identities is: " << maxMetricIDError << " !\n";
    free(MetricIdentities1);
    free(MetricIdentities2);



    //initialise time integrator
    RungeKutta RK(rkorder,rkSSP);

    //INITIALIZE SOLUTION
    // THIS SHOULD NOW BE TAKEN CARE OF ALREADY!!!!!!!!
    // SW_Problem.InitB(1,DGMeshPartition,Nelem,ngl,ngl2,x_phy,y_phy,b);
// THIS SHOULD NOW BE TAKEN CARE OF ALREADY!!!!!!!!



    SW_Problem.CalcBDerivatives(Nelem,ngl,ngl2,g_const,x_phy,y_phy,b,Dmat0,y_eta,y_xi,x_eta,x_xi,Bx,By,J);
//	cout << "i am rank: " << MPI.rank << " and my Nelem_global is " << DGMeshPartition.global_NumElements << "\n";

    SW_Problem.InitQ(1,DGMeshPartition,Nelem,ngl,ngl2,x_phy,y_phy,q,0.0,b, g_const,h_0);

    if(MPI.rank==0)
    {
        cout <<"... finished.\n";

    }

    if ((Testcase==90)||(Testcase==91))
    {
        if(MPI.rank==0)
        {

            cout <<"Now initializing okada displacements...      ";
        }
        okada okadamap(Nelem);

        okadamap.okadamapFull(ngl, x_phy,y_phy, q);

        if(MPI.rank==0)
        {
            cout <<"... finished.\n";

        }
    }

    if(MPI.rank==0)
    {
        cout <<"Allocating Memory on the device...      ";
    }





    if(MPI.rank==0)
    {
        cout <<" ...Memory on Device allocated.\n";
    }




    const dfloat geomface = 1.0/DGBasis.w_GL[0];



    occa_device.initDeviceVariables(N, Nelem,Nfaces,MPI.rank, rkSSP, NEpad,NEsurfpad, Nedgepad,NavgPad, ES, Testcase, epsilon_0, sigma_max, sigma_min, PosPresTOL, geomface, g_const,
                                    PositivityPreserving,
                                    ArtificialViscosity, DiscBottom,h_0,PartialDryTreatment,FrictionTerms);
    //copy all permanent data onto the device
    if(MPI.rank==0)
    {
        cout <<"Copy Necessary Data onto Device...";
    }

    occa_device.copyDeviceVariables(  PositivityPreserving, Nelem,GLw,
                                      normalsX,   normalsY,  Scal,  y_xi, y_eta, x_xi, x_eta, b,  Bx, By,
                                      Dmat, DGBasis.Dstrong, Dhat,  J,  x_phy,  y_phy,  q,  ElementSizes,  gRK,  Qt,
                                      VdmInv, ElemEdgeMasterSlave, ElemEdgeOrientation, ElemToEdge, EdgeData,EdgeReversed );
    if(PartialDryTreatment)
    {
        occa_device.copyPartialDryData(DGBasis.DCentralFD,  DGBasis.DforwardFD,  DGBasis.DbackwardFD   );
    }

    if(MPI.rank==0)
    {
        occa_device.initGlobalVars(Nelem_global,b_global,x_phy_global,y_phy_global,J_global);
    }

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
    free(EdgeReversed );
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
    occa_device.buildDeviceKernels(  KernelVersion,   KernelVersionSTD,   Testcase,   Fluxdifferencing,   NumFlux,  rkSSP,   ArtificialViscosity,   PositivityPreserving, DiscBottom);

    cout <<"rank: " << MPI.rank << " GB allocated on device: " << occa_device.device.bytesAllocated()/(1024.f*1024*1024) << "\n";












    //printing initial conditions




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


    cout << "Minimum Ele Size: " << globalMinEleSize << ".\n";

    occa_device.DGtimeloop(DGMeshPartition.NumElements,
                           Nfaces,
                           MPI,
                           DGMeshPartition,
                           RK,
                           DGBasis,
                           SW_Problem,
                           globalMinEleSize,
                           CFL,
                           DFL,
                           T,
                           Testcase,
                           ArtificialViscosity,
                           rkSSP,
                           NumPlots,
                           NumTimeChecks,
                           g_const,
                           PlotVar,
                           EntropyPlot,
                           PositivityPreserving,
                           DiscBottom);



    // END TIME LOOP

    occa_device.o_q.copyTo(q);



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////ANALYSIS OF SOLUTION ///////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    if(MPI.rank==0)
    {
        gettimeofday(&endTV, NULL);

        double TimeDifferenceMeasured = ((endTV.tv_sec  - startTV.tv_sec) * 1000000u + endTV.tv_usec - startTV.tv_usec) / 1.e6;


        //endTV.tv_sec  - startTV.tv_sec +  (endTV.tv_usec - startTV.tv_usec) / 1000000 ;

        cout <<" MPI Ranks: " << MPI.numtasks<<" time taken in seconds = "<< TimeDifferenceMeasured<< "\n";
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


        SW_Problem.InitQ(0,DGMeshPartition,Nelem_global,ngl,ngl2,x_phy_global,y_phy_global,q_exakt,0,b_global,g_const,h_0);


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

        SW_Problem.InitQ(0,DGMeshPartition,Nelem_global,ngl,ngl2,x_phy_global,y_phy_global,q_exakt,T,b_global,g_const,h_0);



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
            L2Error[ik]= sqrt(fabs(L2Error[ik]));
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











    occa_device.freeOccaVars( rkSSP, PositivityPreserving, ArtificialViscosity );
    if (MPI.rank==0)
    {
        occa_device.freeGlobalVars();
    }

//    o_PackSend.free();
//    o_PackReceive.free();

// done with MPI
    MPI_Finalize();

    return 0;
}






