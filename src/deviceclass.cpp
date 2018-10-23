#include "deviceclass.h"
//#include <limits>
//using std::numeric_limits;

deviceclass::deviceclass(const int rank, const int argc, char *argv[])
{
    if (rank == 0)
    {
        occa::printAvailableDevices();
    }




    if(rank==0)
    {
        cout << "Initializing device...";
    }


    if(!strcmp(argv[1], "Serial"))
        device.setup("mode = Serial");
    else if(!strcmp(argv[1], "OpenMP"))
        device.setup("mode = OpenMP  , schedule = compact, chunk = 10");
    else if(!strcmp(argv[1], "OpenCL"))
        if(rank==0)  //HACK WAY FOR 2 GPUs
        {
            device.setup("mode = OpenCL  , platformID = 0, deviceID = 0");
        }
        else
        {
            device.setup("mode = OpenCL  , platformID = 0, deviceID = 1");
        }
    else if(!strcmp(argv[1], "CUDA"))
        if(rank==0)     //HACK WAY FOR 2 GPUs
        {
            device.setup("mode = CUDA    , deviceID = 0");
        }
        else
        {
            device.setup("mode = CUDA    , deviceID = 1");
        }






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




    if(rank==0)
    {
        cout <<"declaring occa Kernels and Variables... ";
    }







}
void deviceclass:: initDeviceVariables(const int N,
                                       const int Nelem,
                                       const int Nfaces,
                                       const int rank,
                                       const int rkSSP,
                                       const int NEpad,
                                       const int NEsurfpad,
                                       const int Nedgepad,
                                       const int NAvgPad,
                                       const int ES,
                                       const int Testcase,
                                       const dfloat epsilon_0,
                                       const dfloat sigma_max,
                                       const dfloat sigma_min,
                                       const dfloat PosPresTOL,
                                       const dfloat geomface,
                                       const dfloat g_const,
                                       const int PositivityPreserving,
                                       const int ArtificialViscosity,
                                       const int DiscBottom,
                                       const dfloat h_0_in,
                                       const int PartialDry,
                                       const int FrictionTerms,
                                       const int CalcArrivalTimes,
                                       const int CreateTimeSeries)
{
    const int Neq=3;
    const int GradNeq = Neq-1;
    ngl =N+1;
    ngl2=ngl*ngl;
    int nglPad=0;
    if (ngl%4==0)
    {
        nglPad=1;
    }
    h_0=h_0_in;
    calcArrivalTimes=CalcArrivalTimes;
    createTimeSeries = CreateTimeSeries;
    PartialDryTreatment=PartialDry;
    CalcFrictionTerms=FrictionTerms;
    dfloat TOL_PosPres = PosPresTOL;//pow(10.0,-4);
//    dfloat ZeroTOL = pow(10.0,-12);	// double precision
    dfloat ZeroTOL = pow(10.0,-5);
//  	cout << "The range for type float is from "
//       << numeric_limits<dfloat>::min()
//       << " to "
//       << numeric_limits<dfloat>::denorm_min()
//	<< " \n" ;

    const dfloat zero = 0.0;
    const dfloat half = 0.5;
    const dfloat one = 1.0;
    const dfloat fourth = 0.25;
    const dfloat eight = 8.0;
    const dfloat two = 2.0;
    const dfloat onepointfive = 1.5;
    const dfloat halfg = half*g_const;
    const dfloat fourthg = fourth*g_const;
    info.addDefine("procID",rank);
    info.addDefine("NEpad",NEpad);
    info.addDefine("NAvgPad",NAvgPad);
    info.addDefine("Nedgepad",Nedgepad);
    info.addDefine("NEsurfpad",NEsurfpad);
    info.addDefine("ngl",ngl);
    info.addDefine("ngl2",ngl*ngl);
    info.addDefine("Neq",Neq);
    info.addDefine("GradNeq",GradNeq);
    info.addDefine("PI",Constants::PI);
    //info.addDefine("Nelem",Nelem);
    info.addDefine("dfloat", dfloatString);
    info.addDefine("dfloat4", dfloat4String);
    info.addDefine("ES", ES);
    info.addDefine("Testcase", Testcase);

    info.addDefine("eps0",epsilon_0);
    info.addDefine("sigmaMax",sigma_max);
    info.addDefine("sigmaMin",sigma_min);

    info.addDefine("nglPad",nglPad );

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
    cout << "Adding h_0 to kernel as: " << h_0 <<"\n";
    info.addDefine("h_zero",h_0);

    /// const dfloat ManningCoefficient = 0.025;    ///in s/m^(1/3)
    const dfloat ManningCoefficient = 0.025/600.0;    ///in min/km^(1/3)
    const dfloat ManningCoefficient2 = pow(ManningCoefficient,2);
    info.addDefine("ManningCoeff",ManningCoefficient2);
    const dfloat seventhirds = 7.0/3.0;
    info.addDefine("seventhirds",seventhirds);
    ///const dfloat earth_radius= 6.378e6; 	// earth radius in meters
    const dfloat earth_radius= 6.378e3; 	// earth radius in kilometers
    info.addDefine("earth_radius",earth_radius);
    const dfloat OneEightyOverPI = 180.0/M_PI;
    info.addDefine("OneEightyOverPI",OneEightyOverPI);
    /// const dfloat w_angular = 2.0*M_PI/(24.0*3600.0);    //in radians/second
    const dfloat w_angular = 2.0*M_PI/(24.0*60.0);    //in radians/minute
    cout << "Earth's angular velocity is "<< w_angular << "\n";
    info.addDefine("w_angular",w_angular);

    cout << "eps0, sigmaMax, sigmaMin, PosPresTOL "  << epsilon_0 << " " <<  sigma_max << "  "<< sigma_min <<  " " <<TOL_PosPres << "\n";


    const int NoDofs=ngl2*Nelem*Neq;
    const int NoSpaceDofs=ngl2*Nelem;
    const int NoDofsGradient=ngl2*Nelem*(Neq-1);
    const int NoDofsGradientSurf = Nfaces*ngl*(Neq-1);

    if (rkSSP)
    {
        o_Qtmp = device.malloc(NoDofs*sizeof(dfloat));
    }




    o_D  = device.malloc(ngl2*sizeof(dfloat));
    o_Dstrong  = device.malloc(ngl2*sizeof(dfloat));
    o_Dhat  = device.malloc(ngl2*sizeof(dfloat));
    o_VdmInv  = device.malloc(ngl2*sizeof(dfloat));

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
    o_EdgeReversed = device.malloc(Nfaces*sizeof(int));

    o_DBSurf1 = device.malloc(Nfaces*ngl*sizeof(dfloat));
    o_DBSurf2 = device.malloc(Nfaces*ngl*sizeof(dfloat));
    //viscose term

    o_LambdaMax = device.malloc(Nelem*sizeof(dfloat));
    //pospres
    o_EleSizes = device.malloc(Nelem*sizeof(dfloat));

    if (PartialDryTreatment)
    {
        o_DcentralFD  = device.malloc(ngl2*sizeof(dfloat));
        o_DforwardFD  = device.malloc(ngl2*sizeof(dfloat));
        o_DbackwardFD  = device.malloc(ngl2*sizeof(dfloat));
        o_isPartlyDry= device.malloc(Nelem*sizeof(int));
    }

    if (PositivityPreserving == 1)
    {
        o_GLw  = device.malloc(ngl*sizeof(dfloat));
        //o_qAvg = device.malloc(Nelem*5*sizeof(dfloat));

    }
    if (ArtificialViscosity == 1)
    {

        o_qGradientX = device.malloc(NoDofsGradient*sizeof(dfloat));
        o_qGradientY = device.malloc(NoDofsGradient*sizeof(dfloat));

        o_qGradientXL = device.malloc(NoDofsGradientSurf *sizeof(dfloat));
        o_qGradientXR = device.malloc(NoDofsGradientSurf *sizeof(dfloat));
        o_qGradientYL = device.malloc(NoDofsGradientSurf *sizeof(dfloat));
        o_qGradientYR = device.malloc(NoDofsGradientSurf *sizeof(dfloat));
        o_SurfGradientX = device.malloc(NoDofsGradientSurf *sizeof(dfloat));
        o_SurfGradientY = device.malloc(NoDofsGradientSurf *sizeof(dfloat));
        o_SurfacePartsVisc = device.malloc(NoDofsGradientSurf *sizeof(dfloat));
        o_QtVisc = device.malloc(NoDofsGradient*sizeof(dfloat));


        o_ViscPara = device.malloc(Nelem*sizeof(dfloat));
        o_ViscForPlot = device.malloc(Nelem*sizeof(dfloat));

        o_ViscParaL = device.malloc(Nfaces*sizeof(dfloat));
        o_ViscParaR = device.malloc(Nfaces*sizeof(dfloat));

    }

    if(CalcFrictionTerms)
    {

        o_FrictionForPlot= device.malloc(ngl2*Nelem*(Neq-1)*sizeof(dfloat));
    }

    if(calcArrivalTimes)
    {
        dfloat * ArrivalTimings = (dfloat*) calloc(ngl2*Nelem,sizeof(dfloat));
        for (int i =0; i<Nelem*ngl2; i++)
        {
            ArrivalTimings[i]=-1.0;
        }
        o_ArrivalTimings= device.malloc(ngl2*Nelem*sizeof(dfloat));
        o_ArrivalTimings.copyFrom(ArrivalTimings);

    }



}


void deviceclass:: buildDeviceKernels(const int KernelVersion,
                                      const int KernelVersionSTD,
                                      const int Testcase,
                                      const int FluxDifferencing,
                                      const int NumFlux,
                                      const int rkSSP,
                                      const int ArtificialViscosity,
                                      const int PositivityPreserving,
                                      const int DiscBottom )
{


    switch(Testcase)
    {
    case 1:  // THIS INCLUDES DIRICHLET BOUNDARIES FOR PERIODIC CONVERGENCE TEST
    {
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/Dirichlet_ConvTest.okl","CollectEdgeData",info);
        addS = device.buildKernelFromSource("okl/ManufacturedSolutions/S_ConvTest.okl","addS",info);
        break;
    }
    case 8:  // THIS INCLUDES DIRICHLET BOUNDARIES FOR PERIODIC CONVERGENCE TEST
    {
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/Dirichlet_ConvTest3.okl","CollectEdgeData",info);
        addS = device.buildKernelFromSource("okl/ManufacturedSolutions/S_ConvTest3.okl","addS",info);
        break;
    }
    case 88:  // THIS INCLUDES DIRICHLET BOUNDARIES FOR PERIODIC CONVERGENCE TEST
    {
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/Dirichlet_ConvTest4.okl","CollectEdgeData",info);
        addS = device.buildKernelFromSource("okl/ManufacturedSolutions/S_ConvTest4.okl","addS",info);
        break;
    }
    case 89:  // THIS INCLUDES DIRICHLET BOUNDARIES FOR PERIODIC CONVERGENCE TEST
    {
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/Dirichlet_ConvTest5.okl","CollectEdgeData",info);
        addS = device.buildKernelFromSource("okl/ManufacturedSolutions/S_ConvTest5.okl","addS",info);
        break;
    }
    case 6:  // THIS INCLUDES DIRICHLET BOUNDARIES FOR PERIODIC CONVERGENCE TEST
    {
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/InnerOuter.okl","CollectEdgeData",info);
        break;
    }
    case 32:  // Inflow Boundaries for 3 Mound PP test case
    {
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/3MoundInflow.okl","CollectEdgeData",info);
        break;
    }
    case 90:  // THIS INCLUDES DIRICHLET BOUNDARIES FOR PERIODIC CONVERGENCE TEST
    {
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/OceanBoundary.okl","CollectEdgeData",info);
        break;
    }
    case 91:  // THIS INCLUDES DIRICHLET BOUNDARIES FOR PERIODIC CONVERGENCE TEST
    {
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/OceanBoundary.okl","CollectEdgeData",info);
        break;
    }
    default:
        CollectEdgeData=device.buildKernelFromSource("okl/GatherEdgeData/SolidWalls.okl","CollectEdgeData",info);
        break;
    }
    //cout <<"rank: " << MPI.rank << " I AM AT: Collecting Kernels Post TESTCASE" << "\n";
    CollectEdgeData_Bottom=device.buildKernelFromSource("okl/GatherEdgeData/CollectEdgeData_Bottom.okl","CollectEdgeData_Bottom",info);



    if(FluxDifferencing==1)
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

    // Wet Dry Finite Difference Volume kernel for partially dry elements
    if (PartialDryTreatment)
    {
        FindDryElements=   device.buildKernelFromSource("okl/PartialDryTreatment/FindDryElements.okl","FindDryElements",info);
        VolumeKernelPartialDry=device.buildKernelFromSource("okl/PartialDryTreatment/VolumeKernelFD.okl","VolumeKernelPartialDry",info);
    }

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
    case 4:
    {
        calcNumFluxes=device.buildKernelFromSource("okl/RiemannSolvers/calcRoeFlux.okl","calcNumFluxes",info);
        break;
    }
    case 5:
    {
        calcNumFluxes=device.buildKernelFromSource("okl/RiemannSolvers/calcCombinedFluxOcean.okl","calcNumFluxes",info);
        break;
    }
    }

    // calc specific value on edges like jump in b and average h
//    calcEdgeValues          =   device.buildKernelFromSource("okl/DiscontinuousBathimetry/calcEdgeValues.okl","calcEdgeValues",info);

    if (DiscBottom==1)
    {
        // adds additional surface terms due to a possibly discontinuous bottom topography
        calcDiscBottomSurf      =   device.buildKernelFromSource("okl/DiscontinuousBathimetry/calcDiscBottomSurf.okl","calcDiscBottomSurf",info);
        SurfaceKernelDiscBottom       =   device.buildKernelFromSource("okl/DiscontinuousBathimetry/SurfaceKernelDiscBottom.okl","SurfaceKernelDiscBottom",info);
    }
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

    if(CalcFrictionTerms)
    {

        FrictionSource=device.buildKernelFromSource("okl/FrictionSourceTerm/FrictionSource.okl","FrictionSource",info);
    }
    if(calcArrivalTimes)
    {
        setArrivaltimes =   device.buildKernelFromSource("okl/Analysis/ArrivalTimes.okl","ArrivalTimes",info);

    }

    if (ArtificialViscosity)
    {

        CollectEdgeDataGradient =   device.buildKernelFromSource("okl/BR1_Gradient/CollectEdgeDataGradient.okl","CollectEdgeDataGradient",info);
        SurfaceKernelGradient   =   device.buildKernelFromSource("okl/BR1_Gradient/SurfaceKernelGradient.okl","SurfaceKernelGradient",info);
        calcGradient            =   device.buildKernelFromSource("okl/BR1_Gradient/calcGradient.okl","calcGradient",info);
        calcNumFluxesGradient   =   device.buildKernelFromSource("okl/BR1_Gradient/calcNumFluxesGradient.okl","calcNumFluxesGradient",info);
        scaleGradient          =   device.buildKernelFromSource("okl/BR1_Gradient/scaleGradient.okl","scaleGradient",info);
        VolumeKernelViscose     =   device.buildKernelFromSource("okl/ViscoseParts/VolumeKernelViscose.okl","VolumeKernelViscose",info);
        calcNumFluxesViscose    =   device.buildKernelFromSource("okl/ViscoseParts/calcNumFluxesViscose.okl","calcNumFluxesViscose",info);
        SurfaceKernelVisc       =   device.buildKernelFromSource("okl/ViscoseParts/SurfaceKernelVisc.okl","SurfaceKernelVisc",info);
        ShockCapturing          =   device.buildKernelFromSource("okl/ShockCapturing/ShockCapturing.okl","ShockCapturing",info);
    }
    if (PositivityPreserving)
    {
        // calcAvg      		=   device.buildKernelFromSource("okl/Positivity/CalcAvgAndMin.okl","calcAvg",info);			// THIS KERNEL IS JUST TOO SLOW. DONT KNOW HOW TO SPEED IT UP
        preservePosivitity      =   device.buildKernelFromSource("okl/Positivity/PosPres.okl","PosPres",info);
    }


//    MemCopyKernel           =   device.buildKernelFromSource("okl/DG/MemCopyComparison.okl","MemCopyComparison",info);




}

void deviceclass:: copyPartialDryData(const dfloat* DCentralFD, const dfloat* DforwardFD, const dfloat* DbackwardFD)
{

    o_DcentralFD.copyFrom(DCentralFD);
    o_DforwardFD.copyFrom(DforwardFD);
    o_DbackwardFD.copyFrom(DbackwardFD);

}


void deviceclass:: copyTimeSeriesIDs(const int id1,
                                     const int id2,
                                     const int id3,
                                     const int id4)
{

    chennaiID =id1;
    paradipID=id2;
    tuticorinID=id3;
    viskhapatnamID=id4;

}


void deviceclass:: copyDeviceVariables( const int PositivityPreserving, const int Nelem,const dfloat* GLw,
                                        const dfloat * normalsX, const dfloat * normalsY, const dfloat* Scal, const dfloat* y_xi, const dfloat*y_eta, const dfloat*x_xi, const dfloat*x_eta, const dfloat*b, const dfloat* Bx, const dfloat*By,
                                        const dfloat* Dmat, const dfloat*Dstrong, const dfloat*Dhat, const dfloat* J, const dfloat* x_phy, const dfloat* y_phy, const dfloat* q, const dfloat* ElementSizes, const dfloat* gRK, const dfloat* Qt,
                                        const dfloat* VdmInv, const int*ElemEdgeMasterSlave, const int*ElemEdgeOrientation, const int*ElemToEdge, const int*EdgeData, const int*EdgeReversed)
{

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

    o_Dstrong.copyFrom(Dstrong);
    o_Dhat.copyFrom(Dhat);
    o_Jac.copyFrom(J);

    o_x.copyFrom(x_phy);
    o_y.copyFrom(y_phy);
    o_q.copyFrom(q);
    o_EleSizes.copyFrom(ElementSizes);
    o_gRK.copyFrom(gRK);
    o_Qt.copyFrom(Qt);
    o_VdmInv.copyFrom(VdmInv);





    dfloat * qavgtmp = (dfloat*) calloc(Nelem*5,sizeof(dfloat));
    if(PositivityPreserving == 1)
    {
        o_GLw.copyFrom(GLw);

        //o_qAvg.copyFrom(qavgtmp);

    }
    free(qavgtmp);

    o_ElemEdgeMasterSlave.copyFrom(ElemEdgeMasterSlave);
    o_ElemEdgeOrientation.copyFrom(ElemEdgeOrientation);
    o_ElemToEdge.copyFrom(ElemToEdge);
    o_EdgeData.copyFrom(EdgeData);
    o_EdgeReversed.copyFrom(EdgeReversed);

}



void deviceclass:: DGtimeloop(const int Nelem,
                              const int Nfaces,
                              MPI_setup MPI,
                              const MeshPartitioning MeshSplit,
                              RungeKutta RK,
                              basis DGBasis,
                              SW2D SW_Problem,
                              const dfloat globalMinEleSize,
                              const dfloat CFL,
                              const dfloat DFL,
                              const dfloat T,
                              const int Testcase,
                              const int ArtificialViscosity,
                              const int rkSSP,
                              int NumPlots,
                              const int NumTimeChecks,
                              const dfloat g_const,
                              const int PlotVar,
                              const int EntropyPlot,
                              const int PositivityPreserving,
                              const int DiscBottom )
{


    // THESE ARE NEEDED THROUGHOUT RUNTIME FOR MPI COMMUNICATION
    dfloat * q = (dfloat*) calloc(ngl2*Nelem*Neq,sizeof(dfloat));
    dfloat * qL = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));
    dfloat * qR = (dfloat*) calloc(Nfaces*ngl*Neq,sizeof(dfloat));
    dfloat * bL = (dfloat*) calloc(Nfaces*ngl,sizeof(dfloat));
    dfloat * bR = (dfloat*) calloc(Nfaces*ngl,sizeof(dfloat));
    dfloat * ViscParaL = (dfloat*) calloc(Nfaces,sizeof(dfloat));
    dfloat * ViscParaR = (dfloat*) calloc(Nfaces,sizeof(dfloat));

    const int gradientSurfaceDimension = Nfaces*ngl*(Neq-1);
    dfloat * qGradientXL = (dfloat*) calloc(gradientSurfaceDimension,sizeof(dfloat));
    dfloat * qGradientXR = (dfloat*) calloc(gradientSurfaceDimension,sizeof(dfloat));
    dfloat * qGradientYL = (dfloat*) calloc(gradientSurfaceDimension,sizeof(dfloat));
    dfloat * qGradientYR = (dfloat*) calloc(gradientSurfaceDimension,sizeof(dfloat));

//DEBUG VARIABLES
    dfloat * Qt = (dfloat*) calloc(ngl2*Nelem*Neq,sizeof(dfloat));
    dfloat * SurfParts = (dfloat*) calloc(ngl*Nfaces*Neq,sizeof(dfloat));
    dfloat * scal = (dfloat*) calloc(ngl*Nfaces,sizeof(dfloat));
    dfloat * nx = (dfloat*) calloc(ngl*Nfaces,sizeof(dfloat));
    dfloat * ny = (dfloat*) calloc(ngl*Nfaces,sizeof(dfloat));


    //  dfloat * Qx = (dfloat*) calloc(ngl2*Nelem*Neq,sizeof(dfloat));
    //  dfloat * Qy = (dfloat*) calloc(ngl2*Nelem*Neq,sizeof(dfloat));

    //dfloat * QtVisc = (dfloat*) calloc(ngl2*Nelem*Neq,sizeof(dfloat));
    if(MPI.rank==0)
    {
        if (EntropyPlot)
        {
            EntropyOverTime = (dfloat*) calloc(NumPlots,sizeof(dfloat));
            MassOverTime = (dfloat*) calloc(NumPlots,sizeof(dfloat));
            EntropyTimes = (dfloat*) calloc(NumPlots,sizeof(dfloat));
        }

        if (createTimeSeries)
        {
            ChennaiTimeSeries = (dfloat*) calloc(NumPlots,sizeof(dfloat));
            TuticorinTimeSeries = (dfloat*) calloc(NumPlots,sizeof(dfloat));
            VisakhapatnamTimeSeries = (dfloat*) calloc(NumPlots,sizeof(dfloat));
            ParadipTimeSeries = (dfloat*) calloc(NumPlots,sizeof(dfloat));
            TimeSeriesTimes = (dfloat*) calloc(NumPlots,sizeof(dfloat));
        }

    }
    dfloat globalLambdaMax=0.0;
    dfloat maxViscPara=0.0;
    dfloat dt_i=0.0;
    dfloat dt_v=0.0;
    dfloat LocalLambdas[Nelem];
    dfloat rkD[]= {0.0,1.0,0.5};
    dfloat t=0.0;
    dfloat dt=0.0;

    dfloat * ViscPara = (dfloat*) calloc(Nelem,sizeof(dfloat));

    dfloat * mCheckpoints;

    if (NumPlots>0)
    {
        if (Testcase == 32)
        {
            NumPlots=5;
        }
        if (Testcase == 31)
        {
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

        if (Testcase == 32)
        {
            mCheckpoints[1]=8.0;
            mCheckpoints[2]=30.0;
            mCheckpoints[3]=300.0;
            mCheckpoints[4]=900.0;
        }
        if (Testcase == 31)
        {
            mCheckpoints[1]=T/12.0;
            mCheckpoints[2]=T/6.0;
            mCheckpoints[3]=T/4.0;
            mCheckpoints[4]=T;
        }
        if ((Testcase == 90)||(Testcase == 91)||(Testcase == 86)||(Testcase == 87))
        {
            if ((NumPlots==5)&& (T==9600))
            {
                mCheckpoints[1]=1200.0;
                mCheckpoints[2]=2400.0;
                mCheckpoints[3]=4800.0;
                mCheckpoints[4]=9600.0;

            }
            if ((NumPlots==5)&& (T==160))
            {
                mCheckpoints[1]=20.0;
                mCheckpoints[2]=40.0;
                mCheckpoints[3]=80.0;
                mCheckpoints[4]=160.0;
            }
        }
        if ((Testcase == 33)||(Testcase == 36))
        {
            if ((NumPlots==7)&& (T==50))
            {
                mCheckpoints[1]=5.0;

                mCheckpoints[2]=10.0;

                mCheckpoints[3]=20.0;

                mCheckpoints[4]=30.0;

                mCheckpoints[5]=40.0;

                mCheckpoints[6]=50.0;

            }




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



    int plotCount=0;
    if(plotCount<NumPlots)
    {
        if(t>=mCheckpoints[plotCount])
        {
            o_q.copyTo(q);

            if(MPI.rank==0)
            {

                CollectSolution( MPI, MeshSplit, q, q_global);

                if (Testcase == 31)
                {
                    dfloat * q_exakt = (dfloat*) calloc(ngl2*Nelem_global*Neq,sizeof(dfloat));
                    SW_Problem.InitQ(0,MeshSplit,Nelem_global,ngl,ngl2,x_phy_global,y_phy_global,q_exakt,t,b_global,g_const,0.0);
                    PlotSolutionWithExact(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,q_global,b_global,plotCount,q_exakt);
                    free(q_exakt);
                }
                else
                {
                    PlotSolution(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,q_global,b_global,plotCount,h_0);
                }

                if (EntropyPlot)
                {
                    dfloat TotalEntropy=0.0;
                    DGBasis.calcTotalEntropy(g_const,q_global,b_global,J_global,&TotalEntropy);
                    EntropyOverTime[plotCount] = TotalEntropy;
                    dfloat TotalMass=0.0;
                    DGBasis.calcTotalMass(q_global,J_global,&TotalMass);
                    MassOverTime[plotCount] = TotalMass;
                    EntropyTimes[plotCount] = t;
                }
                if (createTimeSeries)
                {
                    int id;
                    id = floor(chennaiID/ngl2)*(Neq-1)*ngl2 + chennaiID;
                    ChennaiTimeSeries[plotCount] = q_global[id ]+b_global[chennaiID]-h_0;
                    id = floor(tuticorinID/ngl2)*(Neq-1)*ngl2 + tuticorinID;
                    TuticorinTimeSeries[plotCount]= q_global[(tuticorinID)*Neq+tuticorinID]+b_global[tuticorinID]-h_0;
                    id = floor(viskhapatnamID/ngl2)*(Neq-1)*ngl2 + viskhapatnamID;
                    VisakhapatnamTimeSeries[plotCount]= q_global[(viskhapatnamID)*Neq+viskhapatnamID]+b_global[viskhapatnamID]-h_0;
                    id = floor(paradipID/ngl2)*(Neq-1)*ngl2 + paradipID;
                    ParadipTimeSeries[plotCount]= q_global[(paradipID)*Neq+paradipID]+b_global[paradipID]-h_0;
                    TimeSeriesTimes[plotCount]=t;

                }


            }
            else
            {
                SendSolution(MPI, MeshSplit,  q);
            }

            plotCount=plotCount+1;

        }
    }


    // Exchange information about the bottom topography only ONCE!
    CollectEdgeData_Bottom(Nfaces,o_EdgeData,o_B,o_bL,o_bR);
    if (MeshSplit.NumProcessors>1)
    {
        o_bL.copyTo(bL);
        o_bR.copyTo(bR);
        CollectEdgeDataMPI_Bonly(MPI, MeshSplit, bL,  bR);
        MPI_Waitall(MeshSplit.NumProcessors,MPI.Recv_b_reqs, MPI.stats);
        MPI_Waitall(MeshSplit.NumProcessors,MPI.Send_b_reqs, MPI.stats);
        o_bL.copyFrom(bL);
        o_bR.copyFrom(bR);
    }



    cout << "about to enter time loop\n";
    if (CalcFrictionTerms)
    {
        maximumFriction = (dfloat*) calloc(ngl2*Nelem_global*(Neq-1),sizeof(dfloat));

    }


    while (t<T)
    {
        FindLambdaMax(Nelem, o_q, o_LambdaMax);





        globalLambdaMax=0.0;
        o_LambdaMax.copyTo(LocalLambdas);
        GetGlobalLambdaMax(MPI,  MeshSplit,LocalLambdas, &globalLambdaMax);
        dt_i = globalMinEleSize/(ngl) * CFL /globalLambdaMax;
        if (Testcase==32)
        {
            if (t==0.0)
            {
                dt_i = 0.0001;
            }
        }
        if ( ArtificialViscosity==1)
        {
            ShockCapturing(Nelem, o_q,o_B,o_VdmInv,o_EleSizes,o_ViscPara,o_ViscForPlot);
            o_ViscPara.copyTo(ViscPara);
            GetGlobalViscParaMax(MPI,  MeshSplit,ViscPara, &maxViscPara);

            dt_v = DFL/(pow(ngl,2)) * pow(globalMinEleSize,2) / maxViscPara;
            dt=fmin(T-t,fmin(dt_i,dt_v));
        }
        else
        {
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
                intermediatetime += rkD[rkstage]*dt;
            }
            else
            {
                intermediatetime += rkB*dt;
            }




            CollectEdgeData(Nfaces,o_EdgeData,o_q,o_x,o_y,o_nx,o_ny, o_bL,o_bR, o_qL, o_qR, intermediatetime);
            if (MeshSplit.NumProcessors>1)
            {
                o_qL.copyTo(qL);
                o_qR.copyTo(qR);
                CollectEdgeDataMPI(MPI, MeshSplit, qL, qR);
            }



// CORRECT VOLUME KERNEL
            VolumeKernel(Nelem, o_Jac,o_Yxi,o_Yeta,o_Xxi,o_Xeta,o_q,o_D,o_Bx,o_By,o_Qt);


            o_Qt.copyTo(Qt);
            device.finish();

//            cout <<"\n q_t Post Volume: \n";
//			for (int ie=0;ie<1;ie++){
//					cout <<"Ele: " << ie <<"\n";
//				for(int j=0;j<ngl;++j){
//					for(int i=0;i<ngl;++i){
//						int id = ie*ngl2*Neq +  j*ngl+i;
//					   cout <<Qt[id+ngl2]<<",  ";
//				  }
//					cout <<"\n";
//				}
//			}


            if (PartialDryTreatment==1)
            {
                FindDryElements(Nelem, o_q, o_isPartlyDry);
                VolumeKernelPartialDry(Nelem, o_Jac,o_Yxi,o_Yeta,o_Xxi,o_Xeta,o_q,o_isPartlyDry,o_DcentralFD,o_DforwardFD,o_DbackwardFD,o_B,o_Qt);
            }



            if (MeshSplit.NumProcessors>1)
            {
                MPI_Waitall(MeshSplit.NumProcessors,MPI.Recv_q_reqs, MPI.stats);
                MPI_Waitall(MeshSplit.NumProcessors,MPI.Send_q_reqs, MPI.stats);
                //o_qL.copyFrom(qL);			// qL is never received, this should be unnecessary!
                o_qR.copyFrom(qR);
                //o_qR.copyFrom(qR,1);
            }

            calcNumFluxes(Nfaces,o_EdgeData,o_EdgeReversed,o_nx,o_ny,o_scal,o_qL,o_qR,o_bL,o_bR,o_SurfaceParts);

            SurfaceKernel(Nelem,o_Jac,o_ElemEdgeMasterSlave,o_ElemEdgeOrientation,o_ElemToEdge, o_SurfaceParts,o_DBSurf1,o_DBSurf2,o_Qt);



            o_Qt.copyTo(Qt);
            device.finish();

//            cout <<"\n q_t Post Surface: \n";
//			for (int ie=0;ie<1;ie++){
//					cout <<"Ele: " << ie <<"\n";
//				for(int j=0;j<ngl;++j){
//					for(int i=0;i<ngl;++i){
//						int id = ie*ngl2*Neq +  j*ngl+i;
//					   cout <<Qt[id+ngl2]<<",  ";
//				  }
//					cout <<"\n";
//				}
//			}


            if (DiscBottom==1)
            {
                calcDiscBottomSurf(Nfaces,o_EdgeReversed,o_qL,o_qR, o_bL,o_bR,o_nx,o_ny,o_scal,o_DBSurf1,o_DBSurf2);
                SurfaceKernelDiscBottom(Nelem,o_Jac,o_ElemEdgeMasterSlave,o_ElemEdgeOrientation,o_ElemToEdge, o_DBSurf1,o_DBSurf2,o_Qt);
            }

            if(CalcFrictionTerms)
            {

                FrictionSource(Nelem,o_y,o_q,o_Qt,o_FrictionForPlot);
            }

            if ( ArtificialViscosity==1)
            {


                // at first RK step we already now o_ViscPara from time step computation!
                if (rkstage>0)
                {
                    ShockCapturing(Nelem, o_q,o_B,o_VdmInv,o_EleSizes,o_ViscPara,o_ViscForPlot);
                }

                calcNumFluxesGradient(Nfaces,o_EdgeData,o_nx,o_ny,o_scal, o_qL, o_qR, o_SurfGradientX,o_SurfGradientY);

                calcGradient(Nelem, o_Jac,o_Yxi,o_Yeta,o_Xxi,o_Xeta,o_q,o_B,o_Dhat,o_qGradientX,o_qGradientY);


                SurfaceKernelGradient(Nelem,o_Jac,o_ElemEdgeMasterSlave,o_ElemEdgeOrientation,o_ElemToEdge, o_SurfGradientX, o_SurfGradientY,o_qGradientX,o_qGradientY);

                scaleGradient(Nelem,o_q,o_qGradientX,o_qGradientY);

                CollectEdgeDataGradient(Nfaces,o_EdgeData,o_qGradientX,o_qGradientY,o_ViscPara,o_ViscParaL,o_ViscParaR, o_qGradientXL, o_qGradientXR,o_qGradientYL, o_qGradientYR);

                if (MeshSplit.NumProcessors>1)
                {
                    o_ViscParaL.copyTo(ViscParaL);
                    o_ViscParaR.copyTo(ViscParaR);
                    o_qGradientXL.copyTo(qGradientXL);
                    o_qGradientXR.copyTo(qGradientXR);
                    o_qGradientYL.copyTo(qGradientYL);
                    o_qGradientYR.copyTo(qGradientYR);

                    // not an OCCA kernel, this is MPI communication for exchanging gradient data
                    CollectViscoseEdgeDataMPI(MPI, MeshSplit, ViscParaL, ViscParaR, qGradientXL,  qGradientXR,qGradientYL,qGradientYR);

                }

                VolumeKernelViscose(Nelem, o_Jac,o_Yxi,o_Yeta,o_Xxi,o_Xeta,o_qGradientX,o_qGradientY,o_Dstrong,o_ViscPara,o_Qt);


                if (MeshSplit.NumProcessors>1)
                {
                    MPI_Waitall(MeshSplit.NumProcessors,MPI.Recv_qX_reqs, MPI.stats);
                    MPI_Waitall(MeshSplit.NumProcessors,MPI.Recv_qY_reqs, MPI.stats);
                    MPI_Waitall(MeshSplit.NumProcessors,MPI.Recv_ViscPar_reqs, MPI.stats);
                    MPI_Waitall(MeshSplit.NumProcessors,MPI.Send_qX_reqs, MPI.stats);
                    MPI_Waitall(MeshSplit.NumProcessors,MPI.Send_qY_reqs, MPI.stats);
                    MPI_Waitall(MeshSplit.NumProcessors,MPI.Send_ViscPar_reqs, MPI.stats);
                    o_ViscParaL.copyFrom(ViscParaL);
                    o_ViscParaR.copyFrom(ViscParaR);
                    o_qGradientXL.copyFrom(qGradientXL);
                    o_qGradientXR.copyFrom(qGradientXR);
                    o_qGradientYL.copyFrom(qGradientYL);
                    o_qGradientYR.copyFrom(qGradientYR);
                }

                calcNumFluxesViscose(Nfaces,o_EdgeData,o_nx,o_ny,o_scal,o_ViscParaL,o_ViscParaR,o_qGradientXL,o_qGradientXR,o_qGradientYL,o_qGradientYR,o_SurfacePartsVisc);

                SurfaceKernelVisc(Nelem,o_Jac,o_ElemEdgeMasterSlave,o_ElemEdgeOrientation,o_ElemToEdge, o_SurfacePartsVisc,o_Qt);


            }

            // add manufactured source term for convergence test
            if ((Testcase==1)||(Testcase==8)||(Testcase==88)||(Testcase==89))
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
                //calcAvg(Nelem,o_EleSizes,o_GLw, o_Jac,o_q,o_qAvg);
                //preservePosivitity(Nelem,o_qAvg,o_q);
                preservePosivitity(Nelem,o_EleSizes,o_GLw, o_Jac,o_q);

            }



        }



        t += dt;

        if (calcArrivalTimes)
        {

            setArrivaltimes(Nelem,t,o_q,o_B,o_ArrivalTimings);
        }
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
                    CollectSolution( MPI, MeshSplit, q, q_global);


                    if (Testcase == 31)
                    {
                        dfloat * q_exakt = (dfloat*) calloc(ngl2*Nelem_global*Neq,sizeof(dfloat));
                        SW_Problem.InitQ(0,MeshSplit,Nelem_global,ngl,ngl2,x_phy_global,y_phy_global,q_exakt,t,b_global,g_const,0.0);
                        PlotSolutionWithExact(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,q_global,b_global,plotCount,q_exakt);

                        free(q_exakt);
                    }
                    else
                    {
                        PlotSolution(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,q_global,b_global,plotCount,h_0);
                    }
                    if(ArtificialViscosity==1)
                    {

                        o_ViscForPlot.copyTo(ViscPara);
                        CollectViscPara(MPI,   MeshSplit, ViscPara, ViscPara_Global);
                        PlotViscoseParameter(Nelem_global, ngl, x_phy_global,y_phy_global, ViscPara_Global, plotCount);
//                        CollectViscosity( MPI, MeshSplit, Qx,Qy, Qx_global, Qy_global);
//                        PlotViscosity(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,Qx_global,Qy_global,plotCount);
//                        CollectViscosity( MPI, MeshSplit, QtVisc,Qy, QtVisc_global, Qy_global);
//                        PlotViscosity(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,QtVisc_global,Qy_global,plotCount);

                    }

                    if (EntropyPlot)
                    {
                        dfloat TotalEntropy=0.0;
                        DGBasis.calcTotalEntropy(g_const,q_global,b_global,J_global,&TotalEntropy);
                        EntropyOverTime[plotCount] = TotalEntropy;
                        dfloat TotalMass=0.0;
                        DGBasis.calcTotalMass(q_global,J_global,&TotalMass);
                        MassOverTime[plotCount] = TotalMass;
                        EntropyTimes[plotCount] = t;
                    }
                    if (createTimeSeries)
                    {
                        ChennaiTimeSeries[plotCount] = q_global[(chennaiID-1)*Neq+chennaiID]+b_global[chennaiID]-h_0;
                        TuticorinTimeSeries[plotCount]= q_global[(tuticorinID-1)*Neq+tuticorinID]+b_global[tuticorinID]-h_0;
                        VisakhapatnamTimeSeries[plotCount]= q_global[(viskhapatnamID-1)*Neq+viskhapatnamID]+b_global[viskhapatnamID]-h_0;
                        ParadipTimeSeries[plotCount]= q_global[(paradipID-1)*Neq+paradipID]+b_global[paradipID]-h_0;
                        TimeSeriesTimes[plotCount]=t;

                    }
                    if (CalcFrictionTerms)
                    {
                        if (t<T)
                        {
                            dfloat * FrictionForPlot = (dfloat*) calloc(ngl2*Nelem_global*(Neq-1),sizeof(dfloat));
                            o_FrictionForPlot.copyTo(FrictionForPlot);
                            DGBasis.UpdateMaximumFriction(FrictionForPlot,maximumFriction);
                            free(FrictionForPlot);
                        }
                        else
                        {
                            dfloat * FrictionForPlot = (dfloat*) calloc(ngl2*Nelem_global*(Neq-1),sizeof(dfloat));
                            o_FrictionForPlot.copyTo(FrictionForPlot);
                            DGBasis.UpdateMaximumFriction(FrictionForPlot,maximumFriction);
                            PlotFriction(Nelem_global,ngl,PlotVar,x_phy_global,y_phy_global,maximumFriction,plotCount);
                            free(FrictionForPlot);
                            free(maximumFriction);
                        }

                    }


                }
                else
                {
                    SendSolution(MPI, MeshSplit,  q);
                    if(ArtificialViscosity==1)
                    {
                        o_ViscForPlot.copyTo(ViscPara);
                        SendViscPara(MPI, MeshSplit,  ViscPara);
                        //            SendViscosity(MPI, MeshSplit,  Qx,Qy);
//                        SendViscosity(MPI, MeshSplit,  QtVisc,Qy);
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


    if (MPI.rank==0)
    {
        if (EntropyPlot)
        {
            PlotEntropy(NumPlots, EntropyTimes,EntropyOverTime);
            PlotMass(NumPlots, EntropyTimes,MassOverTime);

        }
        if (createTimeSeries)
        {
            PlotTimeSeries(NumPlots,TimeSeriesTimes,ChennaiTimeSeries,TuticorinTimeSeries,VisakhapatnamTimeSeries,ParadipTimeSeries);

        }
        if (calcArrivalTimes)
        {
            dfloat * ArrivalTimings = (dfloat*) calloc(ngl2*Nelem_global,sizeof(dfloat));
            o_ArrivalTimings.copyTo(ArrivalTimings);
            PlotArrivalTimings(Nelem_global,ngl,x_phy_global,y_phy_global,ArrivalTimings);
            free(o_ArrivalTimings);

        }
    }
    free(q);
    free(ViscPara);
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

}


void deviceclass:: freeOccaVars(const int rkSSP, const int PositivityPreserving, const int ArtificialViscosity )
{

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
    o_EdgeReversed.free();


    o_DBSurf1.free();
    o_DBSurf2.free();


//viscose term

    o_LambdaMax.free();
//pospres
    o_EleSizes.free();


    if (PositivityPreserving == 1)
    {
        o_GLw.free();
        // o_qAvg.free();

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
    }

}



void deviceclass:: initGlobalVars(const int Nelem_glb,
                                  const dfloat *b_glb,
                                  const dfloat* x_glb,
                                  const dfloat *y_glb,
                                  const dfloat *J_glb )
{
    Nelem_global=Nelem_glb;
    int NoDofs_global=ngl2*Nelem_global*Neq;
    int NoSpaceDofs_global=ngl2*Nelem_global;
    q_global = (dfloat*) calloc(NoDofs_global,sizeof(dfloat));
    b_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));
    J_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));
    x_phy_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));
    y_phy_global = (dfloat*) calloc(NoSpaceDofs_global,sizeof(dfloat));
    ViscPara_Global = (dfloat*) calloc(Nelem_global,sizeof(dfloat));
    std::memcpy(b_global, b_glb, NoSpaceDofs_global*sizeof(dfloat));
    std::memcpy(J_global, J_glb, NoSpaceDofs_global*sizeof(dfloat));
    std::memcpy(x_phy_global, x_glb, NoSpaceDofs_global*sizeof(dfloat));
    std::memcpy(y_phy_global, y_glb, NoSpaceDofs_global*sizeof(dfloat));
}


void deviceclass:: freeGlobalVars()
{
    free(b_global);
    free(q_global);
    free(x_phy_global);
    free(y_phy_global);



}

deviceclass::~deviceclass()
{
    //dtor
}
