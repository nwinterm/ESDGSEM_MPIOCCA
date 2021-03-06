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


#ifndef DEVICECLASS_H
#define DEVICECLASS_H
#include "occa.hpp"
#include "Constants.hpp"
#include "MPI_Communication.h"
#include "MPI_setup.h"
#include "MeshPartitioning.h"
#include "plots.h"
#include "RungeKutta.h"
#include "SW2D.h"
#include "basis.h"

class deviceclass
{
public:
    deviceclass(const int, const int, char *argv[]);
    void initDeviceVariables(const int N,
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
                             const dfloat h_0,
                             const int PartialDry,
                             const int FrictionTerms,
                             const int CalcArrivalTimes,
                             const int CreateTimeSeries,
                             const int CalcMaximumElevation);


    void buildDeviceKernels(const int KernelVersion,
                            const int KernelVersionSTD,
                            const int Testcase,
                            const int FluxDifferencing,
                            const int NumFlux,
                            const int rkSSP,
                            const int ArtificialViscosity,
                            const int PositivityPreserving,
                            const int DiscBottom);
    void copyPartialDryData(const dfloat* DCentralFD,
                            const dfloat* DforwardFD,
                            const dfloat* DbackwardFD);
    void copyTimeSeriesIDs(const int id1,
                           const int id2,
                           const int id3,
                           const int id4,
                           const int id5,
                           const int id6,
                           const int id7);
    void copyDeviceVariables(const int PositivityPreserving,
                             const int Nelem,
                             const dfloat* GLw,
                             const dfloat * normalsX,
                             const dfloat * normalsY,
                             const dfloat* Scal,
                             const dfloat* y_xi,
                             const dfloat*y_eta,
                             const dfloat*x_xi,
                             const dfloat*x_eta,
                             const dfloat*b,
                             const dfloat* Bx,
                             const dfloat*By,
                             const dfloat* Dmat,
                             const dfloat*Dstrong,
                             const dfloat*Dhat,
                             const dfloat* J,
                             const dfloat* x_phy,
                             const dfloat* y_phy,
                             const dfloat* q,
                             const dfloat* ElementSizes,
                             const dfloat* gRK,
                             const dfloat* Qt,
                             const dfloat* VdmInv,
                             const int*ElemEdgeMasterSlave,
                             const int*ElemEdgeOrientation,
                             const int*ElemToEdge,
                             const int*EdgeData,
                             const int*EdgeReversed);
    void freeOccaVars(const int rkSSP,
                      const int PositivityPreserving,
                      const int ArtificialViscosity );
    void initGlobalVars(const int Nelem_glb,
                        const dfloat *b_glb,
                        const dfloat* x_glb,
                        const dfloat *y_glb,
                        const dfloat *J_glb);
    void freeGlobalVars();
    void  DGtimeloop(const int NumElements,
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
                     const int DiscBottom);
    virtual ~deviceclass();

    occa::device device;
    occa::kernelInfo info;
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
    occa::kernel SurfaceKernelDiscBottom;

    occa::kernel FrictionSource;


    occa::kernel preservePosivitity;
    occa::kernel calcAvg;
    occa::kernel FindLambdaMax;

    occa::kernel scaleGradient;
    occa::kernel SurfaceKernelVisc;

    occa::kernel FindDryElements;
    occa::kernel VolumeKernelPartialDry;

    occa::kernel setArrivaltimes;
    occa::kernel setMaximumElevation;

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
    occa::memory o_EdgeReversed;

    occa::memory o_qGradientX, o_qGradientY,o_SurfGradientX,o_SurfGradientY ;
    occa::memory o_qGradientXL,o_qGradientXR,o_qGradientYL,o_qGradientYR;
    occa::memory o_SurfacePartsVisc;
    occa::memory o_EleSizes, o_ViscPara;
    occa::memory o_ViscParaL,o_ViscParaR;

    occa::memory o_DBSurf1,o_DBSurf2;
    occa::memory o_LambdaMax;

    occa::memory o_QtVisc;
    occa::memory o_qAvg;

    occa::memory o_GLw;

    occa::memory o_ViscForPlot;

    occa::memory o_DcentralFD, o_DforwardFD, o_DbackwardFD;

    occa::memory o_isPartlyDry;

    occa::memory o_FrictionForPlot;

    occa::memory o_ArrivalTimings;
    occa::memory o_MaximumElevation;


    dfloat * ChennaiTimeSeries;
    dfloat * TuticorinTimeSeries;
    dfloat * VisakhapatnamTimeSeries;
    dfloat * ParadipTimeSeries;
    dfloat * KochiTimeSeries;
    dfloat * MormugaoTimeSeries;
    dfloat * OkhaTimeSeries;
    dfloat * TimeSeriesTimes;


    dfloat * EntropyOverTime;
    dfloat * MassOverTime;
    dfloat * EntropyTimes;
    dfloat * maximumFriction;
    dfloat * ArrivalTimings;
    dfloat * maximumElevation;


protected:

private:


    int ngl;
    int ngl2;
    int Nelem_global;
    dfloat * b_global;
    dfloat * q_global;
    dfloat * x_phy_global;
    dfloat * y_phy_global;
    dfloat * J_global;
    dfloat * ViscPara_Global;
    dfloat h_0;
    int PartialDryTreatment;
    int CalcFrictionTerms;
    int calcArrivalTimes;
    int createTimeSeries;
    int calcMaximumElevation;


    int chennaiID;
    int paradipID;
    int tuticorinID;
    int viskhapatnamID;
    int kochiID;
    int mormugaoID;
    int okhaID;

};

#endif // DEVICE_H
