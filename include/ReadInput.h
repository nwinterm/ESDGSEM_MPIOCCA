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


#ifndef READINPUT_H
#define READINPUT_H


#include "Constants.hpp"
using namespace Constants;


void ReadInputFile(int *N,
                   string *meshFile,
                   dfloat *CFL,
                   dfloat *DFL,
                   dfloat *T,
                   dfloat *g_const,
                   int *ArtificialViscosity,
                   int *PositivityPreserving,
                   dfloat *PosPresTOL,
                   dfloat *epsilon_0,
                   dfloat *sigma_min,
                   dfloat *sigma_max,
                   int *PlotVar,
                   int *EntropyPlot,
                   int *NumPlots,
                   int *NumTimeChecks,
                   int *Testcase,
                   int *ES,
                   int *NumFlux,
                   int *FluxDifferencing,
                   int *Cartesian,
                   int * rkorder,
                   int * rkSSP,
                   int * NEpad,
                   int * NEsurfpad,
                   int * Nedgepad,
                   int * NavgPad,
                   int * KernelVersion,
                   int * KernelVersionSTD,
                   int * DiscBottom,
                   int * ReadInBottom,
                   int * PartialDryTreatment,
                   int * FrictionTerms,
                   int * ConvertToKM,
                   int * calcArrivalTimes,
                   int * createTimeSeries,
                   int * calcMaximumElevation);

void ReadCartesianData(const int, const int, dfloat *xL,dfloat *xR,dfloat *yL,dfloat *yR, int *NelemX, int *NelemY,bool *PeriodicBD_X,bool *PeriodicBD_Y);
void WriteFullMesh(const int NumNodes, const dfloat *x,const dfloat *y);
void ReadFullMesh(const int NumNodes,const int Nelem, const int N, dfloat *b, dfloat *h_0);

void FindElementID(const int NumNodes, const dfloat *x,const dfloat *y, const dfloat lonToFind, const dfloat latToFind, int* coordID);
#endif // READINPUT_H
