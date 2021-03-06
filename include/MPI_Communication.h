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


#ifndef MPI_Communication_H
#define MPI_Communication_H


#include "Constants.hpp"
#include "MPI_setup.h"
#include "MeshPartitioning.h"
using namespace Constants;

void ShareInputData(MPI_setup, int *,dfloat *,dfloat *, dfloat *,dfloat *,int *,int *,dfloat *,dfloat *,dfloat *,
                    dfloat *,int *,int*,int *,int *,int *, int *,int *, int *,int*, int *, int *,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);




void CollectSolution(MPI_setup, const MeshPartitioning,  const dfloat [], dfloat []);
void SendSolution(MPI_setup, const MeshPartitioning, const dfloat []);
void CollectEdgeDataMPI(MPI_setup, const MeshPartitioning, dfloat [], dfloat []);
void CollectEdgeDataMPI_Bonly(MPI_setup, const MeshPartitioning, dfloat [], dfloat []);

void CollectViscoseEdgeDataMPI(MPI_setup, const MeshPartitioning, dfloat [], dfloat [], dfloat [], dfloat [], dfloat [], dfloat []);


void GetGlobalLambdaMax(MPI_setup, const MeshPartitioning,const dfloat [], dfloat *);
void GetGlobalViscParaMax(MPI_setup, const MeshPartitioning, const dfloat [], dfloat * );
void GetGlobalMinEleSize(MPI_setup, const MeshPartitioning,const dfloat, dfloat *);




void CollectViscPara(MPI_setup, const MeshPartitioning, const dfloat[], dfloat[]);
void SendViscPara(MPI_setup, const MeshPartitioning, const dfloat[]);

void CollectViscosity(MPI_setup, const MeshPartitioning,  const dfloat [],const dfloat [], dfloat [], dfloat []);
void SendViscosity(MPI_setup, const MeshPartitioning, const dfloat [], const dfloat []);
#endif // SW1D_H

