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


#ifndef MESHPARTITIONING_H
#define MESHPARTITIONING_H

#include "Mesh.h"
#include "MPI_setup.h"
#include "Constants.hpp"
#include <stdexcept>
#include <cstring>
using namespace Constants;



class MeshPartitioning
{
public:
    MeshPartitioning(const int, const int);
    int global_NumElements;
    int global_NumEdges;
    int NumElements;
    int NumProcessors;
    int NumEdges;

    int ngl,ngl2;

    dfloat * x_global;
    dfloat * y_global;
    dfloat * xXi_global;
    dfloat * xEta_global;
    dfloat * yXi_global;
    dfloat * yEta_global;
    dfloat * J_global;
    dfloat * nx_global;
    dfloat * ny_global;
    dfloat * scal_global;

    int * MyEdgeInfo;
    int * MyElementLocalToGlobal;
    int * MyEdgesLocalToGlobal;



    int * ElementLocalToGlobal;


    int * MyElementToEdge;
    int * MyElemEdgeMasterSlave;

    int * ElementsPerProc;

    // store edge indices with each processor for MPI communication!
    int * ProcIndex;
    int * CommTags;


    virtual ~MeshPartitioning();
    void DivideBottom(const MPI_setup, const int, const dfloat *,dfloat*);
    void DivideMesh(const Mesh,const MPI_setup );

    void ReceiveBottom(const MPI_setup, dfloat *);
    void ReceiveMesh(const MPI_setup );

    void SortMPIEdges(const MPI_setup);

protected:

private:


//        fmatrix x_GL;

    int * EdgesPerProc;
    int * EdgeLocalToGlobal;
//        int * EdgeGlobalToLocal;


//        imatrix globalEdgeInfo;
    void isInArray(const int, const int[],const int, int* );
    int isInArray2(const int, const int [],const int );
    void SplitRealValuesBetweenProcs(const MPI_setup,const dfloat *, dfloat *);
    void SplitEdgeRealValuesBetweenProcs(const MPI_setup,const dfloat *, dfloat *);

};

#endif // MESHPARTITIONING_H
