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
