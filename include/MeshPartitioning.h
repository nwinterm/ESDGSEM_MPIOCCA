
#ifndef MESHPARTITIONING_H
#define MESHPARTITIONING_H

#include "Mesh.h"
#include "MPI_setup.h"
#include "Constants.hpp"
#include <stdexcept>
using namespace Constants;



class MeshPartitioning
{
    public:
        MeshPartitioning(const int, const int);
        int global_NumElements;
        int global_NumEdges;
        int NumElements;
        int NumProcessors;
        int NumNeighbors;
        int NumEdges;

        int ngl,ngl2;

        fmatrix x_global;
        fmatrix y_global;
        fmatrix xXi_global,xEta_global;
        fmatrix yXi_global,yEta_global;
        fmatrix J_global;
        fmatrix nx_global,ny_global,scal_global;

        imatrix MyElementLocalToGlobal;
        imatrix MyEdgesLocalToGlobal;
        imatrix MyEdgeInfo;

        imatrix ElementGlobalToLocal;
        imatrix ElementLocalToGlobal;
        imatrix EdgeLocalToGlobal;
        imatrix EdgeGlobalToLocal;
        imatrix ElementToEdge;
        imatrix ElemEdgeMasterSlave;

        imatrix ElementsPerProc;

        // store edge indices with each processor for MPI communication!
        imatrix ProcIndex;
        imatrix CommTags;

        // for MPI edge resorting
        imatrix MPIEdges;

        virtual ~MeshPartitioning();

        void DivideMesh(const Mesh ,const MPI_setup );
        void ReceiveMesh(const MPI_setup );
        void ApproximateElementSizes(const int ,const int , const dfloat [], const dfloat [], dfloat []);
        void SortMPIEdges(const MPI_setup);

    protected:

    private:


//        fmatrix x_GL;

        imatrix EdgesPerProc;



//        imatrix globalEdgeInfo;
        void isInArray(const int , const int[]  ,const int , int* );
        void SplitRealValuesBetweenProcs(const MPI_setup,const fmatrix, fmatrix &);
        void SplitEdgeRealValuesBetweenProcs(const MPI_setup ,const fmatrix , fmatrix& );

};

#endif // MESHPARTITIONING_H
