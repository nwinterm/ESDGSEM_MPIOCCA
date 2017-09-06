#ifndef MESH_H
#define MESH_H

#include "Constants.hpp"
#include "ReadInput.h"
#include <stdexcept>
using namespace Constants;



class Mesh
{
    public:
        Mesh(const fmatrix,const int );
        void InitMesh(const string,const bool, const int );
        int m_num_nodes;
        int m_num_edges;
        int m_num_elements;
        int m_order_of_boundary_edges;

        fmatrix x_global;
        fmatrix y_global;
        fmatrix xXi_global,xEta_global;
        fmatrix yXi_global,yEta_global;
        fmatrix J_global;
        fmatrix nx_global,ny_global,scal_global;
        fmatrix NormalsX, NormalsY, Scal;


        imatrix ElementToEdge;
        imatrix EdgeInfo;
        imatrix ElemEdgeMasterSlave;
        imatrix ElemEdgeOrientation;


        virtual ~Mesh();
        void ReadMesh(const string);
        void GenerateMesh();


    protected:

    private:
        int ngl;
        int NelemX;
        int NelemY;
//        fmatrix x_GL;
        fmatrix x_GL;
        fmatrix x_cheby;
        fmatrix w_bary;
        fmatrix Gamma1,Gamma2,Gamma3,Gamma4;
        fmatrix corners;

        fmatrix x_phy,y_phy;
        fmatrix x_xi,x_eta,y_xi,y_eta,J;
        fmatrix x_bndy,y_bndy,scal;
        fmatrix nx,ny;



        void TransfiniteQuadMetrics(const dfloat,const dfloat,dfloat*,dfloat*,dfloat*,dfloat*);
        void TransfiniteQuadMap(const dfloat,const dfloat,dfloat*,dfloat*);
        void QuadMapMetrics(const dfloat,const dfloat,dfloat*,dfloat*,dfloat*,dfloat*);
        void QuadMap(const dfloat,const dfloat,dfloat*,dfloat*);

        void ConstructMappedGeometry(const bool);


        void BarycentricWeights();
        void LagrangeInterpolation(const dfloat,const dfloat[],dfloat*);
        void LagrangeInterpolantDerivative(const dfloat,const dfloat[],dfloat*);
        void DerivativeAt(const fmatrix,const dfloat,dfloat*,dfloat*);
        void EvaluateAt(const fmatrix,const dfloat,dfloat*,dfloat*);

};

#endif // MESH_H
