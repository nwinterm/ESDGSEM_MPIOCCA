#ifndef MESH_H
#define MESH_H

#include "Constants.hpp"
#include "ReadInput.h"
#include <stdexcept>
using namespace Constants;



class Mesh
{
    public:
        Mesh(const dfloat*,const int ,const int);
        void InitMesh(const string,const bool, const int );
        int m_num_nodes;
        int m_num_edges;
        int m_num_elements;
        int m_order_of_boundary_edges;

        dfloat *  x_global;
        dfloat *  y_global;
        dfloat *  xXi_global;
        dfloat * xEta_global;
        dfloat *  yXi_global;
        dfloat * yEta_global;
        dfloat *  J_global;
        dfloat *  nx_global;
        dfloat * ny_global;
        dfloat * scal_global;
        dfloat *  NormalsX;
        dfloat * NormalsY;
        dfloat * Scal;



        int * ElementToEdge;
        int * EdgeInfo;

        int * ElemEdgeMasterSlave;
        int * ElemEdgeOrientation;
//        imatrix ElemEdgeMasterSlave;
//        imatrix ElemEdgeOrientation;

        void ReadMesh(const string);
        void GenerateMesh(const dfloat ,const dfloat ,const dfloat ,const dfloat ,const bool ,const bool );
        void InitDomain(const int , int* ,int *, int *,int *, bool *, bool *, dfloat* , dfloat *, dfloat* , dfloat *);
        virtual ~Mesh();




    protected:

    private:
        int ngl,ngl2;
        int NelemX;
        int NelemY;
//        fmatrix x_GL;
        dfloat *  x_GL;
        dfloat * x_cheby;
        dfloat * w_bary;

        dfloat * x_phy;
        dfloat * y_phy;

        dfloat * x_xi;
        dfloat * x_eta;
        dfloat * y_xi;
        dfloat * y_eta;
        dfloat *J;
        dfloat * x_bndy;
        dfloat * y_bndy;
        dfloat * scal;
        dfloat * nx;
        dfloat *ny;

//        fmatrix x_phy,y_phy;
//        fmatrix x_xi,x_eta,y_xi,y_eta,J;
//        fmatrix x_bndy,y_bndy,scal;
//        fmatrix nx,ny;


        void ConstructMappedGeometry(const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const bool);

        void TransfiniteQuadMetrics(const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat,const dfloat,dfloat*,dfloat*,dfloat*,dfloat*);
        void TransfiniteQuadMap(const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat,const dfloat,dfloat*,dfloat*);
        void QuadMapMetrics(const dfloat * ,const dfloat * ,const dfloat,const dfloat,dfloat*,dfloat*,dfloat*,dfloat*);
        void QuadMap(const dfloat * ,const dfloat * ,const dfloat,const dfloat,dfloat*,dfloat*);



        void BarycentricWeights();
        void LagrangeInterpolation(const dfloat,const dfloat*,dfloat*);
        void LagrangeInterpolantDerivative(const dfloat,const dfloat*,dfloat*);
        void DerivativeAt(const dfloat*,const dfloat*,const dfloat,dfloat*,dfloat*);
        void EvaluateAt(const dfloat*,const dfloat*,const dfloat,dfloat*,dfloat*);

};

#endif // MESH_H
