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


#ifndef MESH_H
#define MESH_H

#include "Constants.hpp"
#include "ReadInput.h"
#include <stdexcept>
using namespace Constants;



class Mesh
{
public:
    Mesh(const dfloat*,const int,const int, const int,const int);
    void InitMesh(const string,const bool, const int);
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

    dfloat *  b_global;

    int * ElementToEdge;
    int * EdgeInfo;

    int * ElemEdgeMasterSlave;
    int * ElemEdgeOrientation;
//        imatrix ElemEdgeMasterSlave;
//        imatrix ElemEdgeOrientation;

    void ReadMesh(const string);
    void GenerateMesh(const dfloat,const dfloat,const dfloat,const dfloat,const bool,const bool );
    void InitDomain(const int, int*,int *, int *,int *, bool *, bool *, dfloat*, dfloat *, dfloat*, dfloat *);
    virtual ~Mesh();




protected:

private:
    int ngl,ngl2;
    int NelemX;
    int NelemY;
    int ConvertToKM;
    int ReadBottom; // This variable is needed to indicate that Ocean Bottom Data must be read in from the mesh file!
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

    dfloat * b_phy;
    dfloat * Gamma1b;
    dfloat * Gamma2b;
    dfloat * Gamma3b;
    dfloat * Gamma4b;
    dfloat * cornersb;
    dfloat b_initial;
    dfloat b_final;
    dfloat interpolated_b;

    dfloat * D;
    dfloat * w_baryLGL;
//        fmatrix x_phy,y_phy;
//        fmatrix x_xi,x_eta,y_xi,y_eta,J;
//        fmatrix x_bndy,y_bndy,scal;
//        fmatrix nx,ny;


    void ConstructMappedGeometry(const dfloat*,const dfloat*,const dfloat*,
                                 const dfloat*,const dfloat*,const dfloat*,const dfloat*,
                                 const dfloat*,const dfloat*,const dfloat*,const dfloat*,
                                 const dfloat*,const dfloat*,const dfloat*,const dfloat*,
                                 const bool);

    void TransfiniteQuadMetrics(const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat,const dfloat,dfloat*,dfloat*,dfloat*,dfloat*);
    void TransfiniteQuadMap(const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat,const dfloat,dfloat*,dfloat*);
    void QuadMapMetrics(const dfloat *,const dfloat *,const dfloat,const dfloat,dfloat*,dfloat*,dfloat*,dfloat*);
    void QuadMap(const dfloat *,const dfloat *,const dfloat,const dfloat,dfloat*,dfloat*);



    void BarycentricWeights();
    void LagrangeInterpolation(const dfloat,const dfloat*,dfloat*);
    void LagrangeInterpolantDerivative(const dfloat,const dfloat*,dfloat*);
    void DerivativeAt(const dfloat*,const dfloat*,const dfloat,dfloat*,dfloat*);
    void EvaluateAt(const dfloat*,const dfloat*,const dfloat,dfloat*,dfloat*);
    void DerivativeAtSingle(const dfloat*,const dfloat,dfloat*);
    void EvaluateAtSingle(const dfloat*,const dfloat,dfloat*);
    void TransfiniteQuadMapSingle(const dfloat*,const dfloat*,const dfloat*,const dfloat*,const dfloat,const dfloat,dfloat*);
    void QuadMapSingle(const dfloat *,const dfloat,const dfloat,dfloat*);


    void BarycentricWeightsLGL();
    void PolynomialDerivativeMatrixLGL();

};

#endif // MESH_H
