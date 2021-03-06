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


#ifndef BASIS_H
#define BASIS_H


#include "Constants.hpp"
#include "SW2D.h"
using namespace Constants;



class basis
{
public:
    basis(const int,const int);
    int N,ngl,ngl2;//,NelemX,NelemY;
    int Nelem;
    int Nelem_global;
//        int Nfaces;
//        int ngl, ngl2;
//        int NoDofs;
//        int NoSpaceDofs;


    dfloat * D;
    dfloat * D_SBP;
    dfloat * DCentralFD;
    dfloat * DbackwardFD;
    dfloat * DforwardFD;
    dfloat * Dstrong;
    dfloat * D0;
    dfloat *  Dhat;
    dfloat *  x_GL;
    dfloat *  w_GL;
    dfloat *  x_Gauss;
    dfloat *  w_Gauss;
    dfloat *  w_bary;
    dfloat *  Vdm;
    dfloat *  VdmInv;
    dfloat *  SubCellMat;

//        dfloat *  L_at_Gauss;
//        void L2Norm(const dfloat[NoDofs],const dfloat[NoDofs],dfloat[Neq]);
//        void LinfNorm(const dfloat[NoDofs],const dfloat[NoDofs],dfloat[Neq]);

    void L2Norm(const dfloat[],const dfloat[],const dfloat[],dfloat[]);
    void calcElementSizes(const dfloat [],dfloat [],dfloat *minEleSize);
    void UpdateMaximumFriction(const dfloat FrictionForPlot[],dfloat maximumFriction[]);
    void LinfNorm(const dfloat[],const dfloat[],dfloat[]);
    void setNelemLocal(const int);
    void setNelemGlobal(const int);
    void  calcEntropyDelta(const dfloat g_const,const dfloat Q[],const dfloat Q_init[],const dfloat b[],const dfloat J[],dfloat *EntropyDelta, dfloat *relEntropyDelta);
    void calcTotalEntropy(const dfloat g_const,const dfloat Q[],const dfloat b[],const dfloat J[],dfloat *TotalEntropy);
    void calcTotalMass(const dfloat Q[],const dfloat J[],dfloat *TotalMass);
    void checkConservation(const dfloat Q[],const dfloat Q_init[],const dfloat J[],dfloat *MassDelta,dfloat *relMassError);
    void ConvertToModal(const dfloat *Q_nodal, dfloat *Q_modal);
    void EvaluteModalPolynomial(const dfloat Q_modal[], dfloat Q_nodal[]);
//        void  CheckWhereItNaNed(const dfloat [],bool*, int*);
//        void  CheckWhereItNaNedTimeDeriv(const dfloat [],bool*, int*);
//        void  EdgesCheckWhereItNaNed(const int, const dfloat [],bool*, int*);

private:
    void LGLNodesAndWeights( );
    void GaussNodesAndWeights( );
    void legendrePolynomialAndDerivative(const int N,const dfloat x,dfloat *L_N,dfloat *dL_N);
    void qAndLEvaluation( const dfloat x_GL,dfloat *q, dfloat *del_q, dfloat *L_N);
    void BarycentricWeights();
    void PolynomialDerivativeMatrix();
    void FiniteDifferenceOperators();
    void ModalTrafoMatrix();
//        void SubCellAverageMatrix();
    void calcEntropyPointwise(const dfloat g_const,const dfloat h,const dfloat hu,const dfloat hv,const dfloat b,dfloat *Entropy);
    void LagrangeInterpolatingPolynomial(const dfloat x0,const int N,const dfloat*,const dfloat*,dfloat Polynomial[]);
};

#endif // BASIS_H
