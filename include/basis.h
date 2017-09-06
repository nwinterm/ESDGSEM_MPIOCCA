#ifndef BASIS_H
#define BASIS_H


#include "Constants.hpp"
#include "SW2D.h"
using namespace Constants;


//fmatrix D;
//fmatrix x_GL;
//fmatrix w_GL;
//fmatrix w_bary;
//
//void initBasis();
//void LGLNodesAndWeights( );
//void qAndLEvaluation( const double x_GL,double *q, double *del_q, double *L_N);
//void BarycentricWeights();
//void PolynomialDerivativeMatrix();




class basis
{
    public:
        basis(const int,const bool);
        int N,ngl,ngl2;//,NelemX,NelemY;
        int Nelem;
        int Nelem_global;
//        int Nfaces;
//        int ngl, ngl2;
//        int NoDofs;
//        int NoSpaceDofs;

        fmatrix D;
        fmatrix D0;
        fmatrix Dhat;
        fmatrix x_GL;
        fmatrix w_GL;
        fmatrix x_Gauss;
        fmatrix w_Gauss;
        fmatrix w_bary;
        fmatrix Vdm;
        fmatrix VdmInv;
        fmatrix SubCellMat;

        fmatrix L_at_Gauss;
//        void L2Norm(const dfloat[NoDofs],const dfloat[NoDofs],dfloat[Neq]);
//        void LinfNorm(const dfloat[NoDofs],const dfloat[NoDofs],dfloat[Neq]);

        void L2Norm(const dfloat[],const dfloat[],const dfloat[],dfloat[]);
        void calcElementSizes(const dfloat [],dfloat [],dfloat *minEleSize);
        void LinfNorm(const dfloat[],const dfloat[],dfloat[]);
        void setNelem(const int,const int);
        void  calcEntropyDelta(const dfloat g_const,const dfloat Q[],const dfloat Q_init[],const dfloat b[],const dfloat J[],dfloat *EntropyDelta);
        void ConvertToModal(const dfloat *Q_nodal, dfloat *Q_modal);
        void EvaluteModalPolynomial(const dfloat Q_modal[], dfloat Q_nodal[]);
        void PosPreservation(const dfloat [],const dfloat [],dfloat []);
        void  CheckWhereItNaNed(const dfloat [],bool*, int*);
        void  CheckWhereItNaNedTimeDeriv(const dfloat [],bool*, int*);
        void  EdgesCheckWhereItNaNed(const int, const dfloat [],bool*, int*);

    private:
        void LGLNodesAndWeights( );
        void GaussNodesAndWeights( );
        void legendrePolynomialAndDerivative(const int N,const dfloat x,dfloat *L_N,dfloat *dL_N);
        void qAndLEvaluation( const dfloat x_GL,dfloat *q, dfloat *del_q, dfloat *L_N);
        void BarycentricWeights();
        void PolynomialDerivativeMatrix();
        void ModalTrafoMatrix();
        void SubCellAverageMatrix();
        void calcEntropyPointwise(const dfloat g_const,const dfloat h,const dfloat hu,const dfloat hv,const dfloat b,dfloat *Entropy);
        void LagrangeInterpolatingPolynomial(const dfloat x0,const int N,const fmatrix x,const fmatrix w,dfloat Polynomial[]);
};

#endif // BASIS_H
