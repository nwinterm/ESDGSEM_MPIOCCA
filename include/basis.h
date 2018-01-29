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
		dfloat * DCentralFD;
		dfloat * DupwindFD;
		dfloat * DdownwindFD;
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
