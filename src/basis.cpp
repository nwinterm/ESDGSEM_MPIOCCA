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



#include "basis.h"
#include "math.h"
#include "cmath"
#include <iostream>



using namespace std;



basis::basis(const int Ninput, const int FluxDifferencing)
{

    N=Ninput;
//NelemX=NelemXinput;
//NelemY=NelemYinput;
//Nelem=NelemX*NelemY;
//Nfaces=(NelemX+1)*NelemY + (NelemY+1)*NelemX; //=2*Nelem + NelemY+NelemX
    ngl=N+1;
    ngl2=ngl*ngl;
//NoDofs=ngl2*Nelem*Neq;
//NoSpaceDofs=ngl2*Nelem;



//const int N=5;

//D.resize(ngl,ngl);
    D = (dfloat*) calloc(ngl2,sizeof(dfloat));
    Dhat = (dfloat*) calloc(ngl2,sizeof(dfloat));
    D0 = (dfloat*) calloc(ngl2,sizeof(dfloat));
    Dstrong = (dfloat*) calloc(ngl2,sizeof(dfloat));





//Dhat.resize(ngl,ngl);
//D0.resize(ngl,ngl);

    x_GL= (dfloat*) calloc(ngl,sizeof(dfloat));
    w_GL= (dfloat*) calloc(ngl,sizeof(dfloat));
    x_Gauss= (dfloat*) calloc(ngl,sizeof(dfloat));
    w_Gauss= (dfloat*) calloc(ngl,sizeof(dfloat));
    w_bary= (dfloat*) calloc(ngl,sizeof(dfloat));
    VdmInv= (dfloat*) calloc(ngl2,sizeof(dfloat));
    Vdm= (dfloat*) calloc(ngl2,sizeof(dfloat));
//L_at_Gauss = (dfloat*) calloc(ngl2,sizeof(dfloat));





    LGLNodesAndWeights();
    GaussNodesAndWeights();

//cout << "LGL Nodes: ";
//for(int i = 0; i < ngl; ++i){
//cout << " " << x_GL[i] << " " ;
//}
//cout << " \n";
//
//cout << "LGL Weights: ";
//for(int i = 0; i < ngl; ++i){
//cout << " " << w_GL[i] << " " ;
//}
//cout << " \n";
//
//cout << "LG Nodes: ";
//for(int i = 0; i < ngl; ++i){
//cout << " " << x_Gauss[i] << " " ;
//}
//cout << " \n";
//
//cout << "LG Weights: ";
//for(int i = 0; i < ngl; ++i){
//cout << " " << w_Gauss[i] << " " ;
//}
//cout << " \n";



    BarycentricWeights();


    PolynomialDerivativeMatrix();



//include surface terms for strong form DG
    for (int i=0; i<ngl; ++i)
    {
        for(int j=0; j<ngl; ++j)
        {
            D0[i*ngl+j]=D[i*ngl+j];
            Dstrong[i*ngl+j]=D[i*ngl+j];
            if (FluxDifferencing)
            {
                D[i*ngl+j]=2.0*D[i*ngl+j];
            }
        }
    }
// Dhat = - M^(-1) * D0^T * M
    for (int i=0; i<ngl; ++i)
    {
        for(int j=0; j<ngl; ++j)
        {
            Dhat[i*ngl+j] = - w_GL[j]*D0[j*ngl+i]/w_GL[i];
        }
    }

//include surface terms for strong form DG
    D[0] = D[0] + 1.0/w_GL[0];
    D[ngl2-1] = D[ngl2-1] - 1.0/w_GL[ngl-1];
    Dstrong[0] = D0[0] + 1.0/w_GL[0];
    Dstrong[ngl2-1] = D0[ngl2-1] - 1.0/w_GL[ngl-1];

    ModalTrafoMatrix();






    DCentralFD = (dfloat*) calloc(ngl2,sizeof(dfloat));

    DforwardFD= (dfloat*) calloc(ngl2,sizeof(dfloat));
    DbackwardFD= (dfloat*) calloc(ngl2,sizeof(dfloat));
    FiniteDifferenceOperators();

     //   D_SBP = (dfloat*) calloc(ngl2,sizeof(dfloat));
//SBP CHECK from david
//for(int i = 0; i < ngl; ++i){
//dfloat Check1=0.0;
//dfloat Check2=0.0;
//        for(int j = 0; j < ngl; ++j){
//		Check1 += D_SBP[i*ngl+j] * 1.0;
//		Check2 += D_SBP[i*ngl+j] * x_GL[j];
//	}
//		cout << "i " << i <<  " Check1: " << Check1 << "\n " ;
//		cout << "i " << i <<  " Check2: " << Check2 << "\n " ;
//}



//cout << "D_SBP: ";
//for(int i = 0; i < ngl; ++i){
//        for(int j = 0; j < ngl; ++j){
//cout << " " << D_SBP[i*ngl+j] << " " ;
//}
//cout << " \n";
//}


//SubCellAverageMatrix();

//dfloat InterpolLagrange[N+1];
//for (int i=0;i<ngl;++i){
//    LagrangeInterpolatingPolynomial(x_Gauss[i],N,x_GL,w_bary,InterpolLagrange);
//    for (int j=0;j<ngl;j++){
//        L_at_Gauss[i*ngl+j] = InterpolLagrange[j];
////        cout << "L_at_Gauss(i,j) " << L_at_Gauss[i*ngl+j] << "\n" ;
//    }
//}


//L_at_Gauss ::  (i,j) = j-tes LagrangePolynom ausgewertet an Gauss Punkt i!!



}




void basis :: FiniteDifferenceOperators()
{

    DCentralFD[0] = (-1.0/(x_GL[1]-x_GL[0]));							// FIRST ENTRY NORMAL UPWIND
    DCentralFD[1] = (1.0/(x_GL[1]-x_GL[0]));							// SECOND ENTRY NORMAL UPWIND
    DCentralFD[ngl2-1] = (1.0/(x_GL[N]-x_GL[N-1]));		//last entry		DOWNWIND
    DCentralFD[ngl2-2] = (-1.0/(x_GL[N]-x_GL[N-1]));		//second last entry		DOWNWIND
    for (int i=1; i<ngl-1; i++)
    {
        const int id = i*ngl+i;
        const int idp1 = id+1;
        const int idm1 = id-1;
        DCentralFD[idp1]=(1.0/(x_GL[i+1]-x_GL[i-1]));
        DCentralFD[idm1]=(-1.0/(x_GL[i+1]-x_GL[i-1]));
    };


    DforwardFD[ngl2-2] = (-1.0/(x_GL[N]-x_GL[N-1]));
    DforwardFD[ngl2-1] = (1.0/(x_GL[N]-x_GL[N-1]));
    for (int i=0; i<ngl-1; i++)
    {
        const int id = i*ngl+i;
        const int idp1 = id+1;
        DforwardFD[id] = (-1.0/(x_GL[i+1]-x_GL[i]));
        DforwardFD[idp1] = (1.0/(x_GL[i+1]-x_GL[i]));
    };

    DbackwardFD[0] = (-1.0/(x_GL[1]-x_GL[0]));
    DbackwardFD[1] = (1.0/(x_GL[1]-x_GL[0]));
    for (int i=1; i<ngl; i++)
    {
        const int id = i*ngl+i;
        const int idm1 = id-1;
        DbackwardFD[id] = (1.0/(x_GL[i]-x_GL[i-1]));
        DbackwardFD[idm1] = (-1.0/(x_GL[i]-x_GL[i-1]));
    };
//    dfloat wFD[ngl];
//    switch(N) {
//    case 1 : {
//		wFD[0] = w_GL[0];
//		wFD[1] = w_GL[1];
//	}
//    case 2 : {
//		wFD[0] = 0.5;
//		wFD[1] = 1.0;
//		wFD[2] = 0.5;
//	}
//    case 3 : {
//		wFD[0] = 0.5 - 0.1*sqrt(5.0);
//		wFD[1] = 0.5 + 0.1*sqrt(5.0);
//		wFD[2] = 0.5 + 0.1*sqrt(5.0);
//		wFD[3] = 0.5- 0.1*sqrt(5.0);
//	}
//}
//
//    D_SBP[0] = -0.5/wFD[0];
//cout << " N: " << N << "DSBP " << D_SBP[0] << "\n";
//    D_SBP[1] = 0.5/wFD[0];
//    D_SBP[ngl2-1] = 0.5/wFD[N];		//last entry		DOWNWIND
//    D_SBP[ngl2-2] = -0.5/wFD[N];		//second last entry		DOWNWIND
//
//    for (int i=1; i<ngl-1; i++)
//    {
//        const int id = i*ngl+i;
//        const int idp1 = id+1;
//        const int idm1 = id-1;
//        D_SBP[idp1]=0.5/wFD[i];
//        D_SBP[idm1]=-0.5/wFD[i];
//    };

};


void basis :: ConvertToModal(const dfloat Q_nodal[], dfloat Q_modal[])
{
    int Neq=3;
    dfloat * Qtmp = (dfloat*) calloc(ngl2*Nelem*Neq,sizeof(dfloat));

    for (int ie=0; ie<Nelem; ie++)
    {
        for (int i=0; i<ngl; i++)
        {
            for(int j=0; j<ngl; j++)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                Q_modal[id]= 0.0;
                id+=ngl2;
                Q_modal[id]= 0.0;
                id+=ngl2;
                Q_modal[id]= 0.0;

            }
        }

    }



// BRAUCHEN HIER Vinv^T * Q * Vinv
    for (int ie=0; ie<Nelem; ie++)
    {
        for (int i=0; i<ngl; i++)
        {
            for(int j=0; j<ngl; j++)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                for (int l=0; l<ngl; l++)
                {
                    int idLoc = ie*ngl2*Neq   +j*ngl+l;
                    int il = i*ngl+l;
//                Qtmp[id] += VdmInv(l+1,i+1) * Q_nodal[idLoc];
//                Qtmp[id+ngl2] += VdmInv(l+1,i+1) * Q_nodal[idLoc+ngl2];
//                Qtmp[id+ngl2+ngl2] += VdmInv(l+1,i+1) * Q_nodal[idLoc+ngl2+ngl2];

                    Qtmp[id] += VdmInv[il] * Q_nodal[idLoc];
                    Qtmp[id+ngl2] += VdmInv[il] * Q_nodal[idLoc+ngl2];
                    Qtmp[id+ngl2+ngl2] += VdmInv[il] * Q_nodal[idLoc+ngl2+ngl2];
                }

            }
        }

    }
//    cout <<"\n Checking VdmInv*Vdm \n";
    for (int ie=0; ie<Nelem; ie++)
    {
        for (int i=0; i<ngl; i++)
        {
            for(int j=0; j<ngl; j++)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                for (int l=0; l<ngl; l++)
                {
                    int idLoc = ie*ngl2*Neq   +l*ngl+i;
                    int jl = j*ngl+l;
//                Q_modal[id]+= Qtmp[idLoc] * VdmInv(l+1,j+1) ;
//                Q_modal[id+ngl2]+= Qtmp[idLoc+ngl2] * VdmInv(l+1,j+1) ;
//                Q_modal[id+ngl2+ngl2]+= Qtmp[idLoc+ngl2+ngl2] * VdmInv(l+1,j+1) ;

                    Q_modal[id]+= Qtmp[idLoc] * VdmInv[jl] ;
                    Q_modal[id+ngl2]+= Qtmp[idLoc+ngl2] * VdmInv[jl] ;
                    Q_modal[id+ngl2+ngl2]+= Qtmp[idLoc+ngl2+ngl2] * VdmInv[jl] ;
                }

            }
        }

    }

    free(Qtmp);

}

void basis :: EvaluteModalPolynomial(const dfloat Q_modal[], dfloat Q_nodal[])
{
    int Neq=3;
    dfloat LN1,dLN1;
    dfloat LN2,dLN2;
    dfloat * Qtmp = (dfloat*) calloc(ngl2*Nelem*Neq,sizeof(dfloat));

    for (int ie=0; ie<Nelem; ie++)
    {
        for (int i=0; i<ngl; i++)
        {
            for(int j=0; j<ngl; j++)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                Q_nodal[id]= 0.0;
                id+=ngl2;
                Q_nodal[id]= 0.0;
                id+=ngl2;
                Q_nodal[id]= 0.0;

            }
        }

    }

//                legendrePolynomialAndDerivative(l,x_GL(i+1),&LN1,&dLN1);
//                LN1 = LN1 * sqrt(l+0.5) ;
//                legendrePolynomialAndDerivative(k,x_GL(j+1),&LN2,&dLN2);
//                LN2 = LN2 * sqrt(k+0.5) ;
//
//                Q_nodal[id]+= Q_modal[idLoc] * LN1 * LN2 ;id+ngl2;idLoc+ngl2;
//                Q_nodal[id]+= Q_modal[idLoc] * LN1 * LN2 ;id+ngl2;idLoc+ngl2;
//                Q_nodal[id]+= Q_modal[idLoc] * LN1 * LN2 ;


    for (int ie=0; ie<Nelem; ie++)
    {
        for (int i=0; i<ngl; i++)
        {
            for(int j=0; j<ngl; j++)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;

                for (int l=0; l<ngl; ++l)   //loop thru Xi
                {
                    int idLoc = ie*ngl2*Neq   +j*ngl+l;

                    int il = i*ngl+l;

//                Qtmp[id]+= Vdm(l+1,i+1) * Q_modal[idLoc]   ;
//                Qtmp[id+ngl2]+= Vdm(l+1,i+1) * Q_modal[idLoc+ngl2]   ;
//                Qtmp[id+ngl2+ngl2]+= Vdm(l+1,i+1) * Q_modal[idLoc+ngl2+ngl2]  ;

                    Qtmp[id]+= Vdm[il] * Q_modal[idLoc]   ;
                    Qtmp[id+ngl2]+= Vdm[il] * Q_modal[idLoc+ngl2]   ;
                    Qtmp[id+ngl2+ngl2]+= Vdm[il] * Q_modal[idLoc+ngl2+ngl2]  ;
                }

            }
        }

    }

    for (int ie=0; ie<Nelem; ie++)
    {
        for (int i=0; i<ngl; i++)
        {
            for(int j=0; j<ngl; j++)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;

                for (int l=0; l<ngl; ++l)   //loop thru Xi
                {
                    int idLoc = ie*ngl2*Neq   +l*ngl+i;
                    int jl = j*ngl+l;

//                Q_nodal[id]+=  Qtmp[idLoc]   * Vdm(l+1,j+1) ;
//                Q_nodal[id+ngl2]+=  Qtmp[idLoc+ngl2]   * Vdm(l+1,j+1) ;
//                Q_nodal[id+ngl2+ngl2]+=  Qtmp[idLoc+ngl2+ngl2]   * Vdm(l+1,j+1) ;

                    Q_nodal[id]+=  Qtmp[idLoc]   * Vdm[jl] ;
                    Q_nodal[id+ngl2]+=  Qtmp[idLoc+ngl2]   * Vdm[jl] ;
                    Q_nodal[id+ngl2+ngl2]+=  Qtmp[idLoc+ngl2+ngl2]   * Vdm[jl] ;

                }
            }

        }
    }



    free(Qtmp);
}

//
//void basis :: SubCellAverageMatrix(){
//
//dfloat InterpolNodes[ngl];
//dfloat SubCellEdges[ngl+1];
//dfloat LagrangeAtInterpolNodes[ngl][ngl];
//dfloat InterpolLagrange[N+1];
//dfloat subCellSize;
//SubCellMat.resize(ngl,ngl);
//for(int i=1;i<=ngl;i++){
//        for(int j=1;j<=ngl;j++){
//    SubCellMat(i,j) = 0.0;}
//}
//
//SubCellEdges[0] = -1.0;
//for (int i=1; i<=ngl;i++){
//    SubCellEdges[i] = SubCellEdges[i-1] + w_GL(i);
//
//}
//
//
////GaussNodesAndWeights();
//for (int i=0;i<ngl;++i){
//    subCellSize = SubCellEdges[i+1] - SubCellEdges[i];
//
//    for (int j=0;j<ngl;j++){
//        InterpolNodes[j] = SubCellEdges[i] + subCellSize * (x_Gauss(j+1) + 1.0) /2.0;
//
//        LagrangeInterpolatingPolynomial(InterpolNodes[j],N,x_GL,w_bary,InterpolLagrange);
//        for (int l=0;l<ngl;l++){
//            LagrangeAtInterpolNodes[l][j] = InterpolLagrange[l];
//
//        }
//    }
//
//    for (int j=0;j<ngl;j++){
//        for (int l=0;l<ngl;l++){
//            SubCellMat(i+1,j+1) +=   subCellSize / (w_GL(i+1) *2.0) * LagrangeAtInterpolNodes[j][l] *w_Gauss(l+1);
//        }
//    }
//
//}
//
////for(int i=1;i<=ngl;i++){
////        for(int j=1;j<=ngl;j++){
////    cout << SubCellMat(i,j) << "  " ;
////
////    }
////cout <<"\n";
////}
//
//};


void basis :: ModalTrafoMatrix()
{
    dfloat LN,dLN;
    dfloat InterpolLagrange[N+1];
    dfloat * InterpolTmpMatrix = (dfloat*) calloc(ngl2,sizeof(dfloat));
    dfloat * VdmTmpMatrix = (dfloat*) calloc(ngl2,sizeof(dfloat));


//GaussNodesAndWeights();
    for (int i=0; i<ngl; ++i)
    {
        InterpolLagrange[i] = 0.0;
        for(int j=0; j<ngl; ++j)
        {
            int ij = i*ngl+j;
            int ji = j*ngl+i;
            legendrePolynomialAndDerivative(j,x_Gauss[i],&LN,&dLN);
            VdmTmpMatrix[ji]=LN * sqrt(j+0.5) * w_Gauss[i];
            legendrePolynomialAndDerivative(j,x_GL[i],&LN,&dLN);
            Vdm[ij] = LN * sqrt(j+0.5);
        }
    }
    for (int i=0; i<=N; i++)
    {

//LagrangeInterpolatingPolynomial(const dfloat x0,const int N,const dfloat x[],const dfloat w[],dfloat Polynomial[])
        LagrangeInterpolatingPolynomial(x_Gauss[i],N,x_GL,w_bary,InterpolLagrange);
        for (int j=0; j<ngl; j++)
        {
            InterpolTmpMatrix[i*ngl+j] = InterpolLagrange[j];
        }
    }
    for (int i=0; i<ngl; ++i)
    {
        for(int j=0; j<ngl; ++j)
        {
            for (int l=0; l<ngl; l++)
            {
                VdmInv[i*ngl+j] += VdmTmpMatrix[i*ngl+l] * InterpolTmpMatrix[l*ngl+j];
            }

        }
    }
//for (int i=1;i<=ngl;++i){
//    cout << "VdmTmpMatrix(:,"<<i<<") = ";
//    for(int j=1;j<=ngl;++j){
//        cout << VdmTmpMatrix(j,i) <<" ";
//    }
//    cout <<" \n";
//}
//
//for (int i=1;i<=ngl;++i){
//    cout << "InterpolTmpMatrix(:,"<<i<<") = ";
//    for(int j=1;j<=ngl;++j){
//        cout << InterpolTmpMatrix(j,i)<<" ";
//    }
//    cout <<" \n";
//}
//
//
//for (int i=1;i<=ngl;++i){
//    cout << "VdmInv(:,"<<i<<") = ";
//    for(int j=1;j<=ngl;++j){
//        cout << VdmInv(j,i)<<" ";
//    }
//    cout <<" \n";
//}
    free(InterpolTmpMatrix);
    free(VdmTmpMatrix);
};

void basis :: legendrePolynomialAndDerivative(const int N_in,const dfloat x,dfloat *L_N,dfloat *dL_N)
{



    dfloat L_N_2=0.0;
    dfloat L_N_1=0.0;
    dfloat dL_N_2=0.0;
    dfloat dL_N_1=0.0;
    *L_N=0.0;
    *dL_N=0.0;



//cout <<"TOL: "<<TOL;
    if (N_in==0)
    {
        *L_N=1.0;
        *dL_N=0.0;
    }
    else if (N_in==1)
    {
        *L_N=x;
        *dL_N=1.0;
    }
    else
    {
        L_N_2=1.0;
        L_N_1=x;
        dL_N_2=0.0;
        dL_N_1=1.0;

        for (int k=2; k<=N_in; k++)
        {
            *L_N=((2.0*k-1.0)/k)*x*L_N_1-(k-1.0)/k*L_N_2;
            *dL_N = dL_N_2 + (2.0*k -1.0)*L_N_1;
            L_N_2 = L_N_1;
            L_N_1 = *L_N;
            dL_N_2 = dL_N_1;
            dL_N_1 = *dL_N;

        };

    };
}




void basis :: setNelemLocal(const int input_Nelem)
{

    Nelem=input_Nelem;

};
void basis :: setNelemGlobal(const int input_Nelem_global)
{


    Nelem_global=input_Nelem_global;
};

void basis :: qAndLEvaluation(dfloat xGL,dfloat *q, dfloat *dq, dfloat *L_N)
{

    dfloat L_Nm2 = 1.0;
    dfloat L_Nm1 = xGL;
    dfloat dL_Nm2 = 0.0;
    dfloat dL_Nm1 = 1.0;
    dfloat dL_N=0.0;
    dfloat L_Np1=0.0;
    dfloat dL_Np1=0.0;

    for (int k=2; k<=N; k++)
    {
        *L_N=(2.0*k-1.0)/k * xGL *L_Nm1 - (k-1.0)/k*L_Nm2;
        dL_N=dL_Nm2 + (2.0*k-1.0)*L_Nm1;
        if (k<N)
        {
            L_Nm2=L_Nm1;
            L_Nm1=*L_N;
            dL_Nm2 = dL_Nm1;
            dL_Nm1=dL_N;
        };
    };

    int k=N+1;
    L_Np1=(2.0*k-1.0)/k *xGL* *L_N - (k-1.0) /k * L_Nm1;
    dL_Np1 = dL_Nm1 + (2.0*k-1.0)* *L_N;
    *q = L_Np1 - L_Nm1;
    *dq = dL_Np1 - dL_Nm1;


};



void basis :: LGLNodesAndWeights()
{

    const int newton_it = 4;
    dfloat q=0.0;
    dfloat del_q=0.0;
    dfloat L_N=0.0;
    dfloat TOL=4.*pow(10.,-16);
    dfloat delta=0.0;

//cout <<"TOL: "<<TOL;

    if (N==1)
    {
        x_GL[0]=-1.0;
        w_GL[0]=1.0;
        x_GL[1]=1.0;
        w_GL[1]=w_GL[0];
    }
    else if (N==2)
    {

        x_GL[0]=-1.0;
        w_GL[0]=2.0/(N*(N+1.0));
        x_GL[ngl-1]=1.0;
        w_GL[ngl-1]=w_GL[0];
        qAndLEvaluation(0.0,&q,&del_q,&L_N);
        x_GL[1] 	=0.0;
        w_GL[1]	= 2.0/(N*(N+1)*pow(L_N,2));

    }
    else
    {
        x_GL[0]=-1.0;
        w_GL[0]=2.0/(N*(N+1.0));
        x_GL[ngl-1]=1.0;
        w_GL[ngl-1]=w_GL[0];

        for (int j=1; j<=floor((N)/2.0); j++)
        {
            x_GL[j] = -cos((j+1.0/4.0)*M_PI/N - 3.0/(8.0*N*M_PI)*1.0/(j+1.0/4.0));

            for (int k=0; k<=newton_it; k++)
            {
                qAndLEvaluation(x_GL[j],&q,&del_q,&L_N);

                delta = -q/del_q;
                x_GL[j] = x_GL[j] + delta;
                if(abs(delta)<=TOL*abs(x_GL[j]))break;
            };


            qAndLEvaluation(x_GL[j],&q,&del_q,&L_N);
            x_GL[ngl-j-1] 	=-x_GL[j];
            w_GL[j]	= 2.0/(N*(N+1)*pow(L_N,2));
            w_GL[ngl-j-1]	=w_GL[j];
        };
        //if the middle point 0.0 is included, also calculated it.
        if(N%2==0)
        {
            qAndLEvaluation(0.0,&q,&del_q,&L_N);
            x_GL[N/2] 	=0.0;
            w_GL[N/2]	= 2.0/(N*(N+1)*pow(L_N,2));
        };

    };

};


void basis :: GaussNodesAndWeights()
{

    const int newton_it = 4;

    dfloat L1=0.0;
    dfloat L2=0.0;
    dfloat TOL=4.*pow(10.,-16);
    dfloat delta=0.0;


    if (N==0)
    {
        x_Gauss[0]=0.0;
        w_Gauss[0]=2.0;
    }
    else if (N==1)
    {
        x_Gauss[0]=-sqrt(1.0/3.0);
        w_Gauss[0]=1.0;
        x_Gauss[1]=-x_Gauss[0];
        w_Gauss[1]=w_Gauss[0];
    }
    else
    {

        for (int j=0; j<floor((N+1.0)/2.0); j++)
        {
            x_Gauss[j] = -cos((2.0*j+1.0)/(2.0*N+2.0)*M_PI);

            for (int k=0; k<=newton_it; k++)
            {
                legendrePolynomialAndDerivative(N+1,x_Gauss[j],&L1,&L2);
                delta = -L1/L2;
                x_Gauss[j] = x_Gauss[j] + delta;
                if(abs(delta)<=TOL*abs(x_Gauss[j]))break;
            };


            legendrePolynomialAndDerivative(N+1,x_Gauss[j],&L1,&L2);

            x_Gauss[ngl-j-1] 	=-x_Gauss[j];
            w_Gauss[j]	= 2.0/((1.0-pow(x_Gauss[j],2))*pow(L2,2));
            w_Gauss[ngl-j-1]	=w_Gauss[j];
        };
        //if the middle point 0.0 is included, also calculated it.
        if(N%2==0)
        {
            legendrePolynomialAndDerivative(N+1,0.0,&L1,&L2);
            x_Gauss[N/2] 	=0.0;
            w_Gauss[N/2]	= 2.0/(pow(L2,2));
        };

    }

};




void basis :: LagrangeInterpolatingPolynomial(const dfloat x0,const int N,const dfloat * x,const dfloat * w,dfloat Polynomial[])
{

    bool xMatchesNode=0;
    dfloat TOL=4.*pow(10.,-16);
    dfloat t =0.0;
    for (int j=0; j<N+1; j++)
    {
        Polynomial[j]=0.0;
        if (abs(x0 - x[j])<TOL)
        {
            xMatchesNode=1;
            Polynomial[j]=1;
        }
    }
    if(!xMatchesNode)
    {
        dfloat s=0.0;
        for(int j=0; j<N+1; j++)
        {
            t = w[j]/(x0-x[j]);
            Polynomial[j]=t;
            s=s+t;
        }
        for(int j=0; j<N+1; j++)
        {
            Polynomial[j] = Polynomial[j]/s;
        }
    }

}


void basis :: BarycentricWeights()
{


    for(int i=0; i<ngl; i++)
    {
        w_bary[i]=1.0;
    };


    for (int j=1; j<ngl; j++)
    {
        for (int k=0; k<j; k++)
        {
            w_bary[k]=w_bary[k]*(x_GL[k]-x_GL[j]);
            w_bary[j]=w_bary[j]*(x_GL[j]-x_GL[k]);
        };
    };

    for (int j=0; j<ngl; j++)
    {
        w_bary[j]=1.0/w_bary[j];
    };

};

void basis :: PolynomialDerivativeMatrix()
{


    for (int i=0; i<ngl; i++)
    {
        for (int j=0; j<ngl; j++)
        {
            if (i!=j)
            {
                const int id = i*ngl+j;
                D[id]=(w_bary[j]/w_bary[i])*(1.0/(x_GL[i]-x_GL[j]));
                D[i*ngl+i]=D[i*ngl+i] -D[id];
            };
        };
    };
};





void basis :: L2Norm(const dfloat Q[],const dfloat Q_exakt[],const dfloat J[],dfloat L2Error[])
{

    L2Error[0]=0.0;
    L2Error[1]=0.0;
    L2Error[2]=0.0;
//cout << "Calculating L2 norm! Nelem:" << Nelem_global <<" J[xid] " << J[2000] << " Q_exakt[id] " << Q_exakt[2000] << " ngl " << ngl <<  " Neq" << Neq <<"\n";
    for (int ie=0; ie<Nelem_global; ie++)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                int xid = ie*ngl2   +j*ngl+i;

                L2Error[0]  +=  pow(Q[id] - Q_exakt[id],2) /J[xid]* w_GL[i]* w_GL[j];
//		if (L2Error[0]<0){
//			cout << "L2Error[0] " << L2Error[0] << " Q[id] " << Q[id] <<" Q_exakt[id] " << Q[id] << " J[xid] " << J[xid] <<"\n";
//		}
                id+=ngl2;
                L2Error[1]  +=  pow(Q[id] - Q_exakt[id],2) /J[xid]* w_GL[i]* w_GL[j];
                id+=ngl2;
                L2Error[2]  +=  pow(Q[id] - Q_exakt[id],2) /J[xid]* w_GL[i]* w_GL[j];
            }
        }
    }


};

void basis :: LinfNorm(const dfloat Q[],const dfloat Q_exakt[],dfloat LinfError[])
{

    LinfError[0]=0.0;
    LinfError[1]=0.0;
    LinfError[2]=0.0;
    int maxErrorElement=-1;
    for (int ie=0; ie<Nelem_global; ie++)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;

		if (fabs(Q[id] - Q_exakt[id])>LinfError[0]){
			maxErrorElement=ie;
			//cout << "estimated value: " << Q[id] << " exact value " << Q_exakt[id] << "\n";
		}
                LinfError[0]  =  max(LinfError[0],fabs(Q[id] - Q_exakt[id]));
                id+=ngl2;
                LinfError[1]  =  max(LinfError[1],fabs(Q[id] - Q_exakt[id]));
                id+=ngl2;
                LinfError[2]  =  max(LinfError[2],fabs(Q[id] - Q_exakt[id]));

            }
        }
    }

	cout << "maximum h error was found in element : " << maxErrorElement << "\n";

};


void basis :: calcElementSizes(const dfloat J[],dfloat EleSize[],dfloat * minEleSize)
{
    dfloat tmpMinEleSize;
    for (int ie=0; ie<Nelem; ie++)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {

                int xid = ie*ngl2   +j*ngl+i;
                EleSize[ie]  +=  1.0 /J[xid]* w_GL[i]* w_GL[j];


            }
        }
        if (ie==0)
        {
            tmpMinEleSize=EleSize[ie];
        }
        else
        {
            tmpMinEleSize = min(tmpMinEleSize,EleSize[ie]);
        }
//        EleSize[ie] = sqrt(EleSize[ie]);
    }
    *minEleSize = sqrt(tmpMinEleSize);
};



void basis :: calcEntropyDelta(const dfloat g_const,const dfloat Q[],const dfloat Q_init[],const dfloat b[],const dfloat J[],dfloat *EntropyDelta,dfloat * relEntropyDelta)
{

    *EntropyDelta=0.0;
    dfloat TotalEntropy_Final = 0.0;
    dfloat TotalEntropy_Init = 0.0;

    for (int ie=0; ie<Nelem_global; ie++)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                int xid = ie*ngl2   +j*ngl+i;

                dfloat E_final=0.0;
                calcEntropyPointwise(g_const,Q[id],Q[id+ngl2],Q[id+ngl2+ngl2],b[xid],&E_final);

                dfloat E_init=0.0;
                calcEntropyPointwise(g_const,Q_init[id],Q_init[id+ngl2],Q_init[id+ngl2+ngl2],b[xid],&E_init);
//            cout <<"Entropy final: "<<E_final <<"\n";
//            cout <<"Entropy init: "<<E_init <<"\n";
                TotalEntropy_Final  +=  E_final /J[xid]* w_GL[i]* w_GL[j];
                TotalEntropy_Init   +=  E_init /J[xid]* w_GL[i]* w_GL[j];

            }
        }
    }

    *EntropyDelta=TotalEntropy_Final - TotalEntropy_Init;
    *relEntropyDelta=*EntropyDelta/TotalEntropy_Init;
};

void basis :: calcTotalEntropy(const dfloat g_const,const dfloat Q[],const dfloat b[],const dfloat J[],dfloat *TotalEntropy)
{

    dfloat r_TotalEntropy= 0.0;

    for (int ie=0; ie<Nelem_global; ie++)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                int xid = ie*ngl2   +j*ngl+i;

                dfloat Entropy=0.0;
                calcEntropyPointwise(g_const,Q[id],Q[id+ngl2],Q[id+ngl2+ngl2],b[xid],&Entropy);

//            cout <<"Entropy final: "<<E_final <<"\n";
//            cout <<"Entropy init: "<<E_init <<"\n";
                r_TotalEntropy  +=  Entropy /J[xid]* w_GL[i]* w_GL[j];


            }
        }
    }

    *TotalEntropy=r_TotalEntropy ;

};

void basis :: calcEntropyPointwise(const dfloat g_const,const dfloat h,const dfloat hu,const dfloat hv,const dfloat b,dfloat *Entropy)
{
    dfloat u;
    dfloat v;
    if (h>pow(10,-12))
    {
        u = hu/h;
        v = hv/h;

    }
    else
    {
        u = 0.0;
        v = 0.0;

    }



    *Entropy=0.5*h*(u*u+v*v) +0.5*g_const*h*h +  g_const*h*b;

};


void basis :: checkConservation(const dfloat Q[],const dfloat Q_init[],const dfloat J[],dfloat *MassDelta, dfloat *relMassError)
{

    *MassDelta=0.0;
    dfloat TotalMass_Final = 0.0;
    dfloat TotalMass_Init = 0.0;

    for (int ie=0; ie<Nelem_global; ie++)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                int xid = ie*ngl2   +j*ngl+i;

//            cout <<"Entropy final: "<<E_final <<"\n";
//            cout <<"Entropy init: "<<E_init <<"\n";
                TotalMass_Final  +=  Q[id] /J[xid]* w_GL[i]* w_GL[j];
                TotalMass_Init   +=  Q_init[id] /J[xid]* w_GL[i]* w_GL[j];

            }
        }
    }

    *MassDelta=TotalMass_Final - TotalMass_Init;
    *relMassError = *MassDelta/TotalMass_Init;
};


void basis :: calcTotalMass(const dfloat Q[],const dfloat J[],dfloat *TotalMass)
{

    dfloat r_TotalMass= 0.0;

    for (int ie=0; ie<Nelem_global; ie++)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                int xid = ie*ngl2   +j*ngl+i;

                r_TotalMass  +=  Q[id] /J[xid]* w_GL[i]* w_GL[j];


            }
        }
    }

    *TotalMass=r_TotalMass ;

};


void basis :: UpdateMaximumFriction(const dfloat FrictionForPlot[],dfloat maximumFriction[])
{


    for (int ie=0; ie<Nelem_global; ie++)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2*(Neq-1)   +j*ngl+i;

                maximumFriction[id]       =     max(abs(FrictionForPlot[id]),maximumFriction[id]);
                maximumFriction[id+ngl2]  =     max(abs(FrictionForPlot[id+ngl2]),maximumFriction[id+ngl2]);

            }
        }
    }


};



