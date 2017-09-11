#include "RungeKutta.h"

RungeKutta::RungeKutta(int RkOrder, int rkSSP)
{
    switch(RkOrder){
case 3:
    rkStages=3;
//    Coeffs.resize(3,3);
     CoeffsA = (dfloat*) calloc(rkStages,sizeof(dfloat));
     CoeffsB = (dfloat*) calloc(rkStages,sizeof(dfloat));
     CoeffsC = (dfloat*) calloc(rkStages,sizeof(dfloat));
    if (rkSSP){
        CoeffsA[0]=1.0;
        CoeffsA[1]=3.0/4.0;
        CoeffsA[2]=1.0/3.0;

        CoeffsB[0]=0.0;
        CoeffsB[1]=1.0/4.0;
        CoeffsB[2]=2.0/3.0;

        CoeffsC[0]=1.0;
        CoeffsC[1]=1.0/4.0;
        CoeffsC[2]=2.0/3.0;


    }else{
        CoeffsA[0]=0.0;
        CoeffsA[1]=-5.0/9.0;
        CoeffsA[2]=-153.0/128.0;

        CoeffsB[0]=0.0;
        CoeffsB[1]=1.0/3.0;
        CoeffsB[2]=3.0/4.0;

        CoeffsC[0]=1.0/3.0;
        CoeffsC[1]=15.0/16.0;
        CoeffsC[2]=8.0/15.0;



    }

    break;
case 4:
    rkStages=5;
//    Coeffs.resize(3,5);
     CoeffsA = (dfloat*) calloc(rkStages,sizeof(dfloat));
     CoeffsB = (dfloat*) calloc(rkStages,sizeof(dfloat));
     CoeffsC = (dfloat*) calloc(rkStages,sizeof(dfloat));
    CoeffsA[0]=0.0;
    CoeffsA[1]=-(567301805773.0/1357537059087.0);
    CoeffsA[2]=-(2404267990393.0/2016746695238.0);
    CoeffsA[3]=-(3550918686646.0/2091501179385.0);
    CoeffsA[4]=-(1275806237668.0/842570457699.0);

    CoeffsB[0]=0.0;
    CoeffsB[1]=(1432997174477.0/9575080441755.0);
    CoeffsB[2]=(2526269341429.0/6820363962896.0);
    CoeffsB[3]=(2006345519317.0/3224310063776.0);
    CoeffsB[4]=(2802321613138.0/2924317926251.0);

    CoeffsC[0]=(1432997174477.0/9575080441755.0);
    CoeffsC[1]=(5161836677717.0/13612068292357.0);
    CoeffsC[2]=(1720146321549.0/2090206949498.0);
    CoeffsC[3]=(3134564353537.0/4481467310338.0);
    CoeffsC[4]=(2277821191437.0/14882151754819.0);
    break;


    }
}

RungeKutta::~RungeKutta()
{
    //dtor
}
