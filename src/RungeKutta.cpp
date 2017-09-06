#include "RungeKutta.h"

RungeKutta::RungeKutta(int RkOrder, bool rkSSP)
{
    switch(RkOrder){
case 3:
    rkStages=3;
    Coeffs.resize(3,3);
    if (rkSSP){
        Coeffs(1,1)=1.0;
        Coeffs(1,2)=3.0/4.0;
        Coeffs(1,3)=1.0/3.0;

        Coeffs(2,1)=0.0;
        Coeffs(2,2)=1.0/4.0;
        Coeffs(2,3)=2.0/3.0;

        Coeffs(3,1)=1.0;
        Coeffs(3,2)=1.0/4.0;
        Coeffs(3,3)=2.0/3.0;


    }else{
        Coeffs(1,1)=0.0;
        Coeffs(1,2)=-5.0/9.0;
        Coeffs(1,3)=-153.0/128.0;

        Coeffs(2,1)=0.0;
        Coeffs(2,2)=1.0/3.0;
        Coeffs(2,3)=3.0/4.0;

        Coeffs(3,1)=1.0/3.0;
        Coeffs(3,2)=15.0/16.0;
        Coeffs(3,3)=8.0/15.0;



    }

    break;
case 4:
    rkStages=5;
    Coeffs.resize(3,5);
    Coeffs(1,1)=0.0;
    Coeffs(1,2)=-(567301805773.0/1357537059087.0);
    Coeffs(1,3)=-(2404267990393.0/2016746695238.0);
    Coeffs(1,4)=-(3550918686646.0/2091501179385.0);
    Coeffs(1,5)=-(1275806237668.0/842570457699.0);

    Coeffs(2,1)=0.0;
    Coeffs(2,2)=(1432997174477.0/9575080441755.0);
    Coeffs(2,3)=(2526269341429.0/6820363962896.0);
    Coeffs(2,4)=(2006345519317.0/3224310063776.0);
    Coeffs(2,5)=(2802321613138.0/2924317926251.0);

    Coeffs(3,1)=(1432997174477.0/9575080441755.0);
    Coeffs(3,2)=(5161836677717.0/13612068292357.0);
    Coeffs(3,3)=(1720146321549.0/2090206949498.0);
    Coeffs(3,4)=(3134564353537.0/4481467310338.0);
    Coeffs(3,5)=(2277821191437.0/14882151754819.0);
    break;


    }
}

RungeKutta::~RungeKutta()
{
    //dtor
}
