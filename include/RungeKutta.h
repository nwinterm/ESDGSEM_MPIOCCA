#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "Constants.hpp"
using namespace Constants;

class RungeKutta
{
public:
    RungeKutta(int,int);
    int rkStages;
    dfloat * CoeffsA;
    dfloat * CoeffsB;
    dfloat * CoeffsC;
    virtual ~RungeKutta();



protected:

private:
};

#endif // RUNGEKUTTA_H
