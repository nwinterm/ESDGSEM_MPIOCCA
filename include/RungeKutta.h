#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "Constants.hpp"
using namespace Constants;

class RungeKutta
{
    public:
        RungeKutta(int,bool);
        int rkStages;
        fmatrix Coeffs;
        virtual ~RungeKutta();



    protected:

    private:
};

#endif // RUNGEKUTTA_H
