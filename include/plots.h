#ifndef plots_H
#define plots_H
#include <iostream>
#include <sstream>
#include <string>
#include "Constants.hpp"
using namespace Constants;



void PlotSolution(const int, const int,const int, const dfloat[],const dfloat[],const dfloat[],const dfloat[],const int);
void PlotViscosity(const int, const int,const int, const dfloat[],const dfloat[],const dfloat[],const dfloat[],const int);
void PlotViscoseParameter(const int , const int , const dfloat [], const dfloat [], const dfloat [], const int );
#endif // plots_H

