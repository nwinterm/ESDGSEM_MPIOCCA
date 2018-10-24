#ifndef plots_H
#define plots_H
#include <iostream>
#include <sstream>
#include <string>
#include "Constants.hpp"
using namespace Constants;



void PlotSolution(const int, const int,const int, const dfloat[],const dfloat[],const dfloat[],const dfloat[],const int,const dfloat);
void PlotFriction(const int, const int,const int, const dfloat[],const dfloat[],const dfloat[],const int);
void PlotSolutionWithExact(const int, const int,const int, const dfloat[],const dfloat[],const dfloat[],const dfloat[],const int,const dfloat[]);
void PlotViscosity(const int, const int,const int, const dfloat[],const dfloat[],const dfloat[],const dfloat[],const int);
void PlotViscoseParameter(const int, const int, const dfloat [], const dfloat [], const dfloat [], const int );
void PlotEntropy(const int, const dfloat [], const dfloat []);
void PlotTimeSeries(const int, const dfloat[], const dfloat[],const dfloat[],const dfloat[],const dfloat[],const dfloat[],const dfloat[],const dfloat[]);
void PlotMass(const int, const dfloat [], const dfloat []);

void PlotArrivalTimings(const int, const int, const dfloat[],const dfloat[],const dfloat[]);
#endif // plots_H

