
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
void PlotMaximumElevation(const int, const int, const dfloat[],const dfloat[],const dfloat[]);

#endif // plots_H

