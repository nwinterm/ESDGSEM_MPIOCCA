#ifndef READINPUT_H
#define READINPUT_H


#include "Constants.hpp"
using namespace Constants;


void ReadInputFile(int *N,
                   string *meshFile,
                   dfloat *CFL,
                   dfloat *DFL,
                   dfloat *T,
                   dfloat *g_const,
                   int *ArtificialViscosity,
                   int *PositivityPreserving,
                   dfloat *PosPresTOL,
                   dfloat *epsilon_0,
                   dfloat *sigma_min,
                   dfloat *sigma_max,
                   int *PlotVar,
                   int *EntropyPlot,
                   int *NumPlots,
                   int *NumTimeChecks,
                   int *Testcase,
                   int *ES,
                   int *NumFlux,
                   int *FluxDifferencing,
                   int *Cartesian,
                   int * rkorder,
                   int * rkSSP,
                   int * NEpad,
                   int * NEsurfpad,
                   int * Nedgepad,
			int * NavgPad,
                   int * KernelVersion,
                   int * KernelVersionSTD,
			int * DiscBottom);

void ReadCartesianData(const int, const int, dfloat *xL,dfloat *xR,dfloat *yL,dfloat *yR, int *NelemX, int *NelemY,bool *PeriodicBD_X,bool *PeriodicBD_Y);

#endif // READINPUT_H
