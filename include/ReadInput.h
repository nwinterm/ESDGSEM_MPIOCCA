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
                   int * DiscBottom,
                   int * ReadInBottom,
                   int * PartialDryTreatment,
                   int * FrictionTerms,
                   int * ConvertToKM,
                   int * calcArrivalTimes,
                   int * createTimeSeries,
                   int * calcMaximumElevation);

void ReadCartesianData(const int, const int, dfloat *xL,dfloat *xR,dfloat *yL,dfloat *yR, int *NelemX, int *NelemY,bool *PeriodicBD_X,bool *PeriodicBD_Y);
void WriteFullMesh(const int NumNodes, const dfloat *x,const dfloat *y);
void ReadFullMesh(const int NumNodes,const int Nelem, const int N, dfloat *b, dfloat *h_0);

void FindElementID(const int NumNodes, const dfloat *x,const dfloat *y, const dfloat lonToFind, const dfloat latToFind, int* coordID);
#endif // READINPUT_H
