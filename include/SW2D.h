#ifndef SW2D_H
#define SW2D_H

#include "MeshPartitioning.h"
#include "Constants.hpp"
using namespace Constants;





void InitQ(const int, const MeshPartitioning,const int,const int,const int,const int,const dfloat[],const dfloat[], dfloat[],const dfloat,const dfloat[],const dfloat);
void InitQNodal(const int,const dfloat,const dfloat, dfloat [],const dfloat,const dfloat,const dfloat);
void InitB(const int, const MeshPartitioning,const int,const int,const int,const int,const dfloat[],const dfloat[],dfloat[]);
void CalcBDerivatives(const int,const int,const int,const dfloat,const dfloat[],const dfloat[],const dfloat[],const dfloat[], const dfloat[], const dfloat[], const dfloat[], const dfloat[], dfloat[],dfloat[], const dfloat[]);


#endif // SW1D_H
