#ifndef MPI_Communication_H
#define MPI_Communication_H


#include "Constants.hpp"
#include "MPI_setup.h"
#include "MeshPartitioning.h"
using namespace Constants;

void ShareInputData(MPI_setup, int *,dfloat *,dfloat *, dfloat *,dfloat *,int *,int *,dfloat *,dfloat *,dfloat *,
                    dfloat *,int *,int *,int *,int *, int *,int *, int *, int *, int *,int*,int*,int*,int*,int*);




void CollectSolution(MPI_setup, const MeshPartitioning,  const dfloat [], dfloat []);
void SendSolution(MPI_setup, const MeshPartitioning, const dfloat []);
void CollectEdgeDataMPI(MPI_setup, const MeshPartitioning, dfloat [], dfloat []);
void CollectEdgeDataMPI_Bonly(MPI_setup, const MeshPartitioning, dfloat [], dfloat []);

void CollectViscoseEdgeDataMPI(MPI_setup, const MeshPartitioning, dfloat [], dfloat [], dfloat [], dfloat [], dfloat [], dfloat []);


void GetGlobalLambdaMax(MPI_setup, const MeshPartitioning,const dfloat [], dfloat *);
void GetGlobalViscParaMax(MPI_setup, const MeshPartitioning, const dfloat [], dfloat * );
void GetGlobalMinEleSize(MPI_setup, const MeshPartitioning,const dfloat, dfloat *);




void CollectViscPara(MPI_setup, const MeshPartitioning, const dfloat[], dfloat[]);
void SendViscPara(MPI_setup, const MeshPartitioning, const dfloat[]);

void CollectViscosity(MPI_setup, const MeshPartitioning,  const dfloat [],const dfloat [], dfloat [], dfloat []);
void SendViscosity(MPI_setup, const MeshPartitioning, const dfloat [], const dfloat []);
#endif // SW1D_H

