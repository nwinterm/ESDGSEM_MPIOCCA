
#ifndef okada_H
#define okada_H

#include "MeshPartitioning.h"
#include "Constants.hpp"
using namespace Constants;


class okada
{
public:
    okada(const int);
    void okadamapFull(const int,const int,const dfloat[],const dfloat[], dfloat[]);

private:
    int NumElements;
dfloat strike_slip (const dfloat x1,
                    const dfloat x2,
                    const dfloat x3,
                    const dfloat y1,
                    const dfloat y2,
                    const dfloat dp,
                    const dfloat dd);
dfloat dip_slip (const dfloat x1,
                 const dfloat x2,
                 const dfloat x3,
                 const dfloat y1,
                 const dfloat y2,
                 const dfloat dp,
                 const dfloat dd);
};






#endif // okada_H
