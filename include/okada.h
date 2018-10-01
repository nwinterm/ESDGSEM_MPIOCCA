
#ifndef okada_H
#define okada_H

#include "MeshPartitioning.h"
#include "Constants.hpp"
using namespace Constants;


class okada
{
public:
    okada(const int, const dfloat);
    void set_ics_from_okada_model();
    void okadamapFull(const int,const int,const dfloat[],const dfloat[], dfloat[]);

private:

};






#endif // okada_H
