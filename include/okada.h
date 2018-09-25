
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
    void okadamapFull(int &xsize, int &ysize,
                            dfloat &xlower,
                            dfloat &xupper,
                            dfloat &ylower,
                            dfloat &yupper);

private:

};






#endif // okada_H
