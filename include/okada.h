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



#ifndef okada_H
#define okada_H

#include "MeshPartitioning.h"
#include "Constants.hpp"
using namespace Constants;


class okada
{
public:
    okada(const int);
    void okadamapFull(const int,const dfloat[],const dfloat[], dfloat[]);

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
