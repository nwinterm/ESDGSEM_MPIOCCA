# ESDGSEM_MPIOCCA

ESDGSEM_MPIOCCA is a numerical solver for the Shallow Water Equations.

It is based on a special Discontinuous Galerkin Spectral Element
Method (DGSEM) and works on unstructured quadrilateral meshes with possibly curved elements. Volume integrals are computed by a flux differencing approach,
leading to an entropy stable and well-balanced scheme.


The main contributers are
* from Mathematical Institute, University of Cologne: Niklas Wintermeyer, Gregor Gassner, Andrew Winters


## License 

ESDGSEM_MPIOCCA is released under the terms of the GNU General Public License v3.0. 
For the full license terms see the included license file [LICENSE.md](LICENSE.md).

## Used libraries

ESDGSEM_MPIOCCA uses several external libraries from open source projects, including:
* [OCCA](https://libocca.org/#/)
* [Open MPI](https://www.open-mpi.org/)


