#ifndef RASCALS_HEADERS3D
#define RASCALS_HEADERS3D

#include <sys/types.h>
#include <sys/stat.h>
#include <execinfo.h>
#include <cxxabi.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <stdio.h>
#include <math.h>

#include <algorithm>
#include <iostream>
#include <ostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>

using namespace std;

#define myMin(a,b) ( (a<b) ? a : b )
#define myMax(a,b) ( (a>b) ? a : b )
#define norm(a) ( (a < 0) ? -a : a)

#define fmatrix matrix<dfloat>
#define imatrix matrix<int>
#define fzero   0

#endif
