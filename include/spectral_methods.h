#ifndef SPECTRAL_METHODS_H
#define SPECTRAL_METHODS_H

#include <igraph/igraph.h>

#include <vector>     // std::vector
#include <iostream>   // std::cout, std::cerr, std::endl
#include <fstream>    // fopen, fclose (to read igraph)
#include <algorithm>  // std::nth_element, std::min_element, std::max_element
#include <math.h>     // pow, sqrt, fabs
#include <getopt.h>   // commandline argument parsing
#include <stdlib.h>   // atoi, atof
#include <sstream>    // stringstream
#include <cmath>      // std::copysign

#include "igraph_ext.h"
#include "utils.h"
#include "math_ext.h"

int spectral_methods(igraph_t&,
					 int,
					 int,
					 igraph_real_t&, 
					 int&);

#endif
