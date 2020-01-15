#ifndef RANDOM_MATRIX_THEORY_H
#define RANDOM_MATRIX_THEORY_H

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

#include "spline.h"
#include "alglib/specialfunctions.h"  // ALGLIB header that contains ChiSquareCDistribution
#include "igraph_ext.h"


int random_matrix_theory(igraph_t&, 
                         igraph_integer_t&, 
                         double&, 
                         double&, 
                         double&, 
                         double&);



#endif