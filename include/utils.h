#ifndef UTILS_H
#define UTILS_H



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

// Results IO
int output_results(std::string&, std::string&);

#endif