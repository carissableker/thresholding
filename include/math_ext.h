#ifndef MATH_EXT_H
#define MATH_EXT_H
 
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

///////////////////////////////////////////////////////////////////////////////
//     Math/Stat functions                                                   //
///////////////////////////////////////////////////////////////////////////////

// Returns the vector of differences between first and 
// last elements of the windows of size n in x
// from igraph_vector_t to  std::vector
int rolling_difference_igraph(igraph_vector_t &, std::vector<double> &, int);


// Returns the vector of differences between first and 
// last elements of the windows of size n in x
// fromstd::vector to  std::vector
int rolling_difference(std::vector<double> &, std::vector<double> &, int);

// Median of a vector
double median(std::vector<double>);

// Mean/average of a vector
double mean(std::vector<double>);

// Standard deviation of a vector
double stddev(std::vector<double>, double dof=1);

// get the exponent to pow value to make a double an int
// Stephen Grady
    
int get_precision(double);

// Range from l to u, incrementing by increment
std::vector<double> range(double, double, double);

// Empirical Cumulative Distribution Function (ecdf)
// given x and t, evaluate eCDF of x at t
std::vector<double> ecdf(std::vector<double>, std::vector<double>);

double sign(double);


// Poisson pdf
double poisson(double, double);

// GOE pdf (actually Wigner-Dyson)
double goe(double, double);


#endif
