#ifndef MATH_EXT_H
#define MATH_EXT_H

#include <vector>     // std::vector
#include <algorithm>  // std::nth_element, std::min/amx_element, size_t,  partial sort
#include <math.h>     // pow, sqrt, fabs
#include <cmath>      // std::copysign
#include <stdexcept>  // std::invalid_argument
#include <fstream>    // ofstream
#include <sstream>    // ostringstream
#include <numeric>    // iota

#include <igraph/igraph.h>

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


double fisher_transform(double, int);


std::vector<size_t> argsort(std::vector<double>, int);

#endif
