#ifndef SPLINE_H
#define SPLINE_H

#include <math.h>     // pow, sqrt, fabs
#include <algorithm>  // std::nth_element, std::min_element, std::max_element
#include <cmath>      // std::copysign
#include <iostream>   // std::cout, std::cerr, std::endl

#include <math_ext.h> 


// Accurate Monotonicity Preserving Cubic Interpolation
// In this case x, y is a cdf, and always monotone increasing
// https://epubs-siam-org.proxy.lib.utk.edu/doi/pdf/10.1137/0904045
std::vector<double> spline(std::vector<double>, std::vector<double>, std::vector<double>);

#endif