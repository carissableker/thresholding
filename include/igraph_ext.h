#ifndef IGRAPH_EXT_H
#define IGRAPH_EXT_H

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
//     IO functions                                                          //
///////////////////////////////////////////////////////////////////////////////

// Read in graph
int read_graph(std::string &, igraph_t&);

///////////////////////////////////////////////////////////////////////////////
//     Graph functions                                                       //
///////////////////////////////////////////////////////////////////////////////

// Threshold graph
// by removing edges with abs weight less than "t" 
// and subsequently vertices with no neighbours 
int threshold_graph(double, igraph_t &);

// Identify largest connected component of the graph and induce
// Also return number of CC
int largest_connected_component(igraph_t &, igraph_t &, igraph_integer_t &);

// Fiedler vector: eigen-vector corresponding to first non-zero eigen-value
// Assume connected graph -> 2nd eigenvector
int Fiedler_vector(igraph_t &, igraph_vector_t &, igraph_real_t &);

// See igraph_get_adjacency for unweighted graph
int get_weighted_adjacency(igraph_t &, igraph_matrix_t &);


#endif
