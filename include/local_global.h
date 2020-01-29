#ifndef LOCAL_GLOBAL_H
#define LOCAL_GLOBAL_H

#include <fstream>    // ofstream
#include <vector>     // std::vector

#include <igraph/igraph.h>

#include "math_ext.h"
#include "spectral_methods.h"


int local_global_pruning(igraph_t&,
				 double,
				 igraph_t&,
				 double&);

int local_global_method(igraph_t&,
				 double,
				 double,
				 double,
                 int,
                 int,
                 std::string&
				 );

#endif