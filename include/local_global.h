#ifndef LOCAL_GLOBAL_H
#define LOCAL_GLOBAL_H

#include <igraph/igraph.h>
#include <vector>     // std::vector


#include "math_ext.h"
#include "spectral_methods.h"


int local_global_method(igraph_t&,
				 double,
				 double,
				 double,
                 int,
                 int,
                 std::string&
				 );

#endif