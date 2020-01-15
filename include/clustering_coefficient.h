#ifndef CLUSTERING_COEFFICIENT_H
#define CLUSTERING_COEFFICIENT_H

#include <igraph/igraph.h>

int graph_clustering_coefficient(igraph_t&, 
								 igraph_integer_t, 
								 igraph_real_t&, 
								 igraph_real_t&);

#endif
