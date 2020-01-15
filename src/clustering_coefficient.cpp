#include "clustering_coefficient.h"

int graph_clustering_coefficient(igraph_t &G, 
								 igraph_integer_t E, 
								 igraph_real_t &clustering_coefficient, 
								 igraph_real_t &clustering_coefficient_r){
    
    igraph_transitivity_undirected(&G, &clustering_coefficient, IGRAPH_TRANSITIVITY_NAN); 

    igraph_t random_G;
    igraph_copy(&random_G, &G);
    igraph_rewire(&random_G, 2*E, IGRAPH_REWIRING_SIMPLE);
    igraph_transitivity_undirected(&random_G, &clustering_coefficient_r, IGRAPH_TRANSITIVITY_NAN); 

    igraph_destroy(&random_G);

    return 0;
}