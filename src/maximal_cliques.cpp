#include "maximal_cliques.h"

int maximal_cliques(igraph_t& G,
					int minimum_cliquesize,
                    igraph_integer_t& clique_count,
                    igraph_integer_t& clique_number){

    igraph_maximal_cliques_count(&G, &clique_count, minimum_cliquesize, 0);
    //igraph_clique_number(&G, &clique_number);

    return 0;
}