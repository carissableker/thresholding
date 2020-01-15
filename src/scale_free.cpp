#include "scale_free.h"

int scale_free_test(igraph_t &G,
					igraph_integer_t V,
					igraph_plfit_result_t &scale_free_result){

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, V); // degrees will go in here. 

    igraph_degree(&G, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

    igraph_power_law_fit(&degrees, &scale_free_result, -1, 0);

    igraph_vector_destroy(&degrees);

    return 0;
}