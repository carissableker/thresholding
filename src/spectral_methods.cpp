#include "spectral_methods.h"

int spectral_methods(igraph_t& G_cc,
                     int windowsize,
                     int minimumpartitionsize,
					 igraph_real_t& second_eigenvalue,
					 int& number_clusters){

	number_clusters = 1; // first

    igraph_vector_t eigenvector;
    igraph_vector_init(&eigenvector, 0);
    std::vector<double> window_differences;

    Fiedler_vector(G_cc, eigenvector, second_eigenvalue);

    // do the sort and step thing with the eigenvector
    igraph_vector_sort(&eigenvector);

	rolling_difference_igraph(eigenvector, window_differences, windowsize);

    double tol = mean(window_differences) + stddev(window_differences)/2.0;
    int cluster_begin = 0;
    int cluster_end = 0;
    bool in_step = false;

    for(int i=0; i<window_differences.size(); i++){

        double d = window_differences[i];

        if(d >= tol){
            // need to enter or stay in a step
            if(in_step == false){
                // enter step and end a cluster
                in_step = true;
                cluster_end = i;
                // end the last cluster,
                // add it to the number of clusters if it is large enough
                if(cluster_end - cluster_begin >= minimumpartitionsize){
                    number_clusters = number_clusters+1;
                }
            }
            // else we're already in the step, so do nothing
        }
        else{
            // not in a step, we're entering or still in a cluster
            if (in_step == true){
                //  entering a cluster
                in_step = false;
                cluster_begin = i;
            }
            //  else already in a cluster so else nothing
        }

    }

    return 0;
}
