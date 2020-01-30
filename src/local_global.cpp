#include "local_global.h"

int local_global_pruning(igraph_t& G,
				 double alpha,
				 igraph_t& new_G,
				 double& mean_k
				 ){

    igraph_integer_t E = igraph_ecount(&G);
    igraph_integer_t V = igraph_vcount(&G);


    igraph_vector_t edge_weights;
    igraph_vector_init(&edge_weights, E);

	igraph_vit_t v_iterator;
	igraph_vit_create(&G, igraph_vss_all(), &v_iterator);
	igraph_integer_t v_id;

	igraph_integer_t e_id;
  	igraph_real_t w;

  	mean_k = 0;

	while (!IGRAPH_VIT_END(v_iterator)) {

		igraph_es_t incident_edges; //selector
		igraph_eit_t e_iterator;

		igraph_integer_t num_neighbours;

		double k;

		v_id = IGRAPH_VIT_GET(v_iterator);

	  	igraph_es_incident(&incident_edges, v_id, IGRAPH_ALL);
	  	igraph_es_size(&G, &incident_edges, &num_neighbours);
	  	igraph_eit_create(&G, incident_edges, &e_iterator);

	  	std::vector<double> incident_weights(num_neighbours);


	  	// get all weights of incident edges
	  	int i=0; // neighbour index for incident_weights
	  	while (!IGRAPH_EIT_END(e_iterator)){
	  		e_id = IGRAPH_EIT_GET(e_iterator);
	  		w = igraph_cattribute_EAN(&G, "weight", e_id);
	  		incident_weights[i] = w;
	  		IGRAPH_EIT_NEXT(e_iterator);
	  		i++;
	  	}

	  	// local threshold for v_id
	  	k = mean(incident_weights) + alpha * stddev(incident_weights);
	  	mean_k += k;

	  	// test each incident edge against local threshold
	  	// add it to the global vector of all edge weights
	  	IGRAPH_EIT_RESET(e_iterator);
	  	while (!IGRAPH_EIT_END(e_iterator)){
	  		e_id = IGRAPH_EIT_GET(e_iterator);
	  		w = igraph_cattribute_EAN(&G, "weight", e_id);
	  		if (w > k){
	  			// add 0.5 to new weights
	  			VECTOR(edge_weights)[e_id] += 0.5;
	  		}
	  		IGRAPH_EIT_NEXT(e_iterator);
	  	}

	  IGRAPH_VIT_NEXT(v_iterator);
	}
	igraph_vit_destroy(&v_iterator);

	mean_k /= V;


	// select edges that passed the local threshold
	// iterate over all edges
	igraph_eit_t all_e_iterator;
  	igraph_eit_create(&G, igraph_ess_all(IGRAPH_EDGEORDER_ID), &all_e_iterator);

    igraph_vector_t edge_indices;
    igraph_vector_init(&edge_indices, 0);

  	while (!IGRAPH_EIT_END(all_e_iterator)){
  		e_id = IGRAPH_EIT_GET(all_e_iterator);
  		w = VECTOR(edge_weights)[e_id];
  		if (w > 0){
			igraph_vector_push_back(&edge_indices, e_id);
		}
  		IGRAPH_EIT_NEXT(all_e_iterator);
  	}
  	igraph_eit_destroy(&all_e_iterator);

  	// induce graph with these edges, delete non-adjacent vertices=True
  	igraph_subgraph_edges(&G, &new_G, igraph_ess_vector(&edge_indices), 1);


	igraph_vector_destroy(&edge_weights);
	igraph_vector_destroy(&edge_indices);

	return 0;
}


int local_global_method(igraph_t& G,
				 double min_alpha,
				 double max_alpha,
				 double alpha_increment,
                 int windowsize,
                 int minimumpartitionsize,
                 std::string& outfile_name
				 ){

    double alpha;
    static const std::vector<double> alpha_vector = range(min_alpha, max_alpha, alpha_increment);
    int num_increments = alpha_vector.size();

	std::ofstream out;
    out.open(outfile_name.c_str());

    // output header
    std::stringstream header;
    header << "alpha";
    header << "\tvertex-count\tedge-count";
    header << "\tmean-k";
    header << "\tconnected-component-count";
    header << "\tlargest-cc-size\t2nd-largest-cc-size";
    header << "\t2nd-eigenvalue\talmost-disconnected-component-count";
    out << header.str();
    out << std::endl;

	igraph_t new_G;
    igraph_t G_cc;

    for (int i_alpha=0; i_alpha < num_increments; i_alpha++){
        alpha = alpha_vector[i_alpha];
        double mean_k;

	  	local_global_pruning(G, alpha, new_G, mean_k);

		igraph_integer_t E = igraph_ecount(&new_G);
    	igraph_integer_t V = igraph_vcount(&new_G);

    	igraph_integer_t cc_count           = -1;
		igraph_integer_t largest_cc_size    = -1;
		igraph_integer_t largest2_cc_size   = -1;

    	largest_connected_component(new_G, G_cc,
    		cc_count, largest_cc_size, largest2_cc_size);


    	int nearly_disconnected_components  = 0;
		igraph_real_t  second_eigenvalue    = std::nan("");

        spectral_methods(G_cc,
            windowsize,
            minimumpartitionsize,
            second_eigenvalue,
            nearly_disconnected_components);

        // message
        std::stringstream message;
        message << alpha;
        message << "\t" << V << "\t" << E;
        message << "\t" << mean_k;
        message << "\t" << cc_count;
        message << "\t" << largest_cc_size   << "\t" << largest2_cc_size;
        message << "\t" << second_eigenvalue << "\t" << nearly_disconnected_components;
        out << message.str();
        out << std::endl;
    }

    igraph_destroy(&new_G);
    igraph_destroy(&G_cc);
    out.close();

    return 0;
}
