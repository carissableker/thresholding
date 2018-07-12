
#include <igraph.h>

#include <vector>     // for std::vector
#include <iostream>   // for std::cout, std::cerr, std::endl
#include <fstream>    // for fopen, fclose (to read igraph)
#include <algorithm>  // for std::nth_element
#include <math.h>     // for pow, sqrt
#include <getopt.h>   // commandline argument parsing
#include <stdlib.h>   // atoi, atof

///////////////////////////////////////////////////////////////////////////////
//     Utility functions                                                     //
///////////////////////////////////////////////////////////////////////////////

// Read in graph
int read_graph(std::string &graph_file_path, igraph_t& G){

    FILE *graph_file;
    graph_file = fopen(graph_file_path.c_str(), "r"); 
    
    // Test if file exists
    if (graph_file == NULL){
        std::cerr << "Error - Unable to open file: " << graph_file_path << std::endl;
        exit(-1);
    }

    // Read in file as graph
    igraph_read_graph_ncol(&G, graph_file, NULL, false, IGRAPH_ADD_WEIGHTS_YES, IGRAPH_UNDIRECTED);

    fclose(graph_file);
    
    return 0;
}

// Threshold graph
// by removing edges with abs weight less than "t" (implemented)
// and subsequently vertices with no neightbours (not implemented)
int threshold_graph(float t, igraph_t &G){
	
    // identify edges to remove
    igraph_vector_t weights;
    igraph_vector_init(&weights, 0);

    igraph_vector_t edge_indices;
    igraph_vector_init(&edge_indices, 0);

    // should this be inside the loop, ie per edge???
    igraph_cattribute_EANV(&G, "weight", igraph_ess_all(IGRAPH_EDGEORDER_ID), &weights);

	igraph_real_t w;
    for (int i=0; i<igraph_ecount(&G); i++){
        w = VECTOR(weights)[i];
        if(w < t &&  w > -t){
            // add this edge index to list of edges to be deleted
            igraph_vector_push_back(&edge_indices, i);
        }     
    }

	std::cout << " About to remove " << igraph_vector_size(&edge_indices) << " edges. ";
    // remove edges 
    igraph_delete_edges(&G, igraph_ess_vector(&edge_indices));
	std::cout << "Removed edges\n";

	igraph_vector_destroy(&weights);
	igraph_vector_destroy(&edge_indices);

    return 0;
}

// Identify largest connected component of the graph and induce
int largest_connected_component(igraph_t &G, igraph_t &G_cc){
    // See also igraph_decompose, but since we only need
    // the largest CC, there is no point inducing all CCs.

	std::cout << " In cc code ";
    igraph_vector_t membership;
    igraph_vector_init(&membership, 0);
    igraph_vector_t csize;
    igraph_vector_init(&csize, 0);
    igraph_integer_t no; 

    igraph_clusters(&G, &membership, &csize, &no, IGRAPH_STRONG);
        
	std::cout << "clusters ";

    // iterate over csize to find largest CC
    int max_cc_size = 0;
    int max_cc_index;
    int this_cc_size;
    for(int i =0; i<no; i++){
        this_cc_size = VECTOR(csize)[i];
        if(this_cc_size > max_cc_size){
            max_cc_index = i;
            max_cc_size = this_cc_size;
        }
     }
    // identified largest CC, now collect its vertices
    igraph_vector_t vertices_in_cc;
    igraph_vector_init(&vertices_in_cc, max_cc_size);
    int j = 0; // index in vertices_in_cc

    for(int i=0; j<max_cc_size; i++){
        if(VECTOR(membership)[i] == max_cc_index){
            VECTOR(vertices_in_cc)[j] = i;
            j += 1;
        }
    }
	
	std::cout << " vertices list";

	std::cout << " number vertices " << igraph_vector_size(&vertices_in_cc);

    // induce the subgraph of the largest CC
    igraph_induced_subgraph(&G, &G_cc, igraph_vss_vector(&vertices_in_cc), IGRAPH_SUBGRAPH_AUTO);
	std::cout << "induced Gcc";

	igraph_vector_destroy(&membership);
	igraph_vector_destroy(&csize);
	igraph_vector_destroy(&vertices_in_cc);

    return 0;
}

// Returns the vector of differences between first and 
// last elements of the windows of size n in x
// from igraph_vector_t to  std::vector
int rolling_difference(igraph_vector_t &x, std::vector<float> &out, int n){

    int len_out = igraph_vector_size(&x) - n + 1;
    out.resize(len_out);

    // for each window start
    for(int ind=0; ind < len_out; ind++) {
        // difference is between first and last elements
        // of the window
        out[ind] = VECTOR(x)[ind+n-1] - VECTOR(x)[ind];
    }
    return 0;
}

// Median of a vector
float median(std::vector<float> v){
    
    int size = v.size();
    if(size == 0){
        return NAN;
    }

    else if(size == 1){
        return v[0];
    }

    else{
        if(size % 2 == 0){
            // even length vector, average of middle 2 numbers
            float o0;
            float o1;
            std::nth_element(v.begin(), v.begin() + size/2, v.end());
            o1 = v[size/2];
            o0 = *std::max_element(v.begin(), v.begin() + size/2 -1);
            return (o0 + o1)/2;
        } 
        else{
            // odd length, middle number
            std::nth_element(v.begin(), v.begin() + (size-1)/2, v.end());
            return v[(size-1)/2];
        }
    }
}

// Mean/average of a vector
float mean(std::vector<float> v){

	float mean = 0;
    float n = v.size();
    for(int i=0; i<n; i++){
        mean = mean + v[i];
    }
    return mean/n;
}

// Standard deviation of a vector
float stdev(std::vector<float> v){
    float n = v.size();
    float v_bar = mean(v);
    float var = 0;

    for(int i=0; i<n; i++){
        var = var + pow( (v[i] - v_bar), 2.0);
    }
    var = var / (n-1);

    return sqrt(var);
}

// Range from l to u, incrementing by increment
std::vector<float> range(float l, float u, float increment){
	// Floating point arithmetic
	float precision = 100.0f;
	int int_increment = static_cast<int>(increment*precision);
	int int_l = static_cast<int>(l*precision);
	int int_u = static_cast<int>(u*precision);

	//std::cout << int_increment << "\t" << int_l << "\t" << int_u << "\n";
	std::vector<float> out;

	for(int int_t=int_l; int_t<=int_u; int_t+=int_increment){
		out.push_back(int_t/precision);
	}
	
	return out;
}


///////////////////////////////////////////////////////////////////////////////
//     Thresholding functions                                                 //
///////////////////////////////////////////////////////////////////////////////

int thresholdSpectral(igraph_t &G_p,
                      float l=0.5,
                      float u=0.99,
                      float increment=0.01,
                      int windowsize=5,
                      int minimumpartitionsize=10){
    
    // results go here
    std::vector<int>  stat_per_t;

    // initialise necessary stuff
    igraph_integer_t E;
    igraph_matrix_t laplacian;
    igraph_matrix_init(&laplacian, 1, 1);
    igraph_vector_t values;
    igraph_vector_init(&values, 0);
    igraph_matrix_t vectors;
    igraph_matrix_init(&vectors, 0, 1);
	igraph_real_t eigen_value;
    igraph_vector_t eigen_vector;
    igraph_vector_init(&eigen_vector, igraph_matrix_nrow(&vectors));


	// connected component
	igraph_t G_cc;
		

    std::vector<float> window_differences(10);

    // nested loop
    float tol;
    int number_clusters;
    int cluster_begin;
    int cluster_end;
    bool in_step;
    float d;
    
	float t;
	std::vector<float> vector_t = range(l, u, increment);
	int num_increments = vector_t.size();
	std::cout << "Number steps: " << num_increments << "\n";

	for(int ii=0;ii<num_increments;ii++){
		std::cout << vector_t[ii] << ", ";
	}
	std::cout << "\n";

    for(int i_t=0; i_t < num_increments; i_t++){
        t = vector_t[i_t];

        std::cout << "Step: " << i_t << ", Threshold: " << t;

        // Threshold step
        threshold_graph(t, G_p); 

        // make sure graph is large enough to continue
        E = igraph_ecount(&G_p); 
		std::cout <<  ", Number edges: " << E;

        if(E < minimumpartitionsize){break;} //not large enough

        // Select the largest connected component (CC)
        largest_connected_component(G_p, G_cc);

		std::cout << " Got G_cc\t";

        if(igraph_ecount(&G_cc) < minimumpartitionsize){break;}

        // Laplacian
        igraph_laplacian(&G_cc, &laplacian, NULL, 0, NULL);
       
		// Eigen decomposition for symmetric matrices using LAPACK
        igraph_lapack_dsyevr(&laplacian, IGRAPH_LAPACK_DSYEV_SELECT, 0, 0, 0, 2, 2, 1e-10, &values, &vectors, 0); 

		// destroy G_cc
		igraph_destroy(&G_cc);

		// output eigen value of interest
		eigen_value = VECTOR(values)[0];
		std::cout << ", 2nd Eigenvalue: " << eigen_value;

        igraph_matrix_get_col(&vectors, &eigen_vector, 0);

        igraph_vector_sort(&eigen_vector);

        // do the step thing with the eigenvector
        rolling_difference(eigen_vector, window_differences, windowsize);

        tol = mean(window_differences) + stdev(window_differences)/2.0;
		
		std::cout << ", TOL: " << tol;

        number_clusters = 1;
        cluster_begin = 0;
        cluster_end = 0;
        in_step = false;

        for(int i=0; i<window_differences.size(); i++){
            d = window_differences[i];

            if(d >= tol){
                // need to enter or stay in a step
                if(in_step == false){
                    // enter step and end a cluster
                    in_step = true;
                    cluster_end = i;
                    // end the last cluster, add it to the number of clusters if it is large enough
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
		std::cout << ", Number clusters: " << number_clusters << "\n";
        stat_per_t.push_back(number_clusters);
    }

	std::cout << "\nDone\n";

    for(int i=0; i<stat_per_t.size(); i++){
        std::cout << vector_t[i] << "\t" << stat_per_t[i] << "\n";
    }

    igraph_destroy(&G_p);
    return 0;
}

int thresholdCliqueDoubling(igraph_t &G_p,
                            float l=0.1,
                            float u=0.99,
                            float increment=0.01,
                            int minimumpartitionsize=3){

    //igraph_t G_p;
    //igraph_copy(&G_p, &G);

    // results go here
    std::vector<float> stat_per_t;

    // initialise necessary stuff
    std::vector<int> clique_count_per_t;

    igraph_integer_t E;
    igraph_integer_t res;

	float t;

	std::vector<float> vector_t = range(l, u, increment);
	int num_increments = vector_t.size();
	std::cout << "Number steps: " << num_increments << "\n";


    for(int i_t=0; i_t < num_increments; i_t++){
        t = vector_t[i_t];

        // Threshold step
        threshold_graph(t, G_p); 

        // make sure graph is large enough to continue
        E = igraph_ecount(&G_p); //std::cout << "Threshold " << t << ", Number edges " << E << "\n";
        
        if(E < minimumpartitionsize){break;} //not large enough
        
        // number of maximal cliques
        igraph_maximal_cliques_count(&G_p, &res, minimumpartitionsize, 0);

        clique_count_per_t.push_back(res);
    }

    // propprtianal difference between thresholds
    stat_per_t.resize(vector_t.size(), 0);

    for(int i = 1; i<vector_t.size(); i++){
        if(clique_count_per_t[i] != 0){
            stat_per_t[i] = (float)clique_count_per_t[i-1] / clique_count_per_t[i]; 
        }
    }

    // print results
    for(int i=0; i<vector_t.size(); i++){
        std::cout << vector_t[i] << "\t" << stat_per_t[i] << std::endl;
    }

    igraph_destroy(&G_p);
    return 0;
}


int thresholdRandomMatrixTheory(){
}

int thresholdPercolation(){
}


///////////////////////////////////////////////////////////////////////////////
//     Commandline arguments                                                 //
///////////////////////////////////////////////////////////////////////////////

void help(std::string prog_name){
    std::cerr <<  "\n Usage: \n";
    std::cerr <<  "   " << prog_name     << " [-OPTIONS]... <GRAPH FILENAME>\n";
    std::cerr <<  "\n Options: \n";
    std::cerr <<  "  -o  --out                      <filename>        path to store results (not implemented)\n";
    std::cerr <<  "                                                         if not given, results are sent to stdout\n";
    std::cerr <<  "  -l  --lower                    <value>            lower bound on thresholds to test (default 0.5)\n";
    std::cerr <<  "  -u  --upper                    <value>            upper bound on thresholds to test (default 0.99)\n";
    std::cerr <<  "  -i  --increment                <value>            threshold increment (default 0.01)\n";
    std::cerr <<  "  -w  --windowsize               <value>            sliding window size for spectral method (default 5)\n";
    std::cerr <<  "  -p  --minimumpartitionsize     <value>            minimum size of graph or subgraph after thresholding (default 3)\n";
    std::cerr <<  "  -m  --method                   [1|2|3]            method  (default = 1)\n";
    std::cerr <<  "  -h  --help                                        Print this help and exit\n";
	std::cerr <<  "\n";
    exit(1);
}

int arguement_parser(int argc, char **argv, 
        // Mandatory arguement definitions
        std::string &infile,  //input file name

        // Here flags (options without arguments) and arguments with defined type
        float &lower,
        float &upper,
        float &increment,
        int &windowsize,
        int &minimumpartitionsize,
        int &method,
        std::string &outfile
        ){

    int next_option;

    const char* const short_options = "hl:u:i:w:p:m:o:" ;
    const struct option long_options[] =
        {    //name,                    has_arg,    flag,        val 
            { "help",                   0,          NULL,        'h' },
            { "lower",                  1,          NULL,        'l' },
            { "upper",                  1,          NULL,        'u'}, 
            { "increment",              1,          NULL,        'i'}, 
            { "windowsize",             1,          NULL,        'w'},
            { "minimumpartitionsize",   1,          NULL,        'p'},
            { "method",                 1,          NULL,        'm'},
            { "outfile",                1,          NULL,        'o' },
            { NULL, 0, NULL, 0 }
        };
 
    // Parse options
    while (1) {
        // Obtain an option
        next_option = getopt_long(argc, argv, short_options, long_options, NULL);
 
        if (next_option == -1)
            break; // No more options. Break loop.
 
        switch (next_option){
 
            case 'h' : // -h or --help 
                help(argv[0]);
                break;
 
            case 'l' : // -l or --lower
                lower=atof(optarg);
                break;

             case 'u' : // -u or --upper
                upper=atof(optarg);
                break;

             case 'i' : // -i or --increment
                increment=atof(optarg);
                break;
            
            case 'w': // -w or --windowsize
                windowsize=atoi(optarg);
                break;

             case 'm' : // -m or --method
                method=atoi(optarg);
                break;
 
            case 'p': // -p or --minimumpartitionsize
                minimumpartitionsize=atoi(optarg);
                break;

            case 'o' : // -o or --outfile
                outfile=optarg;
                break;
 
            case '?' : // Invalid option
                help(argv[0]); // Return help
 
            case -1 : // No more options
                break;
 
            default : // Something unexpected? Aborting
                return(1);
        }
    }
 
    // Mandatory arguements 
    // Current index (optind) must be smaller than the total number of arguments
    if(optind == argc){
        std::cerr << "\n Mandatory argument(s) missing\n";
        help(argv[0]);
    }
    // Iterate over rest of the arguments (i.e. in argv[optind])
    while (optind < argc){
        // only mandatory arguement at this stage is input file name
        infile = argv[optind];        
        optind++;
    }
 
    return 0;
}


///////////////////////////////////////////////////////////////////////////////
//     Main                                                                  //
///////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv){
    // Parse arguements

    // Mandatory arguement definitions
    std::string infile;  //input file name
   // Flags (options without arguments) and arguments with defined type
    float l=0.5;
    float u=0.99;
    float increment=0.01;
    int windowsize=5;
    int minimumpartitionsize=3;
    int method=1;
    std::string outfile_name;

    arguement_parser(argc, argv, infile, l, u, increment, windowsize, minimumpartitionsize, method, outfile_name);

    //std::cout << "infile\t" << infile << "\n";
    //std::cout << "lower\t" << lower << "\n";
    //std::cout << "upper\t" << upper << "\n";
    //std::cout << "increment\t" << increment << "\n";
    //std::cout << "method\t" << method << "\n";
    //std::cout << "outfile\t" << outfile << "\n";
    

	// write results to outfile if it exists, 
	// otherwise write results to std::cout
	std::ostream* outputTarget = &std::cout;

	if(!outfile_name.empty()){
		std::ofstream outfile;
    		outfile.open(outfile_name.c_str());
		outputTarget = &outfile;
	}

    // turn on attribute handling
	// for igraph to handle edge weights
    igraph_i_set_attribute_table(&igraph_cattribute_table);
       
    igraph_t G;
    
	std::cout << "Loading graph ... ";
    read_graph(infile, G);
	std::cout << "done.\n";

    if(method==1){
        thresholdSpectral(G,
                          l=l, 
                          u=u,
                          increment=increment,
                          windowsize=windowsize,
                          minimumpartitionsize=minimumpartitionsize);
    }
    else if(method==2){
        thresholdCliqueDoubling(G,
                                l=l,
                                u=u,
                                increment=increment,
                                minimumpartitionsize=minimumpartitionsize);
    }
    else{
        std::cout << "nothing\n";
    }

    return 0;
}


///////////////////////////////////////////////////////////////////////////////
