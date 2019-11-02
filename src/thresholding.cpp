// made on dev
// Carissa BLeker
// cbleker@vols.utk.edu

#include <igraph/igraph.h>

#include <vector>     // std::vector
#include <iostream>   // std::cout, std::cerr, std::endl
#include <fstream>    // fopen, fclose (to read igraph)
#include <algorithm>  // std::nth_element, std::min_element, std::max_element
#include <math.h>     // pow, sqrt, fabs, M_PI
#include <getopt.h>   // commandline argument parsing
#include <stdlib.h>   // atoi, atof
#include <sstream>    // stringstream
#include <cmath>      // std::copysign

#include "utils.h"
#include "math_ext.h"
#include "igraph_ext.h"
#include "spline.h"


#include "specialfunctions.h"  // ALGLIB header that contains ChiSquareCDistribution

//////////////////////////////////////////////////////////////////////////////
//     Thresholding functions                                                 //
///////////////////////////////////////////////////////////////////////////////
// TODO: DRY

std::string thresholdSpectral(igraph_t &G,
                     double l=0.5,
                     double u=0.99,
                     double increment=0.01,
                     int windowsize=5,
                     int minimumpartitionsize=10){
    
    // compare window size to minimumpartionsize
    if(minimumpartitionsize <= windowsize){
        std::cout << " Warning: cannot have minimumpartitionsize <= windowsize. Using windowsize = 5 and minimumpartitionsize = 10." << std::endl;
        minimumpartitionsize = 10;
        windowsize = 5;
    }

    // initialise necessary stuff
    igraph_integer_t E;         // number edges before threshold
    igraph_integer_t new_E;     // number edges after threshold
    igraph_integer_t V;         // number vertices
    igraph_integer_t V_cc;      // number vertices in LCC

    // other stuff
    igraph_real_t eigenvalue;
    igraph_vector_t eigenvector;
    igraph_vector_init(&eigenvector, 0);
    std::vector<double> window_differences;

    // nested loop
    double tol;
    int number_clusters;
    int cluster_begin;
    int cluster_end;
    bool in_step;
    double d;
    
    // get the threshold increments
    double t;
    static const std::vector<double> t_vector = range(l, u, increment);
    int num_increments = t_vector.size();
    std::cout << " Number steps: " << num_increments << std::endl;

    // results go here
    std::vector<int>  stat_per_t(num_increments);
    std::vector<double> second_eigenvalue_per_t(num_increments);

    // keep track of which thresholds were tested
    std::vector<bool> was_tested_per_t(num_increments, false);

    E = igraph_ecount(&G);

    for(int i_t=0; i_t < num_increments; i_t++){
        t = t_vector[i_t];

        std::cout << "\nStep: " << i_t << ", Threshold: " << t << std::flush;

        // Threshold step
        threshold_graph(t, G); 
        
        // make sure graph is large enough to continue
        V = igraph_vcount(&G);
        new_E = igraph_ecount(&G); 

        std::cout << " " << V << " " << new_E; 
        if(V < minimumpartitionsize){ //not large enough 
            std::cout <<" Graph too small, finished. " << std::flush;
            break;
        } 

        if(new_E < E){
            E = new_E;
        }
        else{
            std::cout << " New number edges is not less than previous number of edges, skipping. " << std::flush;
            continue;
        }

        // Get largest connected component
        igraph_t G_cc;
        largest_connected_component(G, G_cc);

		// make sure largest connected component is large enough to continue
        V_cc = igraph_vcount(&G_cc); 
        std::cout << " V_cc " << V_cc;
        if(V_cc < minimumpartitionsize){
            std::cout << " LCC too small, finished. " << std::flush;
            break;
        }

        Fiedler_vector(G_cc, eigenvector, eigenvalue);

        // destroy G_cc
        igraph_destroy(&G_cc);

        // keep eigenvalue of interest
        second_eigenvalue_per_t[i_t] = eigenvalue;

        // do the sort and step thing with the eigenvector
        igraph_vector_sort(&eigenvector);
		rolling_difference_igraph(eigenvector, window_differences, windowsize);

        tol = mean(window_differences) + stddev(window_differences)/2.0;
        std::cout << " tol: " << tol << std::flush; 

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
        stat_per_t[i_t] = number_clusters;
        was_tested_per_t[i_t] = true;
        std::cout << " Number clusters: " << number_clusters << std::flush;
    }

    std::cout << "\nDone\n" << std::endl;

    // make results into a string
    std::stringstream message;
    message << "threshold\tsecond eigenvalue\tnumber clusters\n";
    for(int i=0; i<stat_per_t.size(); i++){
        if(was_tested_per_t[i]){
            message << t_vector[i] << "\t" << second_eigenvalue_per_t[i] << "\t" << stat_per_t[i] << "\n";
        }
    }

    return message.str();
}

std::string thresholdCliqueDoubling(igraph_t &G,
                     double l=0.1,
                     double u=0.99,
                     double increment=0.01,
                     int minimumpartitionsize=3){

    // initialise necessary stuff
    igraph_integer_t E;            // number edges before threshold
    igraph_integer_t new_E;        // number edges after threshold
    igraph_integer_t V;            // number vertices
    igraph_integer_t clique_count; // number maximal cliques

    // get the threshold increments
    double t;
    static const std::vector<double> t_vector = range(l, u, increment);
    int num_increments = t_vector.size();
    std::cout << "Number steps: " << num_increments << std::endl;

    // results go here
    std::vector<double>  stat_per_t(num_increments); //ratios go here
    std::vector<int>    clique_count_per_t(num_increments);

    // keep track of which thresholds were tested
    std::vector<bool> was_tested_per_t(num_increments, false);

    E = igraph_ecount(&G);

    for(int i_t=0; i_t < num_increments; i_t++){
        t = t_vector[i_t];

        std::cout << "\nStep: " << i_t << ", Threshold: " << t << std::flush;

        // Threshold step
        threshold_graph(t, G); 
        
        // make sure graph is large enough to continue
        V = igraph_vcount(&G);
        new_E = igraph_ecount(&G); 

        if(V < minimumpartitionsize){ //not large enough 
            std::cout <<" Graph too small, finished. " << std::flush;
            break;
        } 

        if(new_E < E){
            E = new_E;
        }
        else{
            std::cout << " New number edges is not less than previous number of edges, skipping. " << std::flush;
            continue;
        }

        std::cout << " Calculating number of maximal cliques. " << std::flush;
        // number of maximal cliques
        igraph_maximal_cliques_count(&G, &clique_count, minimumpartitionsize, 0);

        clique_count_per_t[i_t] = clique_count;
        was_tested_per_t[i_t] = true;
    }

    igraph_destroy(&G);

    // ratio between thresholds: p_t = x_t/x_(t-1), no value for x_m
    int next_clique_num = clique_count_per_t[num_increments-1];
    int clique_num;

    for(int i = num_increments-2; i >= 0; i--){
        if(was_tested_per_t[i]){
            clique_num = clique_count_per_t[i];
            stat_per_t[i] = (double)clique_num / (double)next_clique_num;
            next_clique_num = clique_num;
        }
    }
    std::cout << "\nDone\n" << std::endl;

    // make results into a string
    std::stringstream message;
    message << "threshold\tnumber maximal cliques\tmaximal clique ratio\n";
    for(int i=0; i < num_increments; i++){
        if(was_tested_per_t[i]){
            message << t_vector[i] << "\t" << clique_count_per_t[i] << "\t" << stat_per_t[i] << "\n";
        }
    }

    return message.str();
}

// 
std::string thresholdPercolation(igraph_t &G,
                     double l=0.1,
                     double u=0.99,
                     double increment=0.01,
                     int minimumpartitionsize=3){
    // based on number of connected components and number of vertices 
    // remaining in the graph after threshold

    // initialise necessary stuff
    igraph_integer_t E;            // number edges before threshold
    igraph_integer_t new_E;        // number edges after threshold
    igraph_integer_t V;            // number vertices
    igraph_integer_t cc_count;     // number connected components

    // get the threshold increments
    double t;
    static const std::vector<double> t_vector = range(l, u, increment);
    int num_increments = t_vector.size();
    std::cout << "Number steps: " << num_increments << std::endl;

    // results go here
    std::vector<double>  v_per_t(num_increments); 
    std::vector<int>    cc_count_per_t(num_increments);

    // keep track of which thresholds were tested
    std::vector<bool> was_tested_per_t(num_increments, false);

    E = igraph_ecount(&G);

    for(int i_t=0; i_t < num_increments; i_t++){
        t = t_vector[i_t];

        std::cout << "\nStep: " << i_t << ", Threshold: " << t << std::flush;

        // Threshold step
        threshold_graph(t, G); 
        
        // make sure graph is large enough to continue
        V = igraph_vcount(&G);
        new_E = igraph_ecount(&G); 

        if(V < minimumpartitionsize){ //not large enough 
            std::cout <<" Graph too small, finished. " << std::flush;
            break;
        } 

        if(new_E < E){
            E = new_E;
        }
        else{
            std::cout << " New number edges is not less than previous number of edges, skipping. " << std::flush;
            continue;
        }

        v_per_t[i_t] = V;
        std::cout << " Calculating number of connected components. " << std::flush;
        igraph_clusters(&G, NULL, NULL, &cc_count, IGRAPH_STRONG);
        cc_count_per_t[i_t] = cc_count;
        was_tested_per_t[i_t] = true;
    }

    igraph_destroy(&G);

    std::cout << "\nDone\n" << std::endl;

    // make results into a string
    std::stringstream message;
    message << "threshold\tnumber connected components\tnumber vertices\n";
    for(int i=0; i < num_increments; i++){
        if(was_tested_per_t[i]){
            message << t_vector[i] << "\t" << cc_count_per_t[i] << "\t" << v_per_t[i] << "\n";
        }
    }

    return message.str();
}

// 
std::string thresholdRMT(igraph_t &G,
                     double l=0.1,
                     double u=0.99,
                     double increment=0.01,
                     int minimumpartitionsize=3){
    // RMT 

    // initialise necessary stuff
    igraph_integer_t E;            // number edges before threshold
    igraph_integer_t new_E;        // number edges after threshold
    igraph_integer_t V;            // number vertices


    // get the threshold increments
    double t;
    static const std::vector<double> t_vector = range(l, u, increment);
    int num_increments = t_vector.size();
    std::cout << "Number steps: " << num_increments << std::endl;

    // keep track of which thresholds were tested
    std::vector<bool> was_tested_per_t(num_increments, false);

    E = igraph_ecount(&G);


    // results go here
    std::vector<double> poi_chi_sq_stat_per_t(num_increments);;
    std::vector<double> goe_chi_sq_stat_per_t(num_increments);;

    std::vector<double> poi_chi_sq_pvalue_per_t(num_increments);;
    std::vector<double> goe_chi_sq_pvalue_per_t(num_increments);;

    for(int i_t=0; i_t < num_increments; i_t++){
        t = t_vector[i_t];

        std::cout << "\nStep: " << i_t << ", Threshold: " << t << std::flush;

        // Threshold step
        threshold_graph(t, G); 
        
        // make sure graph is large enough to continue
        V = igraph_vcount(&G);
        new_E = igraph_ecount(&G); 

        if(V < minimumpartitionsize){ //not large enough 
            std::cout <<" Graph too small, finished. " << std::flush;
            break;
        } 

        if(new_E < E){
            E = new_E;
        }
        else{
            std::cout << " New number edges is not less than previous number of edges, skipping. " << std::flush;
            continue;
        }

        // Here
        igraph_matrix_t A;
        igraph_matrix_init(&A, V, V);
        get_weighted_adjacency(G, A);
        
        igraph_vector_t eigenvalues;
        igraph_vector_init(&eigenvalues, V); // eigenvalues will go in here. 

        igraph_lapack_dsyevr(&A, IGRAPH_LAPACK_DSYEV_ALL, 0, 0, 0, 0, 0, 0.0, &eigenvalues, NULL, 0);

        // sort
        igraph_vector_sort(&eigenvalues);

        // make into std::vector and remove duplicates
        // scale eigenvalues: divide by largest absolute value eigenvalue
        // then add 1 (doesn't affect NNSD but fixes spline issues)
        double largest_eigenvalue = fabs(igraph_vector_tail(&eigenvalues));

        std::vector<double> eigenvalues_sorted;

        std::cout << "\n";
        double this_eigenvalue, previous_eigenvalue;

        this_eigenvalue = VECTOR(eigenvalues)[0] / largest_eigenvalue + 1;
        eigenvalues_sorted.push_back(this_eigenvalue);
        previous_eigenvalue = this_eigenvalue;
        for(int i=1; i<igraph_vector_size(&eigenvalues); i++){
            this_eigenvalue = VECTOR(eigenvalues)[i] / largest_eigenvalue + 1;
            if(fabs(previous_eigenvalue - this_eigenvalue) > 0.00001){
                eigenvalues_sorted.push_back(this_eigenvalue);
                std::cout << this_eigenvalue << "\t";
                previous_eigenvalue = this_eigenvalue;
            }
        }

        double n = eigenvalues_sorted.size();
     
        // CDF
        int number_fit_points = floor(0.75*eigenvalues_sorted.size()); //TODO
        double cdf_increment = (eigenvalues_sorted.back() - eigenvalues_sorted[0])/number_fit_points;
        std::vector<double> t = range(0, 2.1, cdf_increment);

        std::vector<double> cdf = ecdf(eigenvalues_sorted, t);

        // Find smooth distribution of eigenvalues by fitting a spline to the CDF and evaluating at eigenvalues values
        std::vector<double> new_cdf = spline(t, cdf, eigenvalues_sorted); 

        // NNSD 
        std::vector<double> NNSD;
        rolling_difference(new_cdf, NNSD, 1);

        int NNSD_size = NNSD.size();

        double NNSD_mean = mean(NNSD);

        std::cout << "\nNNSD_size " << NNSD_size << "\n";
        for (int i=0; i < NNSD_size; i++){
            NNSD[i] = NNSD[i] * n ;
            std::cout << NNSD[i] << "\t";

        }
        std::cout << "\n";

        // Chi2 test on NNSD
        // https://github.com/spficklin/RMTGeneNet/blob/master/threshold/methods/RMTThreshold.cpp
        // https://www.statisticshowto.datasciencecentral.com/goodness-of-fit-test/
        // https://static-content.springer.com/esm/art%3A10.1186%2F1471-2105-8-299/MediaObjects/12859_2006_1671_MOESM3_ESM.pdf
        
        std::sort(NNSD.begin(), NNSD.end());

        // Need to discretise continuous NNSD
        // Use range of [0, 3], and each bin must have 5 counts (observed values)
        double bin_start = 0;
        double bin_end = 0;

        std::vector<double> bin_start_vector;
        std::vector<double> bin_end_vector;
        std::vector<double> bin_count_vector;

        // observed bin frequency, will always be 5 (execpt for last bin)
        float observed_count = 5;
        double expected_count;
        
        int no_bins = 0;
        int NNSD_i = 0; // index of NNSD
        while(NNSD_i < NNSD_size){ // don't know how many bins yet, but going untill end of NNSD

            int this_bin_count = 0;
            bin_start = bin_end; //next bin starts at prev bin end

            int h;
            for(h = NNSD_i; h<NNSD_size; h++){
                if(this_bin_count == observed_count){
                    break;
                }
                else{
                    this_bin_count++;
                    }
            }
            
            if(this_bin_count == observed_count){
                // either still in the range of values, 
                // or reached end but divisivle by 5
                
                if(h == NNSD_size){
                    bin_end = NNSD[h-1];
                }
                else{
                    bin_end = NNSD[h];
                }

                bin_start_vector.push_back(bin_start);
                bin_end_vector.push_back(bin_end);
                bin_count_vector.push_back(this_bin_count);
                no_bins++;
                NNSD_i = h;
            }
            else{ 
                // reached end of last bins without making it until frequency of 5
                // so add this to the previous bin
                // bin start stays the same, bin_end is last value
                bin_end = NNSD[h-1];
                bin_end_vector[no_bins-1] = bin_end;
                bin_count_vector[no_bins-1] = bin_count_vector[no_bins-1] + this_bin_count;

                NNSD_i = h;
            }

        }
        std::cout << "\n\nNumber bins " << no_bins << std::endl;
        std::cout << "\n";
        
        double poi_chi_sq_stat = 0;
        double goe_chi_sq_stat = 0;

        for(int b=0; b<no_bins; b++){

            bin_start = bin_start_vector[b];
            bin_end = bin_end_vector[b];
            observed_count = bin_count_vector[b];

            std::cout << b << "\t" << bin_start << "\t " << bin_end << "\t " << observed_count << "\t ";

            // If Poisson, then 
            expected_count = NNSD_size * ( poisson(bin_start, bin_end) );
            poi_chi_sq_stat += pow(observed_count - expected_count, 2) / expected_count;
            std::cout << expected_count << "\t ";

            // If GOE, then
            expected_count = NNSD_size * ( goe(bin_start, bin_end) );
            goe_chi_sq_stat += pow(observed_count - expected_count, 2) / expected_count;
            std::cout << expected_count << std::endl;
        }

        poi_chi_sq_stat_per_t[i_t] = poi_chi_sq_stat;
        goe_chi_sq_stat_per_t[i_t] = goe_chi_sq_stat;

        // Chi2 test - 
        // df = no_bins - 1
        // p-value = area under right hand tail
        // http://www.alglib.net/specialfunctions/distributions/chisquare.php

        poi_chi_sq_pvalue_per_t[i_t] = alglib::chisquarecdistribution(no_bins -1, poi_chi_sq_stat);
        goe_chi_sq_pvalue_per_t[i_t] = alglib::chisquarecdistribution(no_bins -1, goe_chi_sq_stat);

        was_tested_per_t[i_t] = true;
    }

    igraph_destroy(&G);

    std::cout << "\nDone\n" << std::endl;

    // make results into a string
    std::stringstream message;
    message << "threshold\tPoisson Chi2\tPoisson p-value\tGOE Chi2\tGOE p-value\n";
    for(int i=0; i < num_increments; i++){
        if(was_tested_per_t[i]){
            message << t_vector[i] << "\t" << poi_chi_sq_stat_per_t[i] << "\t" << poi_chi_sq_pvalue_per_t[i];
            message                 << "\t" << goe_chi_sq_stat_per_t[i] << "\t" << goe_chi_sq_pvalue_per_t[i] << "\n";
        }
    }

    return message.str();
}


std::string thresholdDensity(igraph_t &G,
                     double l=0.1,
                     double u=0.99,
                     double increment=0.01,
                     int minimumpartitionsize=3){

    // initialise necessary stuff
    igraph_integer_t E;            // number edges before threshold
    igraph_integer_t new_E;        // number edges after threshold
    igraph_integer_t V;            // number vertices

    double density;

    // get the threshold increments
    double t;
    static const std::vector<double> t_vector = range(l, u, increment);
    int num_increments = t_vector.size();
    std::cout << "Number steps: " << num_increments << std::endl;

    // results go here
    std::vector<double>  stat_per_t(num_increments); // density goes here

    // keep track of which thresholds were tested
    std::vector<bool> was_tested_per_t(num_increments, false);

    E = igraph_ecount(&G);

    for(int i_t=0; i_t < num_increments; i_t++){
        t = t_vector[i_t];

        std::cout << "\nStep: " << i_t << ", Threshold: " << t << std::flush;

        // Threshold step
        threshold_graph(t, G); 
        
        // make sure graph is large enough to continue
        V = igraph_vcount(&G);
        new_E = igraph_ecount(&G); 

        if(new_E < E){
            E = new_E;
        }
        else{
            std::cout << " New number edges is not less than previous number of edges, skipping. " << std::flush;
            continue;
        }

        if(V < minimumpartitionsize){ //not large enough 
            std::cout <<" Graph too small, finished. " << std::flush;
            break;
        } 

        std::cout << " Calculating density. " << std::flush;
        density = 2.0 * E / (V * (V -1));
        stat_per_t[i_t] = density;

        was_tested_per_t[i_t] = true;
    }
    std::cout << "\nDone\n" << std::endl;

    igraph_destroy(&G);

    // make results into a string
    std::stringstream message;
    message << "threshold\tdensity\n";
    for(int i=0; i < num_increments; i++){
        if(was_tested_per_t[i]){
            message << t_vector[i] << "\t" <<  stat_per_t[i] << "\n";
        }
    }

    return message.str();
}


///////////////////////////////////////////////////////////////////////////////
//     Commandline arguments                                                 //
///////////////////////////////////////////////////////////////////////////////

void help(std::string prog_name){
    std::cerr <<  "\n Usage: \n";
    std::cerr <<  "   " << prog_name     << " [-OPTIONS]... <GRAPH FILE PATH>\n";
    std::cerr <<  "\n Graph has to be in .ncol format. \n";
    std::cerr <<  "\n Options: \n";
    std::cerr <<  "  -o  --out                      <filename>         path to store results\n";
    std::cerr <<  "                                                         if not given, results are sent to stdout\n";
    std::cerr <<  "  -l  --lower                    <value>            lower bound on thresholds to test (default 0.5)\n";
    std::cerr <<  "  -u  --upper                    <value>            upper bound on thresholds to test (default 0.99)\n";
    std::cerr <<  "  -i  --increment                <value>            threshold increment (default 0.01)\n";
    std::cerr <<  "  -w  --windowsize               <value>            sliding window size for spectral method (default 5)\n";
    std::cerr <<  "  -p  --minimumpartitionsize     <value>            minimum size of graph or subgraph after thresholding (default 10)\n";
    std::cerr <<  "  -m  --method                   [1|2|3|4]          method  (default = 1)\n";
    std::cerr <<  "                                                         1 - Spectral method\n";
    std::cerr <<  "                                                         2 - Clique ratio (generalised clique doubling)\n";
    std::cerr <<  "                                                         3 - Density\n";
    std::cerr <<  "                                                         4 - Percolation\n";
    std::cerr <<  "                                                         5 - Random matrix theory (not implemented)\n";
    std::cerr <<  "  -h  --help                                        print this help and exit\n";
    std::cerr <<  "\n";
    exit(1);
}

int arguement_parser(int argc, char **argv, 
        // Mandatory arguement definitions
        std::string &infile,  //input file name

        // Here flags (options without arguments) and arguments with defined type
        double &lower,
        double &upper,
        double &increment,
        int &windowsize,
        int &minimumpartitionsize,
        int &method,
        std::string &outfile
        ){

    int next_option;

    const char* const short_options = "hl:u:i:w:p:m:o:" ;
    const struct option long_options[] =
        {    //name,                    has_arg,    flag,        val 
            { "help",                   0,          NULL,        'h'},
            { "lower",                  1,          NULL,        'l'},
            { "upper",                  1,          NULL,        'u'}, 
            { "increment",              1,          NULL,        'i'}, 
            { "windowsize",             1,          NULL,        'w'},
            { "minimumpartitionsize",   1,          NULL,        'p'},
            { "method",                 1,          NULL,        'm'},
            { "outfile",                1,          NULL,        'o'},
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
    double l=0.5;
    double u=0.99;
    double increment=0.01;
    int windowsize=5;
    int minimumpartitionsize=10;
    int method=1;
    std::string outfile_name;

    arguement_parser(argc, argv, infile, l, u, increment, windowsize, minimumpartitionsize, method, outfile_name);

    //std::cout << "infile\t" << infile << "\n";
    //std::cout << "lower\t" << lower << "\n";
    //std::cout << "upper\t" << upper << "\n";
    //std::cout << "increment\t" << increment << "\n";
    //std::cout << "method\t" << method << "\n";
    //std::cout << "outfile\t" << outfile << "\n";

    // check that thresholding range is good
    if(l>=u){
        std::cout << "Error in threshold limits: cannot have l >= u" << std::endl;
        return 0;
    }
    
    // check method num is in range
    // TODO

    // turn on attribute handling
    // for igraph to handle edge weights
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    
    // Load graph
    igraph_t G;
    std::cout << "Loading graph ... " << std::flush;
    read_graph(infile, G);
    std::cout << "done." << std::endl;

    std::string message;

    switch(method){
        case 1 : message = thresholdSpectral(G,
                           l=l, 
                           u=u,
                           increment=increment,
                           windowsize=windowsize,
                           minimumpartitionsize=minimumpartitionsize);
                 break;  
        case 2 : message = thresholdCliqueDoubling(G,
                           l=l,
                           u=u,
                           increment=increment,
                           minimumpartitionsize=minimumpartitionsize);
                 break;
        case 3 : message = thresholdDensity(G, 
                           l=l, 
                           u=u,
                           increment=increment,
                           minimumpartitionsize=minimumpartitionsize);
                 break;
        case 4 : message = thresholdPercolation(G, 
                           l=l, 
                           u=u,
                           increment=increment,
                           minimumpartitionsize=minimumpartitionsize);
                 break;
        case 5 : message = thresholdRMT(G, 
                           l=l, 
                           u=u,
                           increment=increment,
                           minimumpartitionsize=minimumpartitionsize);
                 break;

        default : message =  "";
                       std::cout << "\nNot a valid method selected.\n";
    }

    output_results(outfile_name, message);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
