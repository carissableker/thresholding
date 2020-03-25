// Carissa BLeker
// cbleker@vols.utk.edu

#include <igraph/igraph.h>

#include <vector>     // std::vector
#include <iostream>   // std::cout, std::cerr, std::endl
#include <fstream>    // fopen, fclose (to read igraph), ofstream
#include <algorithm>  // std::nth_element, std::min_element, std::max_element
#include <math.h>     // pow, sqrt, fabs, M_PI
#include <getopt.h>   // commandline argument parsing
#include <stdlib.h>   // atoi, atof
#include <sstream>    // stringstream
#include <cmath>      // std::copysign
#include <math.h>     // isnormal
#include <set>        // sets
#include <string>     // getline
#include <csignal>    // signal

#include "utils.h"
#include "math_ext.h"
#include "igraph_ext.h"

#include "random_matrix_theory.h"
#include "spectral_methods.h"
#include "scale_free.h"
#include "maximal_cliques.h"
#include "local_global.h"
#include "significance.h"
#include "local_rank.h"

int thresholdAnalysis(std::string& outfile_prefix,
                      igraph_t &G,
                      double l,
                      double u,
                      double increment,
                      int windowsize,
                      int minimumpartitionsize,
                      int minimum_cliquesize,
                      double min_alpha,
                      double max_alpha,
                      double alpha_increment,
                      int num_samples,
                      double significance_alpha,
                      double bonferroni_corrected,
                      std::set<int> methods){

    std::string outfile_name;
    igraph_integer_t E = igraph_ecount(&G); // number edges
    igraph_integer_t V = igraph_vcount(&G); // number vertices

    // Max number edges based on number of vertices
    double orig_max_E = 0.5 * V * (V - 1.0);

    std::cout << "Number vertices:  " << V << "\n";
    std::cout << "Number edges:     " << E;
    std::cout << "  (maximum possible number edges " << int(orig_max_E) << ")";
    std::cout << std::endl;
    std::cout << "------------------------------------------------\n\n";

    ///////////////////////////////////////////////////////////////////////
    // non loop methods
    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    // local-global (guzzi2014, rank)
    if(methods.find(2)!=methods.end()){
        outfile_name = outfile_prefix + "local_global.txt";
        local_global_method(G,
                     min_alpha,
                     max_alpha,
                     alpha_increment,
                     windowsize,
                     minimumpartitionsize,
                     outfile_name);
        methods.erase(2);
    }

    if (methods.size() == 0){
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////
    // Threshold loop
    // loop destroys the graph
    ///////////////////////////////////////////////////////////////////////

    // ready the output file
    outfile_name = outfile_prefix + "iterative.txt";
    std::ofstream out;
    out.open(outfile_name.c_str(), std::ofstream::out);
    // is it open
    if (out.fail()) {
        std::cerr << "Error opening file for writing: " << outfile_name << "\n";
        return 0;
    }



    // output header
    std::stringstream header;
    header << "threshold";
    header << "\tvertex-count\tedge-count";
    header << "\tconnected-component-count";
    header << "\tdensity\tdensity-orig-V";
    header << "\tlargest-cc-size\t2nd-largest-cc-size";
    header << "\tclustering-coefficient\trandom-clustering-coefficient";
    header << "\t2nd-eigenvalue\talmost-disconnected-component-count";
    header << "\tmaximal-clique-count\tclique-number";
    header << "\tpoisson-chi2\tpoisson-pvalue";
    header << "\tgoe-chi2\tgoe-pvalue";
    header << "\tscale-free-KS\tscale-free-KS-p-value\tscale-free-alpha";
    out << header.str();
    out << std::endl;

    // get the threshold increments
    double t;
    static const std::vector<double> t_vector = range(l, u, increment);
    int num_increments = t_vector.size();
    std::cout << "Iterative thresholding\n";
    std::cout << "Number steps: " << num_increments << std::endl;

    // Initialise necessary stuff
    int nearly_disconnected_components      = -1;
    igraph_real_t  second_eigenvalue        = std::nan("");

    igraph_integer_t    clique_count        = -1;   // number maximal cliques
    igraph_integer_t    clique_number       = -1;   // maximum clique size

    double  density                         = std::nan("");
    double  density_orig_V                  = std::nan("");

    double  poi_chi_sq_stat                 = std::nan("");
    double  goe_chi_sq_stat                 = std::nan("");

    double  poi_chi_sq_pvalue               = std::nan("");
    double  goe_chi_sq_pvalue               = std::nan("");

    igraph_integer_t     cc_count           = -1;
    igraph_integer_t     largest_cc_size    = -1;
    igraph_integer_t     largest2_cc_size   = -1;

    double  scale_free_pvalue               = std::nan("");
    double  scale_free_KS                   = std::nan("");
    double  scale_free_xmin                 = std::nan("");
    double  scale_free_alpha                = std::nan("");

    igraph_real_t clustering_coefficient    = std::nan("");
    igraph_real_t clustering_coefficient_r  = std::nan("");


    ///////////////////////////////////////////////////////////////////////
    // if methods includes clustering coefficient,
    // then need a copy of G to threshold
    // rewiring at each threshold takes longer (?)
    igraph_t G_random;
    if (methods.find(7) !=methods.end()){
        igraph_copy(&G_random, &G);

        // igraph rewire loses edge weights,
        // need to save them and reassign back ()
        igraph_rewire(&G_random, E, IGRAPH_REWIRING_SIMPLE);

        igraph_vector_t edge_weights;
        igraph_vector_init(&edge_weights, E);

        igraph_cattribute_EANV(&G, "weight", igraph_ess_all(IGRAPH_EDGEORDER_ID), &edge_weights);
        igraph_cattribute_EAN_setv(&G_random, "weight", &edge_weights);

        igraph_vector_destroy(&edge_weights);
    }


    for(int i_t=0; i_t < num_increments; i_t++){
        t = t_vector[i_t];

        std::cout << "\nStep: " << i_t +1  << ", Threshold: " << t << std::flush;

        // Threshold step
        int threshold_status = threshold_graph(t, G);
        V = igraph_vcount(&G);
        E = igraph_ecount(&G);


        if(threshold_status == 1){
            // 1 = all edges removed, stop
            std::cout <<" Graph is empty, finished. " << std::flush;
            break;
        }
        else if( threshold_status == 3){
            // 3 = no edges removed, only skip if not first iteration
            if(i_t > 0){
                std::cout << " No edges removed, skipping. " << std::flush;
                continue;
            }
            else{
                std::cout << "\t\tVertices: " << V << "\tEdges: " << E;
                std::cout << std::flush;
            }
        }
        else if(threshold_status == 2){
            // 2 = some edges are removed, keep going
            // make sure graph is large enough to continue
            if ( (V < minimumpartitionsize) || (E < minimumpartitionsize) ){
                //not large enough
                std::cout <<" Graph too small, finished. " << std::flush;
                break;
            }
            else{
                std::cout << "\t\tVertices: " << V << "\tEdges: " << E;
                std::cout << std::flush;
            }
        }
        else{
            std::cerr << " Something went wrong " << std::endl;
            return -1;
        }

        ///////////////////////////////////////////////////////////////////////
        // Metrics to do by default
        ///////////////////////////////////////////////////////////////////////

        // Density
        density        = (double) E / ( 0.5 * V * (V -1) );
        density_orig_V = (double) E / orig_max_E;


        ///////////////////////////////////////////////////////////////////////
        // Metrics to only do if requested
        ///////////////////////////////////////////////////////////////////////

        for (auto& m : methods) {

                ///////////////////////////////////////////////////////////////
                // None
                if(m==-1){
                    break;
                }

                else if( m==5 || m==8){
                    ///////////////////////////////////////////////////////////////////////
                    //  Largest connected component sizes
                    // (basically percolation)
                    igraph_t G_cc;
                    largest_connected_component(G, G_cc, cc_count,
                        largest_cc_size, largest2_cc_size);

                    // Spectral Methods
                    if(m == 5){
                        if(largest_cc_size >= minimumpartitionsize){
                            spectral_methods(G_cc,
                               windowsize,
                               minimumpartitionsize,
                               second_eigenvalue,
                               nearly_disconnected_components);
                        }
                    }
                    igraph_destroy(&G_cc);
                }

                ///////////////////////////////////////////////////////////////
                // Maximal Clique Number
                else if(m==4){
                    maximal_cliques(G, minimum_cliquesize,
                                    clique_count, clique_number);
                }

                ///////////////////////////////////////////////////////////////
                // Scale free
                else if(m==3){
                    igraph_plfit_result_t scale_free_result;
                    scale_free_test(G, V, scale_free_result);
                    scale_free_pvalue = scale_free_result.p;
                    scale_free_KS = scale_free_result.D;
                    scale_free_xmin = scale_free_result.xmin;
                    scale_free_alpha = scale_free_result.alpha;
                }

                ///////////////////////////////////////////////////////////////
                // Random Matrix Theory
                else if(m==6){
                    random_matrix_theory(G,
                                         V,
                                         poi_chi_sq_stat,
                                         goe_chi_sq_stat,
                                         poi_chi_sq_pvalue,
                                         goe_chi_sq_pvalue);
                }

                ///////////////////////////////////////////////////////////////
                // Clustering coefficient
                else if(m==7){
                    // plain clustering coefficient
                    igraph_transitivity_undirected(&G, &clustering_coefficient, IGRAPH_TRANSITIVITY_NAN);

                    // threshold G_random an d get random clustering coefficient
                    int threshold_status2 = threshold_graph(t, G_random);
                    //std::cout << "\n" <<  threshold_status2 << "   " << igraph_vcount(&G_random) << "  " << igraph_ecount(&G_random) << std::endl;
                    igraph_transitivity_undirected(&G_random, &clustering_coefficient_r, IGRAPH_TRANSITIVITY_NAN);
                }

                ///////////////////////////////////////////////////////////////
                else{
                    //std::cerr << "Unknown method " << m << std::endl;
                    continue;
                }
        }

        ///////////////////////////////////////////////////////////////////////
        // Make results into a string
        std::stringstream message;
        message << t;
        message << "\t" << V                 << "\t" << E;
        message << "\t" << cc_count;
        message << "\t" << density          << "\t" << density_orig_V;
        message << "\t" << largest_cc_size   << "\t" << largest2_cc_size;
        message << "\t" << clustering_coefficient;
        message << "\t" << clustering_coefficient_r;
        message << "\t" << second_eigenvalue;
        message << "\t" << nearly_disconnected_components;
        message << "\t" << clique_count      << "\t" << clique_number;
        message << "\t" << poi_chi_sq_stat   << "\t" << poi_chi_sq_pvalue;
        message << "\t" << goe_chi_sq_stat   << "\t" << goe_chi_sq_pvalue;
        message << "\t" << scale_free_KS     << "\t" << scale_free_pvalue;
        message << "\t" << scale_free_alpha;
        out << message.str();
        out << std::endl;
    }

    if (methods.find(7) !=methods.end()){
        igraph_destroy(&G_random);
    }

    out.close();
    std::cout << "\nDone. \n" << std::endl;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
//     Command-line arguments                                                //
///////////////////////////////////////////////////////////////////////////////

void help(std::string prog_name){
    std::cerr <<  "\n";
    std::cerr <<  "    Usage: \n";
    std::cerr <<  "    " << prog_name     << " [-OPTIONS]... <GRAPH FILE PATH> <OUTPUT FILE PATH> \n\n";
    std::cerr <<  "    Graph has to be in .ncol format. \n";
    std::cerr <<  "    Output file path is the prefix to the results files, which will be of the form: \n";
    std::cerr <<  "        <OUTPUT FILE PATH>.pid.<method_name>.txt\n\n";
    std::cerr <<  "    Options: \n";
    std::cerr <<  "      -l  --lower                  <value>     lower bound on thresholds to test (default 0.5)\n";
    std::cerr <<  "      -u  --upper                  <value>     upper bound on thresholds to test (default 0.99)\n";
    std::cerr <<  "      -i  --increment              <value>     threshold increment (default 0.01)\n";
    std::cerr <<  "      -w  --windowsize             <value>     sliding window size for spectral method (default 5)\n";
    std::cerr <<  "      -p  --minimumpartitionsize   <value>     minimum size of graph or subgraph after threshold (default 10)\n";
    std::cerr <<  "      -n  --num_samples            <value>     number of samples for significance and power calculations (default NULL)\n";
    std::cerr <<  "      -b  --bonferroni_correction              switch to perform bonferroni corrections in significance and power calculations (default FALSE)\n";
    std::cerr <<  "      -c  --minimum_cliquesize     <value>     minimum size of maximal cliques in maximal clique count (default 5)\n";
    std::cerr <<  "      -m  --methods                <value>     comma separated list of methods (defaults to none)\n";
    std::cerr <<  "                                                   0 - all\n";
    std::cerr <<  "                                                   1 - significance and power calculations (only valid for Pearson CC)\n";
    std::cerr <<  "                                                   2 - local-global\n";
    std::cerr <<  "                                                   3 - scale free\n";
    std::cerr <<  "                                                   4 - maximal cliques\n";
    std::cerr <<  "                                                   5 - spectral methods\n";
    std::cerr <<  "                                                   6 - random matrix theory\n";
    std::cerr <<  "                                                   7 - clustering coefficient\n";
    std::cerr <<  "                                                   8 - percolation\n";
    std::cerr <<  "      -h  --help                               print this help and exit\n";
    std::cerr <<  "\n";
    exit(0);
}

int argument_parser(int argc, char **argv,
    // Mandatory argument definitions
    std::string &infile,
    std::string &outfile,
    // Here flags (options without arguments) and arguments with defined type
    double &lower,
    double &upper,
    double &increment,
    int &windowsize,
    int &minimumpartitionsize,
    int &num_samples,
    bool &bonferroni_corrected,
    int &minimum_cliquesize,
    std::string &methods
    ){

    int next_option;

    const char* const short_options = "hl:u:i:w:p:n:bc:m:" ;
    const struct option long_options[] =
        {    //name,                    has_arg,    flag,        val
            { "help",                   0,          NULL,        'h'},
            { "lower",                  1,          NULL,        'l'},
            { "upper",                  1,          NULL,        'u'},
            { "increment",              1,          NULL,        'i'},
            { "windowsize",             1,          NULL,        'w'},
            { "minimumpartitionsize",   1,          NULL,        'p'},
            { "num_samples",            1,          NULL,        'n'},
            { "bonferroni_correction",  0,          NULL,        'b'},
            { "minimum_cliquesize",     1,          NULL,        'c'},
            { "methods",                1,          NULL,        'm'},
            { NULL, 0, NULL, 0 }
        };

    // Parse options
    while (1) {
        // Obtain an option
        next_option = getopt_long(argc, argv,
                                  short_options, long_options, NULL);

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

            case 'w' : // -w or --windowsize
                windowsize=atoi(optarg);
                break;

            case 'p' : // -p or --minimumpartitionsize
                minimumpartitionsize=atoi(optarg);
                break;

            case 'n' : // -p or --minimumpartitionsize
                num_samples=atoi(optarg);
                break;

            case 'b' : // -p or --minimumpartitionsize
                bonferroni_corrected=1;
                break;

            case 'c' : // -p or --minimumpartitionsize
                minimum_cliquesize=atoi(optarg);
                break;

            case 'm' : // -m or --methods
                methods=optarg;
                break;

            case '?' : // Invalid option
                help(argv[0]); // Return help

            case -1 : // No more options
                break;

            default : // Something unexpected? Aborting
                return(1);
        }
    }

    // Mandatory arguments
    // Current index (optind) < than the total number of arguments
    if(optind == argc){
        std::cerr << "\n Mandatory argument(s) missing\n";
        help(argv[0]);
    }
    // Iterate over rest of the arguments (i.e. in argv[optind])
    if (argc - optind != 2){
         std::cerr << "\n Mandatory argument(s) missing\n";
        help(argv[0]);
    }
    infile = argv[optind];
    optind++;
    outfile = argv[optind];

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
//     Main                                                                  //
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv){

    // register signal SIGABRT and signal handler
    std::signal(SIGABRT, signal_handler);

    // Parse arguments
    // Mandatory argument definitions
    std::string infile;  //input file name
    std::string outfile_prefix;

    // Flags (options without arguments) and arguments with defined type
    double l=0.5;
    double u=0.99;
    double increment=0.01;
    int windowsize=5;
    int minimumpartitionsize=10;
    int minimum_cliquesize=5;
    double min_alpha=0;
    double max_alpha=4;
    double alpha_increment=0.1;

    int num_samples=0;
    double significance_alpha=0.01;
    bool bonferroni_corrected=0;

    std::string str_methods="";

    argument_parser(argc, argv, infile, outfile_prefix,
                    l, u, increment, windowsize,
                    minimumpartitionsize, num_samples,
                    bonferroni_corrected, minimum_cliquesize, str_methods);


    // if argument threshold is given,
    // then we are going to threshold the graph according to the
    // argument method and parameter



    // check arguments
    if(outfile_prefix.empty()) {
        std::cerr << "No output prefix specified. " << std::endl;
        return 0;
    }

    // compare window size to minimumpartitionsize
    if(minimumpartitionsize <= windowsize){
        std::cerr << "Warning: cannot have ";
        std::cerr << "minimumpartitionsize <= windowsize. \n";
        std::cerr << "Using windowsize = 5 and minimumpartitionsize = 10.";
        std::cerr << std::endl;
        minimumpartitionsize = 10;
        windowsize = 5;
    }

    // check that threshold range is good
    if(l>=u){
        std::cerr << "Error in threshold limits: ";
        std::cerr << "cannot have l >= u.\n ";
        std::cerr << "Please restart with corrected lower and upper bounds. ";
        std::cerr << std::endl;
        return 0;
    }

    // decode str_methods
    std::set<int> methods = {};
    if(str_methods == "" ){
        methods.insert({-1});
    }
    else if (str_methods == "0" ){
        methods.insert({1, 2, 3, 4, 5, 6, 7});
    }
    else{
        std::istringstream str_methods_stream(str_methods);
        std::vector<int> methods_listed = {};
        std::string item;
        while(std::getline(str_methods_stream, item, ',')){
            methods_listed.push_back(std::stoi(item));
        }
        methods.insert(methods_listed.begin(), methods_listed.end());
    }

    ///////////////////////////////////////////////////////////////////////
    // Get pid to update output file names
    std::string str_pid = get_str_pid();
    outfile_prefix = outfile_prefix + "." + str_pid + ".";

    // check that output file path exists
    std::ofstream out;
    out.open(outfile_prefix.c_str(), std::ofstream::out);
    // is it open
    if (out.fail()) {
        std::cerr << "Error opening file for writing: " << outfile_prefix << "\n";
        return 0;
    }

    std::cout << "\n";
    std::cout << "------------------------------------------------\n";
    std::cout << "input graph file:      "  << infile << "\n";
    std::cout << "output file prefix:    "  << outfile_prefix << "\n";
    std::cout << "lower threshold:       "  << l << "\n";
    std::cout << "upper threshold:       "  << u << "\n";
    std::cout << "threshold increment:   "  << increment << "\n";
    std::cout << "------------------------------------------------\n";


    ///////////////////////////////////////////////////////////////////////////
    // Type I error (false positive rate) and
    // Type II error (false negative rate) control
    // only for PearsonCC
    // Have to have n - number of samples (not number of variables)
    if(methods.find(1)!=methods.end()){
        std::string outfile_name;
        outfile_name = outfile_prefix + "statistical_errors.txt";
        control_statistical_errors(significance_alpha,
                                  num_samples,
                                  0, //E
                                  bonferroni_corrected,
                                  outfile_name);
        methods.erase(1);
    }

    if (methods.size() == 0){
        return 0;
    }

    // turn on attribute handling
    // for igraph to handle edge weights
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // Load graph
    igraph_t G;
    std::cout << "Loading graph ... " << std::flush;
    read_graph(infile, G, IGRAPH_ADD_WEIGHTS_YES);
    std::cout << "done." << std::endl;

    int status;
    status = thresholdAnalysis(outfile_prefix,
                               G,
                               l,
                               u,
                               increment,
                               windowsize,
                               minimumpartitionsize,
                               minimum_cliquesize,
                               min_alpha,
                               max_alpha,
                               alpha_increment,
                               num_samples,
                               significance_alpha,
                               bonferroni_corrected,
                               methods);

    igraph_destroy(&G);
    return status;
}

///////////////////////////////////////////////////////////////////////////////

