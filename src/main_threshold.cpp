// made on dev
// Carissa BLeker
// cbleker@vols.utk.edu

#include <igraph/igraph.h>

#include <vector>     // std::vector
#include <iostream>   // std::cout, std::cerr, std::endl
#include <fstream>    // fopen, fclose (to read igraph), ofstream
#include <getopt.h>   // commandline argument parsing
#include <stdlib.h>   // atoi, atof
#include <sstream>    // stringstream
#include <string>     // getline

#include "utils.h"
#include "math_ext.h"
#include "igraph_ext.h"

#include "local_global.h"
#include "local_rank.h"

///////////////////////////////////////////////////////////////////////////////

int threshold(std::string& outfile,
              igraph_t &G,
              std::string method,
              double absolute,
              double local_global_alpha,
              int rank){

    igraph_integer_t E = igraph_ecount(&G); // number edges
    igraph_integer_t V = igraph_vcount(&G); // number vertices

    std::cout << "Number vertices:  " << V << "\n";
    std::cout << "Number edges:     " << E << std::endl;

    if (method == "absolute"){
        threshold_graph(absolute, G);
        write_graph(outfile, G);
    }
    else if (method == "local-global"){
        igraph_t new_G;
        double mean_k;
        local_global_pruning(G, local_global_alpha, new_G, mean_k);
        write_graph(outfile, new_G);
    }
    else if (method=="rank"){
        igraph_t new_G;
        double mean_k;
        local_global_pruning(G, local_global_alpha, new_G, mean_k);
        write_graph(outfile, new_G);
   }
    else {
        std::cerr << "something went wrong...";
    }

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
    std::cerr <<  "      -l  --lower                 <value>     lower bound on thresholds to test (default 0.5)\n";
    std::cerr <<  "      -u  --upper                 <value>     upper bound on thresholds to test (default 0.99)\n";
    std::cerr <<  "      -i  --increment             <value>     threshold increment (default 0.01)\n";
    std::cerr <<  "      -w  --windowsize            <value>     sliding window size for spectral method (default 5)\n";
    std::cerr <<  "      -p  --minimumpartitionsize  <value>     minimum size of graph or subgraph after threshold (default 10)\n";
    std::cerr <<  "      -m  --methods               <value>     comma separated list of methods (defaults to none)\n";
    std::cerr <<  "                                                  0 - all\n";
    std::cerr <<  "                                                  1 - maximal cliques\n";
    std::cerr <<  "                                                  2 - scale free\n";
    std::cerr <<  "                                                  3 - spectral methods\n";
    std::cerr <<  "                                                  4 - random matrix theory\n";
    std::cerr <<  "                                                  5 - clustering coefficient\n";
    std::cerr <<  "                                                  6 - local-global\n";
    std::cerr <<  "      -h  --help                              print this help and exit\n";
    std::cerr <<  "\n";
    exit(0);
}

int argument_parser(int argc, char **argv,
    // Mandatory argument definitions
    std::string &infile,
    std::string &outfile,
    // Here flags (options without arguments) and arguments with defined type
    std::string &method,
    double &absolute,
    double &local_global_alpha,
    int &rank){


    const char* const short_options = "ha:l:u:r:" ;
    const struct option long_options[] =
        {    //name,                    has_arg,    flag,        val
            { "help",                   0,          NULL,        'h'},
            { "absolute",               1,          NULL,        'a'},
            { "local-global",           1,          NULL,        'l'},
            { "rank",                   1,          NULL,        'r'},
            { NULL, 0, NULL, 0 }
        };

    // Parse option
    int option;
    option = getopt_long(argc, argv,
                              short_options, long_options, NULL);

    switch (option){

        case 'h' : // -h or --help
            help(argv[0]);
            break;

        case 'a' : // -a or --absolute
            absolute=atof(optarg);
            method="absolute";
            break;

         case 'l' : // -l or --local-global
            local_global_alpha=atof(optarg);
            method="local-global";
            break;

         case 'r' : // -r or --rank
            rank=atoi(optarg);
            method="rank";
            break;

        case '?' : // Invalid option
            help(argv[0]); // Return help
            break;

        case -1 : // No more options
            break;

        default : // Something unexpected? Aborting
            return(1);
    }

    // Mandatory arguments
    // Current index (optind) < than the total number of arguments
    if(optind == argc){
        std::cerr << "\n Mandatory argument(s) missing\n";
        help(argv[0]);
    }
    // Iterate over rest of the arguments (i.e. in argv[optind])
    infile = argv[optind];
    optind++;
    outfile = argv[optind ];

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
//     Main                                                                  //
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv){

    // Parse arguments
    // Mandatory argument definitions
    std::string infile;  //input file name
    std::string outfile;

    // Flags (options without arguments) and arguments with defined type
    std::string method;
    double absolute;
    double local_global_alpha;
    int    rank;

    argument_parser(argc, argv, infile, outfile,
                    method, absolute, local_global_alpha, rank);

    // check arguments
    if(outfile.empty()) {
        std::cerr << "No output file name specified. " << std::endl;
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////

    std::cout << "\n";
    std::cout << "------------------------------------------------\n";
    std::cout << "input graph file:      "  << infile << "\n";
    std::cout << "output graph file      "  << outfile << "\n";
    std::cout << "threshold method       "  << method << "\n";
    std::cout << "------------------------------------------------\n";

    // turn on attribute handling
    // for igraph to handle edge weights
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // Load graph (names = true)
    igraph_t G;
    std::cout << "Loading graph ... " << std::flush;
    read_graph(infile, G, IGRAPH_ADD_WEIGHTS_YES, true);
    std::cout << "done." << std::endl;

    int status;
    status = threshold(outfile, G, method, absolute, local_global_alpha, rank);
    return status;
}

///////////////////////////////////////////////////////////////////////////////

