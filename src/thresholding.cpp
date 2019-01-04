// Carissa BLeker
// cbleker@vols.utk.edu

#include <igraph/igraph.h>

#include <vector>     // std::vector
#include <iostream>   // std::cout, std::cerr, std::endl
#include <fstream>    // fopen, fclose (to read igraph)
#include <algorithm>  // std::nth_element, std::min_element, std::max_element
#include <math.h>     // pow, sqrt, fabs
#include <getopt.h>   // commandline argument parsing
#include <stdlib.h>   // atoi, atof
#include <sstream>    // stringstream
#include <cmath>      // std::copysign

///////////////////////////////////////////////////////////////////////////////
//     IO functions                                                          //
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

// write message to outfile if it exists, 
// otherwise write results to std::cout
int output_results(std::string& outfile_name, std::string& message){

    if(!outfile_name.empty()){
        std::ofstream outfile;
        outfile.open(outfile_name.c_str());
        outfile << message; 
        outfile.close();
    }
    else{
        std::cout << message;
    }

    return 0;
}


///////////////////////////////////////////////////////////////////////////////
//     Graph functions                                                       //
///////////////////////////////////////////////////////////////////////////////

// Threshold graph
// by removing edges with abs weight less than "t" 
// and subsequently vertices with no neighbours 
int threshold_graph(double t, igraph_t &G){

    // identify edges to remove
    igraph_vector_t weights;
    igraph_vector_init(&weights, 0);

    igraph_vector_t edge_indices;
    igraph_vector_init(&edge_indices, 0);

    // should this be inside the loop, ie per edge???
    igraph_cattribute_EANV(&G, "weight", igraph_ess_all(IGRAPH_EDGEORDER_ID), &weights);

    igraph_real_t w;
    for (long int i=0; i<igraph_ecount(&G); i++){
        w = VECTOR(weights)[i];
        if(w < t &&  w > -t){
            // add this edge index to list of edges to be deleted
            igraph_vector_push_back(&edge_indices, i);
        }     
    }
    
    // no edges to remove, don't change G
    if(igraph_vector_size(&edge_indices) == 0){
        return 0;
    }
    
    // removing all edges, so just delete all vertices
    else if(igraph_ecount(&G) <= igraph_vector_size(&edge_indices)){
        igraph_delete_vertices(&G, igraph_vss_all());
        return 0;
    }
    
    // remove edges 
    igraph_delete_edges(&G, igraph_ess_vector(&edge_indices));

    std::cout << " Removed " << igraph_vector_size(&edge_indices) << " edges and ";

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&edge_indices);

    // identify and remove degree 0 vertices
    igraph_vector_t vertex_degrees;
    igraph_vector_init(&vertex_degrees, 0);

    igraph_vector_t vertex_indices;
    igraph_vector_init(&vertex_indices, 0);

    igraph_degree(&G, &vertex_degrees, igraph_vss_all(), IGRAPH_ALL, false);

    for(long int i=0; i<igraph_vcount(&G); i++){
        if(VECTOR(vertex_degrees)[i] == 0){
            igraph_vector_push_back(&vertex_indices, i);
        }
    }

    if(igraph_vector_size(&vertex_indices) > 0){
        // remove them
        igraph_delete_vertices(&G, igraph_vss_vector(&vertex_indices));
    }
    
    std::cout << igraph_vector_size(&vertex_indices) << " vertices. " << std::flush;

    igraph_vector_destroy(&vertex_degrees);
    igraph_vector_destroy(&vertex_indices);

    return 0;
}

// Identify largest connected component of the graph and induce
int largest_connected_component(igraph_t &G, igraph_t &G_cc){
    // See also igraph_decompose, but since we only need
    // the largest CC, there is no point inducing all CCs.

    igraph_vector_t membership;
    igraph_vector_init(&membership, 0);
    igraph_vector_t csize;
    igraph_vector_init(&csize, 0);
    igraph_integer_t no; 
    
    igraph_clusters(&G, &membership, &csize, &no, IGRAPH_STRONG);
        
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

    if(max_cc_size < 2){
        return 0;
    }

    // check if there is more than on CC with max_cc_size
    // TODO what to do if there is more that 1?
    int num_max_cc = 0;
    for(int i=0; i<no; i++){
        this_cc_size = VECTOR(csize)[i];
        if(this_cc_size == max_cc_size){
            num_max_cc++;
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
    
    // induce the subgraph of the largest CC
    igraph_induced_subgraph(&G, &G_cc, igraph_vss_vector(&vertices_in_cc), IGRAPH_SUBGRAPH_AUTO);

    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&csize);
    igraph_vector_destroy(&vertices_in_cc);

    return 0;
}

// Fiedler vector: eigen-vector corresponding to first non-zero eigen-value
// Assume connected graph -> 2nd eigenvector
int Fiedler_vector(igraph_t &G, igraph_vector_t &eigenvector, igraph_real_t &eigenvalue){

    std::cout << " Attempting to get the Fiedler vector and value. " << std::flush;
    // make sure G has edges and that G is connected
    if(igraph_ecount(&G) < 1){
        std::cout << " Fielder failed: no edges " << std::endl;
        return 0;
    }

    igraph_bool_t is_connected;
    igraph_is_connected(&G, &is_connected, IGRAPH_STRONG);
    if(!is_connected){
        std::cout <<" Fielder failed: not connected " << std::endl;
        return 0;
    }
    
    // dimension of laplacian = num vertices
    igraph_integer_t V = igraph_vcount(&G);

    igraph_matrix_t laplacian;
    igraph_matrix_init(&laplacian, V, V);

    // init eigen values and vectors
    igraph_vector_t values;
    igraph_matrix_t vectors;
    
    igraph_vector_init(&values, V); // first two eigenvalues will go in here. 
                                    // Size is V since https://github.com/igraph/igraph/issues/1109
    igraph_matrix_init(&vectors, V, 2); // number vertices by 2 eigenvectors

    // Laplacian of G
    igraph_laplacian(&G, &laplacian, NULL, false, NULL);

    // Eigen decomposition for symmetric matrices using LAPACK
    igraph_lapack_dsyevr(&laplacian, IGRAPH_LAPACK_DSYEV_SELECT, 0, 0, 0, 1, 2, 0.0, &values, &vectors, 0); 

    // should be 0.0 (or there abouts)
    std::cout << " 1st eigenvalue: " << VECTOR(values)[0] << std::flush;

    // set eigenvalue and eigenvector of interest
    eigenvalue = VECTOR(values)[1];
    igraph_matrix_get_col(&vectors, &eigenvector, 1);

    // remove laplacian, vectors and values
    igraph_matrix_destroy(&laplacian);
    igraph_vector_destroy(&values);
    igraph_matrix_destroy(&vectors);

    return 0;
}

// See igraph_get_adjacency for unweighted graph
int get_weighted_adjacency(igraph_t &G, igraph_matrix_t &Adj){
    igraph_eit_t edgeit;                             // interate over edges
    long int V = igraph_vcount(&G);               // number vertices / matrix size    

    long int from, to;
    igraph_integer_t ffrom, fto;

    igraph_matrix_resize(&Adj, V, V);
    igraph_matrix_null(&Adj); // set all entries to zero

    // create the edge iterator
    igraph_eit_create(&G, igraph_ess_all(IGRAPH_EDGEORDER_ID), &edgeit);

    // edge weights
    igraph_vector_t weights;
    igraph_vector_init(&weights, 0);
    igraph_cattribute_EANV(&G, "weight", igraph_ess_all(IGRAPH_EDGEORDER_ID), &weights);
    igraph_real_t w;


    while (!IGRAPH_EIT_END(edgeit)){
        long int edge=IGRAPH_EIT_GET(edgeit);
        igraph_edge(&G, (igraph_integer_t) edge, &ffrom, &fto);
        from=ffrom;
        to=fto;
        w = VECTOR(weights)[edge];

        MATRIX(Adj, from, to) = w; 
        
        if (from != to){
            MATRIX(Adj, to, from) = w;
        }

        IGRAPH_EIT_NEXT(edgeit);
     }

    // clean
    igraph_eit_destroy(&edgeit);


    return 0;
}




///////////////////////////////////////////////////////////////////////////////
//     Math/Stat functions                                                   //
///////////////////////////////////////////////////////////////////////////////

// Returns the vector of differences between first and 
// last elements of the windows of size n in x
// from igraph_vector_t to  std::vector
int rolling_difference(igraph_vector_t &x, std::vector<double> &out, int n){

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
double median(std::vector<double> v){
    
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
            double o0;
            double o1;
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
double mean(std::vector<double> v){
    double mean = 0.0;
    int n = v.size();
    for(int i=0; i<n; i++){
        mean = mean + v[i];
    }
    return mean/n;
}

// Standard deviation of a vector
double stddev(std::vector<double> v, double dof=1){
    double n = v.size();
    double v_bar = mean(v);
    double var = 0;

    for(int i=0; i<n; i++){
        var = var + pow( (v[i] - v_bar), 2.0);
    }

    var = var / (n-dof);
    return sqrt(var);
}

// get the exponent to pow value to make a double an int
// Stephen Grady
int get_precision(double k){
    
    int exponent=0;
    bool doneWith0=false;

    std::ostringstream convert;
    convert<<k;
    std::string str=convert.str();

    for (int i = str.size()-1; i > -1; --i){
        if(str[i]=='.'){
            return exponent;
        }
        // first nonzero value, start counting places
        if(str[i]!='0'){
            doneWith0=true;
        }
        if(doneWith0){
            exponent++;
        }
    }
    return 0;
}


// Range from l to u, incrementing by increment
std::vector<double> range(double l, double u, double increment){
    // doubleing point arithmetic
    double precision = pow(10, get_precision(increment));

    int int_increment = static_cast<int>(increment*precision);
    int int_l = static_cast<int>(l*precision);
    int int_u = static_cast<int>(u*precision);

    std::vector<double> out;

    for(int int_t=int_l; int_t<=int_u; int_t+=int_increment){
        out.push_back(int_t/precision);
    }
    
    return out;
}


// Empirical Cumulative Distribution Function (ecdf)
// given x and t, evaluate eCDF of x at t
std::vector<double> ecdf(std::vector<double> x, std::vector<double> t){

    if( !std::is_sorted(x.begin(), x.end()) ){
        std::sort(x.begin(), x.end());
    }
    int len = t.size();
    int n = x.size();   //Number observations

    // x is sorted, so for each x_i in x, F(x_i) = index of x + 1
    // Need to evaluate F(t_i)

    double min_x = x[0];
    double max_x = x.back();

    std::vector<double> F(len);
    double t_i;
    double sum_i;
    for(int i=0; i<len; ++i){
        t_i = t[i];
        if(t_i<min_x){
            F[i] = 0;
        } else if(t_i>max_x){
            F[i] = 1;
        } else {
            sum_i = std::upper_bound(x.begin(), x.end(), t_i) - x.begin(); 
            F[i] = (sum_i)/n;
        }
    }
    return F;
}


double sign(double x){
    if(x == 0){
        return 1.0;
    }
    else{
        return std::copysign(1.0, x);
    }
}

// Accurate Monotonicity Preserving Cubic Interpolation
// In this case x, y is a cdf, and always monotone increasing
// https://epubs-siam-org.proxy.lib.utk.edu/doi/pdf/10.1137/0904045
std::vector<double> spline(std::vector<double> x, std::vector<double> y, std::vector<double> t){
    // x - monotone increasing
    // evaluate at t
    // note: using i=0, ..., n-1 instead of 1, ..., n as in paper

    int n = x.size();
    std::vector<double> dx(n-1);
    std::vector<double> S(n);
    double dy;
    for(int i=0; i < n-1; i++){
        dx[i] = x[i+1] - x[i];
        dy = y[i+1] - y[i];
        S[i] = dy / dx[i];
        std::cout << i << "\t" << x[i] << "\t" << x[i+1] <<  "\tS[i] "<< S[i] <<  std::endl;
    }

    S[n-1] = S[n-2];

    std::cout << "2" << std::endl << std::flush;

    // dy/dx or f_dot approximation
    std::vector<double> f_dot(n);
    double f_dot_i;
    double S_i_min, S_i_max;
    for(int i=0; i < n; i++){
        //j = i+1/2

        // boundaries
        if(i==0){
            f_dot_i = ( (2*dx[i] + dx[i+1])*S[i] - dx[i]*S[i+1] )/( dx[i] + dx[i+1] );
        }
        else if(i==n-1){
            f_dot_i = ( (2*dx[i-1] + dx[i-2])*S[i-1] - dx[i-1]*S[i-2] )/( dx[i-1] + dx[i-2] );
        }
        else{
            // Use Fritsch-Butland
            S_i_min = std::min(S[i], S[i-1]);
            S_i_max = std::max(S[i], S[i-1]); 
        
            std::cout << S_i_min << "\t" << S_i_max << std::endl;
            f_dot_i = 3*S_i_min * S_i_max / (S_i_max + 2*S_i_max); // ( dx[i-2]*S[i-1] + dx[i-1]*S[i-2] ) / ( x[i] - x[i-2]  );
        }
        
        // Hyman filter: eqn 2.6
        double sigma = sign(f_dot_i);
        std::cout << "i=" << i <<  "\tf_dot_i= " << f_dot_i << "\tsigma= " << sigma << std::endl;
        if(sigma > 0.0){
            f_dot_i = std::min( std::max(0.0, f_dot_i),  3*std::min(fabs(S[i-1]), fabs(S[i])) );
        } 
        else{
            f_dot_i = std::max( std::min(0.0, f_dot_i), -3*std::min(fabs(S[i-1]), fabs(S[i])) );
        }

        f_dot[i] = f_dot_i;
    }
    
    std::cout << "3" << std::endl << std::flush;

    // evalute
    int t_n = t.size();
    std::vector<double> P(t_n);
    int i, j;
    double t_j; 
    double c1, c2, c3, c4, p;

    for(j=0; j<t_n; j++){

        // x[i] <= t[j] < x[i+1]
        t_j = t[j];
        i = std::lower_bound(x.begin(), x.end(), t_j) - x.begin() -1;

        std::cout << "j=" << j << "\ti=" << i << "\tt_i=" << t_j << "\tx_i=" << x[i] << std::endl << std::flush;

        c1 = y[i];
        c2 = f_dot[i];
        c3 = (3*S[i] - f_dot[i+1] - 2*f_dot[i]) / (dx[i]);
        c4 = (2*S[i] - f_dot[i+1] - f_dot[i]) / pow(dx[i], 2.0); 
        p = c1 + (t_j - x[i])*c2 + pow(t_j - x[i], 2.0)*c3 + pow(t_j - x[i], 3.0)*c4;

        P[j] = p;
    }

    return P;
}



///////////////////////////////////////////////////////////////////////////////
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
		rolling_difference(eigenvector, window_differences, windowsize);

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

    // results go here
    std::vector<double>  stat_per_t(num_increments); 

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

        // Here
        igraph_matrix_t A;
        igraph_matrix_init(&A, V, V);
        std::cout << " made adj matrix "  << std::flush;
        get_weighted_adjacency(G, A);
        
        igraph_vector_t eigenvalues;
        igraph_vector_init(&eigenvalues, V); // eigenvalues will go in here. 

        std::cout << " get eigenvalues" << std::flush; 
        igraph_lapack_dsyevr(&A, IGRAPH_LAPACK_DSYEV_ALL, 0, 0, 0, 0, 0, 0.0, &eigenvalues, NULL, 0);

        // sort
        igraph_vector_sort(&eigenvalues);

        // make into std::vector and remove duplicates
        // also keep largest absolute value eigenvalue
        double largest_eigenvalue = fabs(VECTOR(eigenvalues)[0]);
        
        std::vector<double> eigenvalues_sorted;
        eigenvalues_sorted.push_back(VECTOR(eigenvalues)[0]);

        double this_eigenvalue;
        for(int i=1; i<igraph_vector_size(&eigenvalues); i++){
            this_eigenvalue = VECTOR(eigenvalues)[i];
            if(VECTOR(eigenvalues)[i-1] != this_eigenvalue){
                eigenvalues_sorted.push_back(this_eigenvalue);
                if(fabs(this_eigenvalue) > largest_eigenvalue){
                    largest_eigenvalue = fabs(this_eigenvalue);
                }
            }
        }

        //std::cout << "\n";
        //int zzz = eigenvalues_sorted.size();
        //for (int i =0; i < zzz; i++){
        //    std::cout << eigenvalues_sorted[i] << "\n";
        //}

        // scale eigenvalues: divide by largest absolute value eigenvalue
        for(int i=0; i<eigenvalues_sorted.size(); i++){
            eigenvalues_sorted[i] = eigenvalues_sorted[i]/largest_eigenvalue;
        } 

        int number_fit_points = floor(0.75*eigenvalues_sorted.size()); //TODO
        std::cout << "\t" << number_fit_points << "\t"; 
        
        // equally spaced values between range of eigenvalues for CDF
        double increment = (eigenvalues_sorted.back() - eigenvalues_sorted[0])/number_fit_points;
        std::vector<double> t = range(eigenvalues_sorted[0], eigenvalues_sorted.back(), increment);

        std::vector<double> F = ecdf(eigenvalues_sorted, t);
        
        //std::cout << "\n";
        //int zz = eigenvalues_sorted.size();
        //for (int i =0; i < zz; i++){
        //    std::cout << eigenvalues_sorted[i] << "\n";
        //}
        //std::cout << "\n";
        //int z = F.size();
        //for (int i =0; i < z; i++){
        //    std::cout << t[i] << " " << F[i] << "\n";
        //}

        // Find smooth distribution of eigenvalues by fitting a spline to the CDF and evaluating at t values
        spline = 
        new_eigenvalues = 

        




        //std::vector<double> NNSD;
        //rolling_difference(new_eigenvalues, NNSD, 1);




        // Chi2 test
        



        stat_per_t[i_t] = V;
        was_tested_per_t[i_t] = true;
    }

    igraph_destroy(&G);

    std::cout << "\nDone\n" << std::endl;

    // make results into a string
    std::stringstream message;
    message << "threshold\tstat\n";
    for(int i=0; i < num_increments; i++){
        if(was_tested_per_t[i]){
            message << t_vector[i] << "\t" << stat_per_t[i] << "\n";
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
