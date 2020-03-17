#include "random_matrix_theory.h"

int random_matrix_theory(igraph_t& G,
                         igraph_integer_t &V,
                         double  &poi_chi_sq_stat,
                         double  &goe_chi_sq_stat,
                         double  &poi_chi_sq_pvalue,
                         double  &goe_chi_sq_pvalue
                         ){

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

    double this_eigenvalue, previous_eigenvalue;

    this_eigenvalue = VECTOR(eigenvalues)[0] / largest_eigenvalue + 1;
    eigenvalues_sorted.push_back(this_eigenvalue);
    previous_eigenvalue = this_eigenvalue;
    for(int i=1; i<igraph_vector_size(&eigenvalues); i++){
        this_eigenvalue = VECTOR(eigenvalues)[i] / largest_eigenvalue + 1;
        if(fabs(previous_eigenvalue - this_eigenvalue) > 0.00001){
            eigenvalues_sorted.push_back(this_eigenvalue);
            //std::cout << this_eigenvalue << "\n";
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

    //std::cout << "\nNNSD_size " << NNSD_size << "\n";
    for (int i=0; i < NNSD_size; i++){
        NNSD[i] = NNSD[i] * n ;
    //    std::cout << "\n" << i << "\t" << NNSD[i];
    }

    // Chi2 test on NNSD
    // https://github.com/spficklin/RMTGeneNet/blob/master/threshold/methods/RMTThreshold.cpp
    // https://www.statisticshowto.datasciencecentral.com/goodness-of-fit-test/
    // https://static-content.springer.com/esm/art%3A10.1186%2F1471-2105-8-299/MediaObjects/12859_2006_1671_MOESM3_ESM.pdf

    std::sort(NNSD.begin(), NNSD.end());

    // Need to discretise continuous NNSD
    // Use histogram equalization:
    // observed bin frequency, will always be 5 (except for last bin)

    double bin_start = 0;
    double bin_end = 0;

    std::vector<double> bin_start_vector;
    std::vector<double> bin_end_vector;
    std::vector<double> bin_count_vector;

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
            //std::cout << "\n h\t" << h << "\tbin end\t" << bin_end << "\t" << NNSD[h-1] << "\n";
            bin_end_vector[no_bins-1] = bin_end;
            bin_count_vector[no_bins-1] = bin_count_vector[no_bins-1] + this_bin_count;

            NNSD_i = h;
        }
    }

    int dof = no_bins -1;

    if(dof > 1){

        poi_chi_sq_stat = 0;
        goe_chi_sq_stat = 0;

        for(int b=0; b<no_bins; b++){

            bin_start = bin_start_vector[b];
            bin_end = bin_end_vector[b];
            observed_count = bin_count_vector[b];

            //std::cout << b << "\t" << bin_start << "\t " << bin_end << "\t " << observed_count << "\t ";

            // If Poisson, then
            expected_count = NNSD_size * ( poisson(bin_start, bin_end) );
            poi_chi_sq_stat += pow(observed_count - expected_count, 2) / expected_count;
            //std::cout << expected_count << "\t ";

            // If GOE, then
            expected_count = NNSD_size * ( goe(bin_start, bin_end) );
            goe_chi_sq_stat += pow(observed_count - expected_count, 2) / expected_count;
            //std::cout << expected_count << std::endl;
        }

        // p-value = area under right hand tail
        // http://www.alglib.net/specialfunctions/distributions/chisquare.php

        if(std::isnormal(poi_chi_sq_stat)){
            poi_chi_sq_pvalue = alglib::chisquarecdistribution(dof, poi_chi_sq_stat);
        }
        if(std::isnormal(goe_chi_sq_stat)){
            goe_chi_sq_pvalue = alglib::chisquarecdistribution(dof, goe_chi_sq_stat);
        }
    }

    return 0;
}

