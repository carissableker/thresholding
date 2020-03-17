#include "significance.h"




double control_statistical_errors(double significance_alpha,
								  int num_samples,
							   	  int E,
						   	      bool bonferroni_corrected,
				                  std::string& outfile_name
				                  ){
	if(num_samples <= 4){
    	std::cerr << "**************************************************\n";
		std::cerr << "Warning: Cannot compute significance errors with num_samples <= 4.\n";
    	std::cerr << "**************************************************\n";

		return 0;
	}


	std::ofstream out;
    out.open(outfile_name.c_str(), std::ofstream::out);
    // is it open
    if (out.fail()) {
        std::cerr << "Error opening file for writing: " << outfile_name << "\n";
        return 0;
    }

	///////////////////////////////////////////////////////////////////////
	// Control Type I
    double k = num_samples-2;
    // calculate max p-value
    double max_pvalue= significance_alpha;
    if(bonferroni_corrected){
        max_pvalue = significance_alpha / double(E);
    }
    //std::cout << k << "\t" << 1-max_pvalue/2;

    // backwards calculate t-statistic using inverse student's T CDF
    double t_crit = alglib::invstudenttdistribution(k, 1-max_pvalue/2);

    // calculate correlation  value (i.e. the threshold based on alpha)
    double r_crit = t_crit / sqrt( pow(t_crit, 2) + k);

    std::stringstream header_type1;
    header_type1 << "# Critical Pearson correlation for obtaining alpha significance of ";
    header_type1 << significance_alpha << " with sample size of " << num_samples << " is ";
    header_type1 << r_crit;
    out << header_type1.str();
    out << std::endl;

	///////////////////////////////////////////////////////////////////////
	// Control Type II

    std::stringstream header_type2;
    header_type2 << "# Power analysis \n";
    header_type2 << "r";
    header_type2 << "\t" << "power";
    out << header_type2.str();
    out << std::endl;

	// This code based on pwr.r.test (sort-of)
	// Under H_a: rho=r,
	// z_r is std norm with mean = z_rc and std = 1/sqtr(n-3)

	double z_crit  = fisher_transform(r_crit, num_samples);

	// Range from l to u, incrementing by increment
	std::vector<double> r_range = range(0.01, 0.99, 0.01);

    for(int i=0; i<r_range.size(); i++){
    	double r_alternative=r_range[i];
    	double z_r = fisher_transform(r_alternative, num_samples);

    	double power = 1 - alglib::normaldistribution( ( z_crit - z_r) * sqrt(double(num_samples) - 3) )
    					 + alglib::normaldistribution( (-z_crit - z_r) * sqrt(double(num_samples) - 3) );

        std::stringstream message;
        message << r_alternative;
        message << "\t" << power;
        out << message.str();
        out << std::endl;
    }

    out.close();
    return 0;
}