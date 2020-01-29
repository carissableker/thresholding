#ifndef SIGNIFICANCE_H
#define SIGNIFICANCE_H

#include <fstream>    // ofstream

#include "alglib/specialfunctions.h"  // ALGLIB header that contains invstudenttdistribution

#include "math_ext.h"



double control_statistical_errors(double,
								  int,
							   	  int,
						   	      bool,
				                  std::string&
				                  );


#endif