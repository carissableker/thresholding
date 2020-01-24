#ifndef UTILS_H
#define UTILS_H

#include <vector>     // std::vector
#include <iostream>   // std::cout, std::cerr, std::endl
#include <fstream>    // fopen, fclose (to read igraph)
#include <stdlib.h>   // atoi, atof
#include <sstream>    // stringstream
#include <unistd.h>    // getpid()

// Results IO
int output_results(std::string&, std::string&);

// get pid and return as a string
std::string get_str_pid();


#endif