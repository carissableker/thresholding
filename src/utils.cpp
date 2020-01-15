#include <utils.h> 

/*
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
*/

// incrementally write to an already stream
int output_results(std::ofstream outfile, std::string& message){
    outfile << message; 
    return 0;
}


