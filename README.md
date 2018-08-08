# thresholding


Repository for thresholding code. 


## To compile:

Depends on [igraph](igraph.org/c/).  

    make

## To use:

    ./threshold --help
     
	Usage: 
       ./threshold [-OPTIONS]... <GRAPH FILE PATH>

     Graph has to be in .ncol format. 

     Options: 
      -o  --out                      <filename>         path to store results
                                                            if not given, results are sent to stdout
      -l  --lower                    <value>            lower bound on thresholds to test (default 0.5)
      -u  --upper                    <value>            upper bound on thresholds to test (default 0.99)
      -i  --increment                <value>            threshold increment (default 0.01)
      -w  --windowsize               <value>            sliding window size for spectral method (default )5
      -p  --minimumpartitionsize     <value>            minimum size of graph or subgraph after thresholding (default 5)
      -m  --method                   [1|2|3|4]          method  (default = 1)
                                                             1 - Spectral method
                                                             2 - Clique ratio (generalised clique doubling)
		                                          			 3 - Density
                                                             4 - Percolation (not implemented)
                                                             5 - Random matrix theory (not implemented)
      -h  --help                                        print this help and exit


