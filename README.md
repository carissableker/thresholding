# thresholding


Repository for thresholding code. 


## To compile:

Depends on [igraph](igraph.org/c/).  

    make

## To use:

    ./threshold --help

    Usage: 
      ./threshold [-OPTIONS]... <GRAPH FILENAME>

    GRAPH is in .ncol format. 

     Options: 
      -o  --out                      <filename>        path to store results
                                                             if not given, results are sent to stdout
      -l  --lower                    <value>            lower bound on thresholds to test (default 0.5)
      -u  --upper                    <value>            upper bound on thresholds to test (default 0.99)
      -i  --increment                <value>            threshold increment (default 0.01)
      -w  --windowsize               <value>            sliding window size for spectral method (default 5)
      -p  --minimumpartitionsize     <value>            minimum size of graph or subgraph after thresholding (default 5)
      -m  --method                   [1|2|3]            method  (default = 1) Spectral methohs only. 
