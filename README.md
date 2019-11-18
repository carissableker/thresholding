# thresholding

Repository for thresholding code. 

## Installation

Compiling also installs [igraph](igraph.org/c/)-0.7.1, so you will need to have
igraph dependency (libxml2) installed beforehand (see below). 

    git clone git@github.com:carissableker/thresholding.git
    cd thresholding

The dev branch is the only branch that works at the moment:

    git checkout dev
    git pull origin dev

To compile: 

    make

Executable will be in `bin` folder. 

### Optional instructions for libxml2:

If libxml2 is not already installed, the following is one way of installing it. 
Another (preferable and easier) option is `apt-get`. 

    wget ftp://xmlsoft.org/libxml2/libxml2-<version>.tar.gz
    tar -xzvf libxml2-<version>.tar.gz 
    ./configure  --without-python --prefix=/my_optional path/
    make
    make install

And if you used a prefix: 

    echo 'export PATH="/my_optional_path/bin:$PATH"
    export LD_LIBRARY_PATH="/my_optional_path/lib:$LD_LIBRARY_PATH"
    export LD_RUN_PATH="/my_optional_path/lib:$LD_RUN_PATH"' >> ~/.bashrc
    source ~/.bashrc

## Usage

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
      -w  --windowsize               <value>            sliding window size for spectral method (default 5)
      -p  --minimumpartitionsize     <value>            minimum size of graph or subgraph after thresholding (default 5)
      -h  --help                                        print this help and exit


