# Graph Thresholding

Repository for thresholding code, based on dissertation research. We use
a number of methods to analyse a weighted graph in order to suggest an
optimal threshold for the graph.

The components of this repository are:
 - C++ threshold anlysis code
 - script to output a histogram of the edge weights of a graph
 - Python3 Jupyter notebook to analyse the results of the thresholding analysis code
 - script to output a threshold graph


## Installation
---

### C++ program

Compiling also installs [igraph](igraph.org/c/)-0.7.1, so you will need to have
igraph dependency (libxml2) installed beforehand (see [below](#optionalinstructionsforlibxml2)).

    git clone git@github.com:carissableker/thresholding.git
    cd thresholding

To compile:

    make

Executable will be in `bin` folder.


### Scripts

The bash script `absolute_global_threshold` uses awk, and should be executable.

The bash script `edge_weight_histogram` also uses awk to calculate the bin counts,
and uses Python3 to produce the vizualization.

### Jupter notebook




## Usage
---

### 1. Graph analysis

It is recommended to first view a histogram of the edge weights:

```bash
$ ./bin/edge_weight_histogram <GRAPH FILE PATH> <OUTPUT FILE PATH> <BIN WIDTH> <OUTPUT FILE PATH>
```
This will create a tab separated file with bin counts in `<OUTPUT FILE PATH>.tsv.`
and a SVG figure of the histogram at `<OUTPUT FILE PATH>.svg`.
The histogram can be used to decide on lower and upper bounds on the thresholds.

For analysis of the graph in terms of threshold options, run the `threshold` program.

```bash
$ ./bin/threshold --help

Usage:
./bin/threshold [-OPTIONS]... <GRAPH FILE PATH> <OUTPUT FILE PATH>

Graph has to be in .ncol format.
Output file path is the prefix to the results files, which will be of the form:
    <OUTPUT FILE PATH>.pid.<method_name>.txt

Options:
  -l  --lower                 <value>     lower bound on thresholds to test (default 0.5)
  -u  --upper                 <value>     upper bound on thresholds to test (default 0.99)
  -i  --increment             <value>     threshold increment (default 0.01)
  -w  --windowsize            <value>     sliding window size for spectral method (default 5)
  -p  --minimumpartitionsize  <value>     minimum size of graph or subgraph after thresholding (default 10)
  -m  --methods               <value>     comma seperated list of methods (defaults to none)
                                              0 - do all methods
                                              1 - maximal clicques
                                              2 - scale free
                                              3 - spectral methods
                                              4 - random matrix theory
                                              5 - clustering coefficient
                                              6 - local-global (guzzi2014)
                                              7 - statistical significance
                                              8 - local (rank)
  -h  --help                              print this help and exit
```

### 2. Analysis of results

The output of graph analysis is a number of files containing statistics and
metrics of the graph. A Python3 Jupyter notebook is supplied to
interactively analyse these files. This notebook can be found at ./


### 3. Threshold

For simple thresholding, the fastest is probably to use the `./absolute_global_threshold` script.

```bash
Usage:  absolute_global_threshold <input.ncol> <t> <output.ncol>
```

This thresholds the graph at |t|, i.e. removes all edges with absolute value smaller or equal to t.


Alternatively, `threshold` can also be used to threshold the graph. (See usage above).


## To do:
---

* Change Chi2 test of GOE to KS test for continuous distributions
* Add Spearman significance
* Make increment loop optional
* Consider positive and negative values separately
* Consider support for Windows/Mac


## Other
---

### Optional instructions for libxml2

If libxml2 is not already installed, the following is one way of installing it.
Another (preferable and easier) option is `apt-get`.

```bash
    wget ftp://xmlsoft.org/libxml2/libxml2-<version>.tar.gz
    tar -xzvf libxml2-<version>.tar.gz
    ./configure  --without-python --prefix=/my_optional path/
    make
    make install
```

And if you used a prefix:

```bash
    echo 'export PATH="/my_optional_path/bin:$PATH"
    export LD_LIBRARY_PATH="/my_optional_path/lib:$LD_LIBRARY_PATH"
    export LD_RUN_PATH="/my_optional_path/lib:$LD_RUN_PATH"' >> ~/.bashrc
    source ~/.bashrc
```