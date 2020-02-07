# Graph Thresholding

Repository for thresholding code, based on dissertation research. We use
a number of methods to analyse a weighted graph in order to suggest an
optimal threshold for the graph.

The components of this repository are:
 - C++ threshold analysis code
 - C++ threshold code
 - Bash threshold script
 - script to output a histogram of the edge weights of a graph
 - Python3 Jupyter notebook to analyse the results of the thresholding analysis code

---
## Installation
### C++ program

Compiling also installs [igraph](igraph.org/c/)-0.8.0, so you will need to have
igraph dependency (libxml2) installed beforehand
(see [below](#optionalinstructionsforlibxml2)).

    git clone git@github.com:carissableker/thresholding.git
    cd thresholding

To compile:

    make

Executables will be in `bin` folder.

### Scripts

The bash script `absolute_global_threshold` uses awk, and should be executable.

The bash script `edge_weight_histogram` also uses awk to calculate the
bin counts, and uses Python3 (including seaborn, matplotlib, and pandas) to produce the visualization.

### Jupyter notebook

TODO

---
## Usage

Graphs need to be in `.ncol` format to be read correctly
[igraph_read_graph_ncol](https://igraph.org/c/doc/igraph-Foreign.html#igraph_read_graph_ncol).
The format is defined by
[Large Graph Layout](http://lgl.sourceforge.net/#FileFormat).
In this application, it is a simple weighted, whitespace separated edge list. Vertex names cannot contain whitespace.

  vertex<sub>1</sub>⇥vertex<sub>2</sub>⇥weight<sub>1,2</sub> <br>
  vertex<sub>1</sub>⇥vertex<sub>3</sub>⇥weight<sub>1,3</sub> <br>
  &middot; <br>
  &middot; <br>
  &middot; <br>

### 1. Graph analysis

It is recommended to first view a histogram of the edge weights:

```bash
$ ./bin/edge_weight_histogram <GRAPH FILE PATH> <OUTPUT FILE PATH> <BIN WIDTH>
```
This will create a tab separated file with bin counts in
`<OUTPUT FILE PATH>.tsv.`
and a SVG figure of the histogram at `<OUTPUT FILE PATH>.svg`.
The histogram can be used to decide on lower and upper bounds on the thresholds.

For analysis of the graph, run the `analysis` program.

```bash
$ ./bin/analysis --help

    Usage:
    ./bin/analysis [-OPTIONS]... <GRAPH FILE PATH> <OUTPUT FILE PATH>

    Graph has to be in .ncol format.
    Output file path is the prefix to the results files, which will be of the form:
        <OUTPUT FILE PATH>.pid.<method_name>.txt

    Options:
      -l  --lower                  <value>     lower bound on thresholds to test (default 0.5)
      -u  --upper                  <value>     upper bound on thresholds to test (default 0.99)
      -i  --increment              <value>     threshold increment (default 0.01)
      -w  --windowsize             <value>     sliding window size for spectral method (default 5)
      -p  --minimumpartitionsize   <value>     minimum size of graph or subgraph after threshold (default 10)
      -n  --num_samples            <value>     number of samples for significance and power calculations (default NULL)
      -b  --bonferroni_correction              switch to perform bonferroni corrections in significance and power calculations (default FALSE)
      -c  --minimum_cliquesize     <value>     minimum size of maximal cliques in maximal clique count (default 5)
      -m  --methods                <value>     comma separated list of methods (defaults to none)
                                                   0 - all
                                                   1 - maximal cliques
                                                   2 - scale free
                                                   3 - spectral methods
                                                   4 - random matrix theory
                                                   5 - clustering coefficient
                                                   6 - local-global
                                                   7 - significance and power calculations (only valid for Pearson CC)
      -h  --help                               print this help and exit
```

### 2. Analysis of results

The output of graph analysis is a number of files containing statistics and
metrics of the graph. A Python3 Jupyter notebook is supplied to
interactively analyse these files. This notebook can be found at ./TODO


### 3. Threshold

For simple thresholding, the fastest is probably to use
the `./bin/absolute_global_threshold` script.

```bash
$ ./bin/absolute_global_threshold <input.ncol> <t> <output.ncol>
```

This thresholds the graph at |t|, i.e. removes all edges with absolute value smaller or equal to t.

Alternatively, `threshold` can also be used to threshold the graph

```bash
$ ./bin/threshold  --help

    Usage:
    ./bin/threshold [-OPTIONS]... <GRAPH FILE PATH> <OUTPUT FILE PATH>

    Graph has to be in .ncol format.

    One of the following options have to be given:
    Options:
      -a  --absolute              <value>     Threshold graph at absolute of <value>
      -l  --local-global          <value>     Use local-global method to threshold with alpha = <value>
      -r  --rank                  <value>     Use top <value> ranked edges per vertex to threshold graph
      -h  --help                              print this help and exit
```

---
## Example:

```bash
$ cd example
$ ls
ecoliGraph.ncol
```




```bash
$ ../bin/edge_weight_histogram ecoliGraph.ncol ecoliGraph-edge-hist 0.01
Edge distribution at ecoliGraph-edge-hist.tsv and ecoliGraph-edge-hist.svg
$ head ecoliGraph-edge-hist.tsv
bin_start bin_end bin_count
-1.00000  -0.99000  974
-0.99000  -0.98000  2718
-0.98000  -0.97000  4375
-0.97000  -0.96000  5813
```

![example-histogram](./doc/figures/ecoliGraph-edge-hist.svg?raw=true "Example histogram")


This graph has already been thresholded at 0.5 for speed and file size.


```bash

```
```bash

```


---
## TODOs:

### Algorithmic
* Change Chi2 test of GOE to KS test for continuous distributions
* Add Spearman significance
* Consider positive and negative values separately

### Implementation wise
* Test cases
* Improve memory allocation and deallocation
* Make increment loop optional
* Consider support for Windows/Mac
* Generally improve C++ code (with help from a C++ developer)
* Fix igraph dependency in Makefile
* ?Other input/output formats

---
## Other

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