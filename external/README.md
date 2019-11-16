# External libries

## igraph 

Approximate method:

    wget https://igraph.org/nightly/get/c/igraph-0.7.1.tar.gz -O ./external/igraph-0.7.1.tar.gz
    tar -zxvf external/igraph-0.7.1.tar.gz
    git add external/igraph-0.7.1/*
    find ./external/igraph-0.7.1/ -type f -exec git update-index --assume-unchanged '{}' \;

