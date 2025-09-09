# RAHL

A tool for indexing undirected edge-weighted graphs, such as road networks, to speed up distance queries.
It consists of the following main files:

* road_network.h / road_network.cpp: core library
* util.h / util.cpp: library with additional tools
* main.cpp: run experiments on one or more graphs and summarize results

## Compile & Run

Files can be compiled using `make`.
Experiments can then be run with

```
./cut [beta] graph_1.gr ... graph_n.gr
```

The balance parameter `beta` is optional, with 0.15 as default.
Graph files for testing, as well as a description of the expected file format, can be found at http://www.diag.uniroma1.it/~challenge9/download.shtml.
