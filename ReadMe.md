The following is the code to run Wormhole, a fast algorithm to compute estimated distances between vertex pairs in graphs that exhibit a core-periphery structure. This is an accompaniment to the manuscript _A Sublinear Algorithm for Approximate Shortest Paths in Large Networks_ by Sabyasachi Basu, Nadia Kōshima, Talya Eden, Omri Ben-Eliezer, C. Seshadhri (currently in submission).

# Setup

Navigate to the root directory and run the following commands on terminal:

```
make
mkdir bin
mkdir graphs
mkdir input
mkdir results
```

# Running the algorithm

## Sanitization and CGraph Creation

Graph datasets sourced from the Internet are often haphazardly and inefficiently numbered. The objective of sanitization is to get rid of redundant vertices and repeated edges. In our implmentation, we only remove repeated edges. Our implementation is for undirected graphs. We run `sanitize.cpp` to this effect.

The input to `sanitize` is a file containing the graph as an edgelist. The output of the sanitization procedure is stored as a `.edges` file. The first line has the number of vertices (`nVertices`) and the number of edges(`nEdges`), separated by whitespace. The first line is followed by the edge list in a typical COO format.

```
<nVertices> <nEdges>
<v1> <v2>
<v1> <v3>
...
```

To sanitize an edgelist `INPUTNAME.EXTENSION` stored under the `graphs/` directory, the command is
`./exe/sanitize INPUTNAME .EXTENSION`;
the output is saved as `graphs/INPUTNAME-sanitized.edges`.

To create the CGraph from a sanitized edgelist named `graphs/GRAPHNAME.edges`, run
`./exe/createCGraph GRAPHNAME`
The output stores the graph as a CSR matrix in `bin/GRAPHNAME_CGraph.bin`

## Generate Input

To generate a sample of edge queries, run
`./exe/generateInput GRAPHNAME`
By default, it creates 10000 random queries. Each query is a vertex pair `s t`, and reflexive queries of the form `s s` are strictly avoided. The output is stored as `input/GRAPHNAME-input.txt`.

## Running Wormhole

```
./exe/generateL0Seed GRAPHNAME seed
./exe/sample L0 GRAPHNAME_seed_SIZE 0 10000
```

The output of the first line is a set of files of the form `bin/GRAPHNAME_side_SIZE_L0.bin` and `bin/GRAPHNAME_side_SIZE_L1.bin`; where `SIZE` refers to the percentage of the graph chosen in the inner ring; we call this the size. As set up, this size is in the set $\{0.5, 1.0, 1.5, \ldots, 10.0\}$; the user may change these by tweaking the appropriate parameters in `exe/generateL0seed.cpp`. The output of the second line is stored as a set of estimated distances and estimated runtimes in microseconds, stored under the `results/` folder.

# Datasets

All our datasets are publicly available networks taken from repositories on the internet. We list the datasets used, and point the reader to the appropriate website for downloading these files. The files are from one of two sources: the SNAP large networks database [[1]](#1), and the KONECT project [[2]](#2). The following table has the links to the appropriate pages.

<div align="center">

| _Network_       | _Source_                                                     |
| --------------- | ------------------------------------------------------------ |
| caida           | [KONECT](http://konect.cc/networks/as-caida20071105/)        |
| epinions        | [SNAP](https://snap.stanford.edu/data/soc-Epinions1.html)    |
| slashdot        | [SNAP](https://snap.stanford.edu/data/soc-Slashdot0811.html) |
| dblp            | [SNAP](https://snap.stanford.edu/data/com-DBLP.html)         |
| higgs-twitter   | [SNAP](https://snap.stanford.edu/data/higgs-twitter.html)    |
| skitter         | [SNAP](https://snap.stanford.edu/data/as-Skitter.html)       |
| large-dblp      | [KONECT](http://konect.cc/networks/dblp_coauthor/)           |
| soc-pokec       | [SNAP](https://snap.stanford.edu/data/soc-Pokec.html)        |
| soc-livejournal | [SNAP](https://snap.stanford.edu/data/soc-LiveJournal1.html) |
| soc-orkut       | [SNAP](https://snap.stanford.edu/data/com-Orkut.html)        |
| wikipedia       | [KONECT](http://konect.cc/networks/wikipedia_link_en/)       |
| soc-twitter     | [SNAP](https://snap.stanford.edu/data/twitter-2010.html)     |

</div>

# References

<a id="1">[1]</a>
Jure Leskovec and Andrej Krevl. 2014. SNAP Datasets: Stanford Large
Network Dataset Collection. http://snap.stanford.edu/data.

<a id="2">[2]</a>
Jérôme Kunegis. 2013. KONECT: The Koblenz Network Collection. In
_Proceedings of the 22nd International Conference on World Wide Web (WWW ’13 Companion)_. Association for Computing
Machinery, 1343–1350.
