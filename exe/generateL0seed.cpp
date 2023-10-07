#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Graph.h"
#include "Escape/L0Graph.h"
#include "Escape/IndexedBinaryHeap.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <chrono>

// exe/generateL0seed [graph_name] [max_seed/seed]
using namespace Escape;
using namespace std::chrono;
int main(int argc, char *argv[]) {
    std::string  graph_name = argv[1];
    std::string init_type = argv[2];

    std::string delimiter = "_";
    std::size_t d_pos = init_type.find(delimiter);
    bool high_degree_seed = d_pos != std::string::npos;

    std::string graph_filename = graph_name + ".edges";
    std::string graph_folder = "graphs/";
    std::string graph_loc = graph_folder + graph_filename;
    std::string input_folder = "inputs/";
    std::string results_folder = "results/";
    std::string type = "_" + init_type;
    std::ofstream graph_file(results_folder + graph_name + type + ".txt");

    // Graph g;
    // if (loadGraph(graph_loc.c_str(), g, 1, IOFormat::escape))
    // exit(1);
    // printf("Loaded graph\n");
    // CGraph cg = makeCSR(g);
    // cg.sortById();
    // printf("Converted to CSR\n");
    CGraph cg;
    cg.loadGraphFromFile(graph_name);
    std::cout<<cg.nEdges<<std::endl;
    std::cout<<cg.nVertices<<std::endl;
    
    VertexIdx vcount = cg.nVertices;
    int numL0s = 4;
    VertexIdx* L0_sizes = new VertexIdx[numL0s];
    for(int s = 0; s < numL0s; s++){
        float percent = (float)s*0.005;
        L0_sizes[s] = (VertexIdx) (percent*vcount);
        printf("L0 size: %ld %f\n", L0_sizes[s], percent);
    }
    for(int L0_idx = 1; L0_idx < numL0s; L0_idx++){
        auto start = std::chrono::high_resolution_clock::now();
        L0Graph L0 = L0Graph(cg, L0_sizes[L0_idx], high_degree_seed);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = duration_cast<std::chrono::microseconds>(end - start);
        double duration_count = (double) duration.count();
        float percent = L0_idx*0.5;
        std::string num_text = "_" + std::to_string(percent);
        std::string L0_size_str = num_text.substr(0, num_text.find(".")+2);
        L0.checkForBadL0();
        L0.writeGraphToFile(graph_name + type  + L0_size_str);
        L0.print_size(graph_file, graph_name + type  + L0_size_str);
        graph_file << "initialization time for " << percent << "%: " <<duration_count <<"\n";
    }
    graph_file.close();
    return 0;
}