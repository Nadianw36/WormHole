#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Graph.h"
#include "Escape/L0Graph.h"
#include "Escape/Config.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <chrono>

using namespace Escape;
using namespace std::chrono;
int main(int argc, char *argv[])
{
    std::string graph_name = argv[1];
    std::string graph_filename = graph_name + ".edges";
    std::string graph_loc = GRAPH_FOLDER + graph_filename;
    std::ofstream graph_file(RESULTS_FOLDER + graph_name + ".txt");

    Graph g;
    if (loadGraph(graph_loc.c_str(), g, 1, IOFormat::escape))
        exit(1);
    // printf("Loaded graph\n");
    CGraph cg = makeCSR(g);
    cg.sortById();
    // printf("Converted to CSR\n");
    std::cout << cg.nEdges << std::endl;
    std::cout << cg.nVertices << std::endl;

    VertexIdx vcount = cg.nVertices;
    int num_sizes = 4;
    VertexIdx *L0_sizes = new VertexIdx[num_sizes];
    for (int s = 0; s < num_sizes; s++)
    {
        float percent = (float)s * 0.005;
        L0_sizes[s] = (VertexIdx)(percent * vcount);
        // printf("L0 size: %ld %f\n", L0_sizes[s], percent);
    }
    // start_timer
    auto start = std::chrono::high_resolution_clock::now();

    std::priority_queue<std::pair<EdgeIdx, VertexIdx>> p;
    bool *L0 = new bool[vcount];
    bool *L1 = new bool[vcount];
    std::fill_n(L0, vcount, false);
    std::fill_n(L1, vcount, false);
    for (VertexIdx v = 0; v < vcount; v++)
    {
        EdgeIdx degrees = cg.degree(v);
        p.push(std::make_pair(degrees, v));
    }
    // end timer 1
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - start);
    double duration_count = (double)duration.count();
    for (int L0_idx = 1; L0_idx < 6; L0_idx++)
    {
        start = std::chrono::high_resolution_clock::now();
        for (VertexIdx i = L0_sizes[L0_idx - 1]; i < L0_sizes[L0_idx]; i++)
        {
            // printf("%ld \n", i);
            VertexIdx vl0 = p.top().second;
            // printf("%ld \n", vl0);
            p.pop();
            L0[vl0] = true;
            L1[vl0] = false;
            // printf("    added to L0\n", vl0);
            for (EdgeIdx e = cg.offsets[vl0]; e < cg.offsets[vl0 + 1]; e++)
            {
                VertexIdx nbor = cg.nbors[e];
                if (!L0[nbor])
                {
                    L1[nbor] = true;
                }
            }
            // printf("    added neighbors to L1\n", vl0);
        }
        auto end = std::chrono::high_resolution_clock::now();
        duration = duration_cast<std::chrono::microseconds>(end - start);
        duration_count += (double)duration.count();
        float percent = L0_idx * 0.5;
        std::string L0_size_str = "_" + std::to_string(percent);
        std::ofstream L0File(BIN_FOLDER + graph_name + L0_size_str + "_L0.bin", std::ios::out | std::ios::binary);
        std::ofstream L1File(BIN_FOLDER + graph_name + L0_size_str + "_L1.bin", std::ios::out | std::ios::binary);
        L0File.write((char *)L0, sizeof(bool) * (int64_t)vcount);
        L1File.write((char *)L1, sizeof(bool) * (int64_t)vcount);
        L0File.close();
        L1File.close();
        graph_file << "initialization time for " << percent << "%: " << duration_count << "\n";
    }

    for (int L0_idx = 1; L0_idx < num_sizes; L0_idx++)
    {
        // printf("pruning L0\n");
        float percent = L0_idx * 0.5;
        std::string L0_size_str = "_" + std::to_string(percent);
        L0Graph L0 = L0Graph(cg, graph_name + L0_size_str);
        VertexIdx size = L0.connectComponentsL0();
        L0.writeGraphToFile(graph_name + L0_size_str);
        L0.print_size(graph_file, graph_name + L0_size_str);
    }
    graph_file.close();
    return 0;
    // run connected components after generation
}