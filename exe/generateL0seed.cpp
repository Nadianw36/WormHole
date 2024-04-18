#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Graph.h"
#include "Escape/L0Graph.h"
#include "Escape/IndexedBinaryHeap.h"
#include "Escape/Config.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <chrono>

// exe/generateL0seed [graph_name] [max_seed/seed]
using namespace Escape;
using namespace std;
using namespace chrono;
int main(int argc, char *argv[])
{
    string graph_name = argv[1];
    string init_type = argv[2];

    string delimiter = "_";
    size_t d_pos = init_type.find(delimiter);
    bool high_degree_seed = d_pos != string::npos;

    string type = "_" + init_type;
    string graph_results_folder = RESULTS_FOLDER + graph_name + "/" + L0_FOLDER;
    ofstream graph_file(graph_results_folder + graph_name + type + ".txt");

    checkSetupFor(graph_name);
    std::cout << std::setprecision(3);

    CGraph cg;
    cg.loadGraphFromFile(graph_name);
    cout << cg.nVertices << " " << cg.nEdges << endl;

    VertexIdx vcount = cg.nVertices;

    int numL0s = 21;
    VertexIdx *L0_sizes = new VertexIdx[numL0s];
    for (int s = 0; s < numL0s; s++)
    {
        float percent = (float)s * 0.005;
        L0_sizes[s] = (VertexIdx)(percent * vcount);
        graph_file << "core size " << percent << "%: " << L0_sizes[s] << " vertices\n";
        cout << "core size " << percent << "%: " << L0_sizes[s] << " vertices\n";
    }

    for (int L0_idx = 1; L0_idx < numL0s; L0_idx++)
    {
        auto start = chrono::high_resolution_clock::now();
        L0Graph L0 = L0Graph(cg, L0_sizes[L0_idx], high_degree_seed);
        auto end = chrono::high_resolution_clock::now();
        auto duration = duration_cast<chrono::nanoseconds>(end - start);
        long long int duration_count = (double)duration.count();
        float percent = L0_idx * 0.5;
        string num_text = "_" + to_string(percent);
        string L0_size_str = num_text.substr(0, num_text.find(".") + 2);
        L0.checkForBadL0();
        L0.writeGraphToFile(graph_name + type + L0_size_str);
        L0.print_size(graph_file, graph_name + type + L0_size_str);
        graph_file << "initialization time for " << percent << "%: " << duration_count << " nanoseconds\n";
        std::cout << "initialization time for " << percent << "%: " << duration_count << " nanoseconds\n";
    }
    graph_file.close();
    return 0;
}