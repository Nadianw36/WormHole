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
    const char *exeDescription =
      "Arguments: \n"
      "-g   Graph name. The input CSR Graph should be named {graph_name}_CGraph.bin \n"
      "-gi  Graph input relative location. The edgelist should be found in the folder "
      "/bin/{input_loc}/{graph_name}_CGraph.bin. L0 bin\n"
      "-go  Core output relative location. The core boolean lists will be found in "
      "/bin/{graph_output_loc}/L0 .\n"
      "-o   Output relative location. The files contianing core statistics will be"
      "found in /results/{output_loc}/L0."
      "Any non existing subfolders will be created.\n";

    if (foundHelpFlag(argc, argv, exeDescription))
        return 0;
    std::string graph_name = getCmdOption(argc, argv, "-g", "need to provide graph name after flag -g");
    std::string output_loc = getCmdOption(argc, argv, "-o", "need to provide output subfolder after flag -o");
    std::string graph_input_loc = getCmdOption(argc, argv, "-gi", "need to provide graph input subfolder after flag -gi");
    std::string graph_output_loc = getCmdOption(argc, argv, "-go", "need to provide graph output subfolder after flag -go");

    string graph_results_folder = RESULTS_FOLDER + graph_output_loc + "/";
    ofstream graph_file(graph_results_folder + graph_name  + "_core-init-times.txt");

    checkL0SetupFor(graph_name);
    std::cout << std::setprecision(3);

    CGraph cg;
    cg.loadGraphFromFile(graph_name);
    cout << cg.nVertices << " " << cg.nEdges << endl;

    VertexIdx vcount = cg.nVertices;

    int numL0s = 11;
    float increment = 0.01;
    VertexIdx *L0_sizes = new VertexIdx[numL0s];
    for (int s = 0; s < numL0s; s++)
    {
        float percent = (float)s * increment;
        L0_sizes[s] = (VertexIdx)(percent * vcount);
        graph_file << "core size " << percent << "%: " << L0_sizes[s] << " vertices\n";
        cout << "core size " << percent << "%: " << L0_sizes[s] << " vertices\n";
    }

    for (int L0_idx = 1; L0_idx < numL0s; L0_idx++)
    {
        auto start = chrono::high_resolution_clock::now();
        L0Graph L0 = L0Graph(cg, L0_sizes[L0_idx], false);
        auto end = chrono::high_resolution_clock::now();
        auto duration = duration_cast<chrono::nanoseconds>(end - start);
        long long int duration_count = (double)duration.count();
        float percent = L0_idx * (increment * 100);
        string num_text = "_" + to_string(percent);
        string L0_size_str = num_text.substr(0, num_text.find(".") + 2);
        L0.checkForBadL0();
        L0.writeGraphToFile(graph_output_loc  + "/" + graph_name + L0_size_str);
        L0.print_size(graph_file, graph_name + L0_size_str);
        graph_file << "initialization time for " << percent << "%: " << duration_count << " nanoseconds\n";
        std::cout << "initialization time for " << percent << "%: " << duration_count << " nanoseconds\n";
    }
    graph_file.close();
    return 0;
}