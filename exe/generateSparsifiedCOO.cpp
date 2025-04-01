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

// exe/generateL0COO graph_name [-b bin_path] [-o output_path] [-WH]

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
      "-gCore  Core input relative location. The core boolean lists will be found in "
      "/bin/{graph_output_loc}/L0 .\n"
      "-o   Output relative location. The files contianing core statistics will be"
      "found in /graphs/{output_loc}. "
      "Edges L0->L0 and L0->L1 will be in {graph_name}_inter-core-COO.txt and"
      " all other edges will be in {graph_name}_outside-core-COO.txt"
      "Any non existing subfolders will be created.\n";

    if (foundHelpFlag(argc, argv, exeDescription))
        return 0;
    std::string graph_name = getCmdOption(argc, argv, "-g", "need to provide graph name after flag -g");
    std::string output_loc = getCmdOption(argc, argv, "-o", "need to provide output subfolder after flag -o");
    std::string graph_input_loc = getCmdOption(argc, argv, "-gi", "need to provide graph input subfolder after flag -gi");
    std::string graph_output_loc = getCmdOption(argc, argv, "-gCore", "need to provide Core input relative location -gCore");
    std::string graph_name = argv[1];

    // checkL0SetupFor(graph_name);

    CGraph cg;
    cg.loadGraphFromFile(graph_input_loc, graph_name);
    cout << cg.nEdges << endl;
    cout << cg.nVertices << endl;

    const char* colour[20]
        = { "0-5", "1-0", "1-5", "2-0","2-5", "3-0", "3-5","4-0", "4-5","5-0", "5-5","6-0", "6-5","7-0", "7-5",
        "8-0", "8-5", "9-0", "9-5", "10-0"};
    for (int i = 0; i < 20; i++)
    {
        std::string graph_L0_name = graph_name + "_" + colour[i];
        L0Graph L0 = L0Graph(cg, graph_L0_name);
        L0.checkForBadL0();
        L0.writeSparsifiedC00(graph_output_loc + "/" + graph_L0_name);
    } 
    

    // L0.checkForBadL0();
    // float sparsified[] = {0.3, 0.4, 0.5, 0.6, 0.7};
    // for (int i = 0; i < 5; i++)
    // {
    //     L0.writeSparsifiedC00(graph_L0_name, sparsified[i], prune);
    // }

    return 0;
}