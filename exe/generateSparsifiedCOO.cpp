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
    std::string graph_name = argv[1];

    checkL0SetupFor(graph_name);

    CGraph cg;
    cg.loadGraphFromFile(graph_name);
    cout << cg.nEdges << endl;
    cout << cg.nVertices << endl;

    const char* colour[20]
        = { "0-5", "1-0", "1-5", "2-0","2-5", "3-0", "3-5","4-0", "4-5","5-0", "5-5","6-0", "6-5","7-0", "7-5",
        "8-0", "8-5", "9-0", "9-5", "10-0"};
    for (int i = 0; i < 20; i++)
    {
        std::string graph_L0_name = graph_name + "_seed_" + colour[i];
        L0Graph L0 = L0Graph(cg, graph_L0_name);
        L0.checkForBadL0();
        L0.writeSparsifiedC00(graph_L0_name);
    } 
    

    // L0.checkForBadL0();
    // float sparsified[] = {0.3, 0.4, 0.5, 0.6, 0.7};
    // for (int i = 0; i < 5; i++)
    // {
    //     L0.writeSparsifiedC00(graph_L0_name, sparsified[i], prune);
    // }

    return 0;
}