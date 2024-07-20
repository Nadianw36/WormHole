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
    std::string graph_L0_name = argv[1];
    bool prune = true;

    std::string delimiter = "_";
    std::size_t d_pos = graph_L0_name.find(delimiter);
    if (d_pos != std::string::npos)
        graph_name = graph_L0_name.substr(0, d_pos);

    checkSetupFor(graph_name);

    CGraph cg;
    cg.loadGraphFromFile(graph_name);
    cout << cg.nEdges << endl;
    cout << cg.nVertices << endl;

    L0Graph L0 = L0Graph(cg, graph_L0_name);
    L0.checkForBadL0();
    L0.writeSparsifiedC00(graph_L0_name);

    // L0.checkForBadL0();
    // float sparsified[] = {0.3, 0.4, 0.5, 0.6, 0.7};
    // for (int i = 0; i < 5; i++)
    // {
    //     L0.writeSparsifiedC00(graph_L0_name, sparsified[i], prune);
    // }

    return 0;
}