#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Graph.h"
#include "Escape/L0Graph.h"
#include "Escape/Config.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <chrono>
#include <string>
#include <cstring>
#include <fstream>

using namespace Escape;
using namespace std::chrono;

// exe/filname [BiBFS/L0] [graphname/graphname_seed_{n}] [start_from] [end_at]
// note: graphname cannot contain any instances of _
int main(int argc, char *argv[])
{

  std::string graph_name = argv[1];
  std::string graph_L0_name = argv[1];

  std::string delimiter = "_";
  std::size_t d_pos = graph_L0_name.find(delimiter);
  if (d_pos != std::string::npos)
    graph_name = graph_L0_name.substr(0, d_pos);

  std::ofstream graph_file(RESULTS_FOLDER + graph_L0_name + ".txt");

  CGraph cg;
  cg.loadGraphFromFile(graph_name);
  // printf("loading L0 and L1\n");
  L0Graph L0;
  L0 = L0Graph(cg, graph_L0_name);
  L0.print_size(graph_file, graph_L0_name);
  return 0;
}