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

using namespace Escape;

// generat
int main(int argc, char *argv[])
{
  EdgeIdx new_nVertices = 1275997901;
  std::string graph_name = argv[1];
  std::string L0_size = argv[2];
  std::string graph_L0_name = graph_name + "_seed_" + L0_size;
  std::string graph_filename = graph_name + ".edges";
  std::string graph_loc = GRAPH_FOLDER + graph_filename;
  std::ofstream graph_file(RESULTS_FOLDER + graph_name + "L0_optimizedd" + ".txt");

  CGraph cg;
  cg.loadGraphFromFile(graph_name);
  L0Graph L0;
  // printf("loading L0 and L1\n");
  L0 = L0Graph(cg, graph_L0_name);
  cg.sanitizeForL0(graph_L0_name, L0.L0, L0.L1, 1275997901);

  cg.saveGraphToFile(graph_L0_name);
  // printf("has it update? %ld \n", cg.nEdges);
  L0.print_size(graph_file, graph_L0_name + "_optimized");

  return 0;
}