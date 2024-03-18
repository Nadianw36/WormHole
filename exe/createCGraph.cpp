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
  std::string graph_name = argv[1];
  // printf("hiii\n");
  checkSetupFor(graph_name);
  std::string graph_filename = graph_name + ".edges";
  std::string graph_loc = GRAPH_FOLDER + graph_filename;
  Graph g;
  if (loadGraph(graph_loc.c_str(), g, 1, IOFormat::escape))
    exit(1);
  // printf("Loaded graph\n");
  CGraph cg = makeCSR(g);
  cg.sortById();
  // printf("Converted to CSR\n");

  cg.saveGraphToFile(graph_name);

  return 0;
}