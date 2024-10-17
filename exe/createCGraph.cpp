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

// createCGraph -o [output] -i [input_folder] -g [graph_name]
// createCGraph -h
int main(int argc, char *argv[])
{
  const char *exeDescription =
      "Arguments: \n"
      "-g   Graph name. The input edgelist should be named {graph_name}.edges \n"
      "-i   Input relative location. The edgelist should be found in the folder "
      "/graphs/{input_loc}. \n"
      "-o   Output relative location. The binary file containing the CGraph will "
      "found in /bin/{output_loc}/{graph_name}_CGraph.bin. Any non existing subfolders will be created.";

  if (foundHelpFlag(argc, argv, exeDescription))
    return 0;
  std::string graph_name = getCmdOption(argc, argv, "-g", "need to provide graph name after flag -g");
  std::string input_loc = getCmdOption(argc, argv, "-i", "need to provide input subfolder after flag -i");
  std::string output_loc = getCmdOption(argc, argv, "-o", "need to provide output subfolder after flag -o");

  std::string graph_filename = graph_name + ".edges";
  std::string graph_loc = GRAPH_FOLDER + input_loc + "/" + graph_filename;

  Graph g;
  if (loadGraph(graph_loc.c_str(), g, 1, IOFormat::escape))
    exit(1);
  printf("Loaded graph\n");
  CGraph cg = makeCSR(g);
  cg.sortById();
  printf("Converted to CSR\n");
  cg.saveGraphToFile(output_loc, graph_name);
  return 0;
}