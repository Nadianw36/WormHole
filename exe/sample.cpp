#include "Escape/GraphIO.h"
#include "Escape/Graph.h"
#include "Escape/BiBFSGraph.h"
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

int countQueries(bool *visited, VertexIdx nVertices) { return std::count(visited, visited + (int64_t)nVertices, true); };

// exe/filname [graphname/graphname_seed_{n}] [start_from] [end_at]
// note: graphname cannot contain any instances of _
int main(int argc, char *argv[])
{
  const char *exeDescription =
      "Arguments: \n"
      "-g   Graph name. The input CSR Graph should be named {graph_name}_CGraph.bin \n"
      "-gi  Graph input relative location. The edgelist should be found in the folder "
      "/bin/{input_loc}/{graph_name}_CGraph.bin. \n"
      "-i   The input filename and relative location. The queries list should be  "
      "/inputs/{input_loc}.\n"
      "-o   Output relative location. The file containing the results will be "
      "found in /results/{output_loc}/{graph_name}_distances.txt."
      "Any non existing subfolders will be created.\n"
      "-s   Optional. Start query number, zero indexed, inclusive. Default is 0.\n"
      "-t   Optional. End query number, zero indexed, exclusive. Default is 10,000.\n";

  if (foundHelpFlag(argc, argv, exeDescription))
    return 0;
  std::string graph_name = getCmdOption(argc, argv, "-g", "need to provide graph name after flag -g");
  std::string input_loc = getCmdOption(argc, argv, "-i", "need to provide input subfolder after flag -i");
  std::string output_loc = getCmdOption(argc, argv, "-o", "need to provide output subfolder after flag -o");
  std::string graph_input_loc = getCmdOption(argc, argv, "-gi", "need to provide graph input subfolder after flag -gi");

  char *char_start = getOptionalCmdOption(argc, argv, "-s");
  char *char_end = getOptionalCmdOption(argc, argv, "-t");
  string str_start = (char_start != nullptr) ? char_start : "0";
  string str_end = (char_end != nullptr) ? char_end : "10000";

  std::string results_folder = RESULTS_FOLDER + output_loc + "/";
  create_directories(results_folder);

  std::string results_prefix = str_start + "_" + str_end + "_";
  std::string results_path = results_folder + results_prefix + graph_name;

  CGraph cg;
  cg.loadGraphFromFile(graph_input_loc, graph_name);

  BiBFSGraph bg = BiBFSGraph(cg);

  int startFrom = stoi(str_start);
  int ROUNDS = stoi(str_end);

  std::ofstream distancesFile(results_path + "_distances.txt");
  std::ofstream failedSampleFile(results_path + "_failed.txt");

  std::ifstream inputFile(INPUT_FOLDER + input_loc);
  std::cout << INPUT_FOLDER + input_loc << endl;
  printf("starting round\n");
  VertexIdx v1, v2;
  for (int round = 0; round < ROUNDS; round++)
  {
    inputFile >> v1 >> v2;
    // m
    if (round < startFrom || v1 == v2)
      continue;

    int distance = bg.BidirectionalBFS(v1, v2);

    if (distance == -1)
    {
      distancesFile << -1 << std::endl;
      failedSampleFile << (int64_t)v1 << " " << (int64_t)v2 << " " << round << "\n";
      continue;
    }

    distancesFile << distance << "\n";
  }
  distancesFile.close();

  return 0;
}