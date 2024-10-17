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
string getInputName(string command, string graphName, string graphL0Name, string start, string finish)
{

  string inputFolder = getSubfolderName("L0-BiBFS");
  string inputPath = RESULTS_FOLDER + graphName + "/" + inputFolder;

  string inputFilenamePrefix = inputPath + start + "_" + finish + "_" + graphL0Name + "_L0-BiBFS_";

  if (!command.compare("L0-BiBFS-RandomL0"))
    return inputFilenamePrefix + "core-queries_random-L0.txt";
  else if (!command.compare("L0-BiBFS-HighestDegL0"))
    return inputFilenamePrefix + "core-queries_highest-deg-L0.txt";
  else
    return inputFilenamePrefix + "core-queries_random-L0-random-L1.txt";
}

// exe/queryOnCore [graphname_seed_{size}]
// note: graphname cannot contain any instances of _
int main(int argc, char *argv[])
{
  std::string graph_name = argv[1];
  std::string graph_L0_name = argv[1];

  std::string delimiter = "_";
  std::size_t d_pos = graph_L0_name.find(delimiter);

  checkL0SetupFor(graph_name);

  std::string str_start = argv[2];
  std::string str_finish = argv[3];

  if (d_pos != std::string::npos)
    graph_name = graph_L0_name.substr(0, d_pos);

  CGraph cg;
  cg.loadGraphFromFile(graph_name);
  BiBFSGraph bg = BiBFSGraph(cg);
  int startFrom = stoi(str_start);
  int ROUNDS = stoi(str_finish);

  L0Graph L0 = L0Graph(cg, graph_L0_name);
  L0.checkForBadL0();

  string types[3] = {"L0-BiBFS-RandomL0", "L0-BiBFS-RandomL0FromRandomL1", "L0-BiBFS-HighestDegL0"};
  for (int i = 0; i < 3; i++)
  {

    std::string command = types[i];

    std::string results_graph_sub_folder = getSubfolderName(command);
    std::cout << results_graph_sub_folder << std::endl;

    std::string results_folder = RESULTS_FOLDER + graph_name + "/" + results_graph_sub_folder;

    std::string results_prefix = str_start + "_" + str_finish + "_";
    std::string results_suffix = "_" + command;
    std::string results_path = results_folder + results_prefix + graph_L0_name + results_suffix;

    std::ofstream graph_file(results_folder + graph_name + results_suffix + ".txt");

    VertexIdx vcount = cg.nVertices;

    std::ofstream distancesFile(results_path + "_distances.txt");
    std::ofstream pathsFile(results_path + "_paths.txt");
    std::ofstream runtimesFile(results_path + "_runtimes_nano.txt");
    std::ofstream queriesFile(results_path + "_queries.txt");

    std::string inputFilename = getInputName(command, graph_name, graph_L0_name, str_start, str_finish);
    std::ifstream inputFile(inputFilename);
    std::cout << "reading queries from " << inputFilename << std::endl;

    // printf("starting round\n");
    std::chrono::high_resolution_clock::time_point start_clock;
    std::chrono::high_resolution_clock::time_point end_clock;
    VertexIdx v1, v2;
    int totalRuntimeThroughCore = 0;

    while (inputFile >> v1 >> v2)
    {
      if (v1 < 0 || v2 < 0)
      { // doesn't traverse through core
        distancesFile << -1 << "\n";
        runtimesFile << -1 << "\n";
        queriesFile << -1 << "\n";
        pathsFile << "\n";
        continue;
      }

      L0.reset();
      start_clock = std::chrono::high_resolution_clock::now();
      int distance = L0.distance_L0_to_L0(v1, v2); // for now don't run with PLL
      end_clock = std::chrono::high_resolution_clock::now();
      bool *visited = L0.getVisited();

      vector<VertexIdx> path = L0.recoverPath(v1, v2);

      for (const auto &p : path)
        pathsFile << p << " ";
      pathsFile << "\n";

      distancesFile << distance << "\n";
      auto duration = duration_cast<std::chrono::nanoseconds>(end_clock - start_clock);
      runtimesFile << (long long int)(duration.count()) << "\n";
      totalRuntimeThroughCore += (long long int)(duration.count());

      visited[v1] = false;
      visited[v2] = false;

      queriesFile << countQueries(visited, bg.nVertices) << "\n";
    }
    printf("For queries sampled as %s, all queries on core ran in %lld nanoseconds\n", command.c_str(), totalRuntimeThroughCore);
    distancesFile.close();
    pathsFile.close();
    runtimesFile.close();
    queriesFile.close();
  }
  return 0;
}