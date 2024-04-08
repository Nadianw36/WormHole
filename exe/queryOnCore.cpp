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
string getInputName(string command, string graphName, string graphL0Name)
{

  string inputFolder = getSubfolderName("L0-BiBFS");
  string inputPath = RESULTS_FOLDER + graphName + "/" + inputFolder;

  string inputFilenamePrefix = inputPath + "0_10000_" + graphL0Name + "_L0-BiBFS_";

  if (!command.compare("L0-BiBFS-RandomL0"))
    return inputFilenamePrefix + "core-queries_random-L0.txt";
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

  checkSetupFor(graph_name);

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

  string types[2] = {"L0-BiBFS-RandomL0", "L0-BiBFS-RandomL0FromL1"};
  for (int i = 0; i < 2; i++)
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

    std::vector<int> numQueries;
    std::vector<int> distances;
    std::vector<int> runtimes;

    std::ofstream distancesFile(results_path + "_distances.txt");
    std::ofstream pathsFile(results_path + "_paths.txt");
    std::ofstream runtimesFile(results_path + "_runtimes.txt");
    std::ofstream queriesFile(results_path + "_queries.txt");

    std::ifstream inputFile;
    inputFile.open(getInputName(command, graph_name, graph_L0_name));
    std::cout << getInputName(command, graph_name, graph_L0_name) << std::endl;

    // printf("starting round\n");
    std::chrono::high_resolution_clock::time_point start_clock;
    std::chrono::high_resolution_clock::time_point end_clock;

    for (int round = 0; round < ROUNDS; round++)
    {
      VertexIdx v1;
      VertexIdx v2;
      inputFile >> v1 >> v2;

      if (round < startFrom || v1 == v2)
        continue;
      // printf("query %ld %d %ld %d \n", v1, L0.L0[v1], v2, L0.L0[v2]);
      if (round % 1000 == 0 && round > startFrom)
      {
        // printf("round %i\n", round);

        for (const auto &e : distances)
          distancesFile << e << "\n";
        for (const auto &e : runtimes)
          runtimesFile << e << "\n";
        for (const auto &e : numQueries)
          queriesFile << e << "\n";

        distances.clear();
        runtimes.clear();
        numQueries.clear();
      }

      if (v1 == v2)
        continue;

      if (v1 < 0 || v2 < 0)
      { // doesn't traverse through core
        distances.push_back(-1),
            runtimes.push_back(-1);
        numQueries.push_back(-1);
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

      distances.push_back((double)(distance));
      auto duration = duration_cast<std::chrono::microseconds>(end_clock - start_clock);
      runtimes.push_back((double)(duration.count()));

      visited[v1] = false;
      visited[v2] = false;

      numQueries.push_back(countQueries(visited, bg.nVertices));
    }
    // printf("ending rounds\n");
    for (const auto &e : distances)
      distancesFile << e << "\n";
    for (const auto &e : runtimes)
      runtimesFile << e << "\n";
    for (const auto &e : numQueries)
      queriesFile << e << "\n";
    distancesFile.close();
    pathsFile.close();
    runtimesFile.close();
  }
  return 0;
}