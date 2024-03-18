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
using namespace std::chrono;

// exe/L0 [action] [graphname] [L0_size (ie seed_{L0})] [start_from] [end_at]
int main(int argc, char *argv[])
{
  std::string graph_name = argv[1];
  std::string results_folder = RESULTS_FOLDER + graph_name + "/";
  std::ofstream graph_file(results_folder + graph_name + "_PLL" + ".txt");

  CGraph cg;
  cg.loadGraphFromFile(graph_name);
  std::vector<std::pair<int, int>> edgelist;

  for (VertexIdx i = 0; i < cg.nVertices; i++)
  {
    for (EdgeIdx e = cg.offsets[i]; e < cg.offsets[i + 1]; e++)
    {
      VertexIdx nbor = cg.nbors[e];
      std::pair<int, int> edge((int64_t)i, (int64_t)nbor);
      edgelist.push_back(edge);
    }
  }
  // printf("Loaded graph\n");
  // printf("creating pll graph\n");
  auto start_init = high_resolution_clock::now();
  PrunedLandmarkLabeling<> pll;
  pll.ConstructIndex(edgelist);
  auto end_init = high_resolution_clock::now();
  auto duration_init = duration_cast<microseconds>(end_init - start_init);
  graph_file << "initializing time: " << (double)(duration_init.count()) << "\n";

  pll.PrintSize(graph_file);

  // check how experiments should run; don't include a third argument to run normally (without picking samples beforehand)
  int startFrom = atoi(argv[2]);
  int ROUNDS = atoi(argv[3]);

  // initialize varibles for sampling
  PrunedLandmarkLabeling<> pll_graph;

  std::vector<int> actualDistances;
  std::vector<int> actualRuntimes;

  // make files
  // printf("created variables\n");

  std::fstream nodesTouchedFile;
  std::fstream nodesTouchedWithoutSamplesFile;

  std::ofstream queriesFile;
  std::ofstream queriesAccumulatedWithoutSampleFile;
  std::ofstream queriesAccumulatedFile;
  std::ofstream estimatedDistancesFile;
  std::ofstream estimateRuntimesFile;
  std::ofstream estimateSampleFile;

  std::ofstream actualDistancesFile;
  std::ofstream actualRuntimesFile;
  std::ofstream actualSampleFile;

  actualDistances.clear();
  actualRuntimes.clear();

  actualDistancesFile.open(RESULTS_FOLDER + graph_name + "_PLL_distances.txt");
  actualRuntimesFile.open(RESULTS_FOLDER + graph_name + "_PLL_runtimes.txt");
  actualSampleFile.open(RESULTS_FOLDER + graph_name + "_PLL_failed.txt");

  std::ifstream inputFile;
  inputFile.open(INPUT_FOLDER + graph_name + "_input.txt");
  // printf("starting rounds\n");
  for (int round = 0; round < ROUNDS; round++)
  {
    VertexIdx v1;
    VertexIdx v2;
    inputFile >> v1 >> v2;
    if (round < startFrom)
      continue;

    if (round % 1000 == 0)
    {
      // printf("round %i\n", round);
      for (const auto &e : actualDistances)
        actualDistancesFile << e << "\n";
      for (const auto &e : actualRuntimes)
        actualRuntimesFile << e << "\n";

      actualDistances.clear();
      actualRuntimes.clear();
      // printf("uploaded through round %i\n", round - 1);
    }

    if (v1 == v2)
      continue;

    std::chrono::high_resolution_clock::time_point start_actual;
    int actual;
    // int actual = pll_graph.QueryDistance(v1, v2);
    std::chrono::high_resolution_clock::time_point end_actual;
    start_actual = std::chrono::high_resolution_clock::now();
    actual = pll.QueryDistance(v1, v2);
    end_actual = std::chrono::high_resolution_clock::now();
    if (actual == -1 || actual >= cg.nEdges)
    {
      actualSampleFile << (int64_t)v1 << " " << (int64_t)v2 << "\n";
      actual = -1;
    }
    actualDistances.push_back((double)(actual));
    auto duration_actual = duration_cast<std::chrono::microseconds>(end_actual - start_actual);
    actualRuntimes.push_back((double)(duration_actual.count()));

    std::chrono::high_resolution_clock::time_point start_estimate;
    SamplingInformation e;
    std::chrono::high_resolution_clock::time_point end_estimate;
  }
  for (const auto &e : actualDistances)
    actualDistancesFile << e << "\n";
  for (const auto &e : actualRuntimes)
    actualRuntimesFile << e << "\n";

  actualDistancesFile.close();
  actualRuntimesFile.close();
  actualSampleFile.close();
  return 0;
}