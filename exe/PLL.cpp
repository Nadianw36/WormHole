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

// exe/PLL [graph_name] [start] [finish]
int main(int argc, char *argv[])
{
  std::string graph_name = argv[1];
  checkL0SetupFor(graph_name);

  std::string results_graph_sub_folder = RESULTS_FOLDER + graph_name + "/" + getSubfolderName("PLL");
  std::cout << results_graph_sub_folder << std::endl;

  std::ofstream graph_file(results_graph_sub_folder + graph_name + "_PLL" + ".txt");

  CGraph cg;
  cg.loadGraphFromFile(graph_name);
  std::vector<std::pair<int, int>> edgelist;
  printf("creating edgelist\n");
  for (VertexIdx i = 0; i < cg.nVertices; i++)
  {
    for (EdgeIdx e = cg.offsets[i]; e < cg.offsets[i + 1]; e++)
    {
      VertexIdx nbor = cg.nbors[e];
      std::pair<int, int> edge((int64_t)i, (int64_t)nbor);
      edgelist.push_back(edge);
    }
  }
  printf("Loaded graph\n");
  printf("creating pll graph\n");
  auto start_init = high_resolution_clock::now();
  PrunedLandmarkLabeling<> pll;
  pll.ConstructIndex(edgelist);
  auto end_init = high_resolution_clock::now();
  auto duration_init = duration_cast<nanoseconds>(end_init - start_init);

  graph_file << "finished creating graph in: " << (long long int)(duration_init.count()) << " nanoseconds\n";
  printf("finished creating graph in: %lld nanoseconds\n", (long long int)(duration_init.count()));

  pll.PrintSize(graph_file);

  int startFrom = atoi(argv[2]);
  int ROUNDS = atoi(argv[3]);

  std::ofstream actualDistancesFile(results_graph_sub_folder + graph_name + "_PLL_distances.txt");
  std::ofstream actualRuntimesFile(results_graph_sub_folder + graph_name + "_PLL_runtimes_nano.txt");
  std::ofstream actualSampleFile(results_graph_sub_folder + graph_name + "_PLL_failed.txt");

  std::ifstream inputFile(INPUT_FOLDER + graph_name + "_input.txt");

  long long totalRuntime = 0;
  int totalRounds = 0;
  printf("starting rounds\n");
  for (int round = 0; round < ROUNDS; round++)
  {
    VertexIdx v1;
    VertexIdx v2;
    inputFile >> v1 >> v2;
    if (round < startFrom || v1 == v2)
      continue;

    std::chrono::high_resolution_clock::time_point start_actual = std::chrono::high_resolution_clock::now();
    int distance = pll.QueryDistance(v1, v2);
    std::chrono::high_resolution_clock::time_point end_actual = std::chrono::high_resolution_clock::now();

    if (distance == -1 || distance >= cg.nEdges)
    {
      actualSampleFile << (int64_t)v1 << " " << (int64_t)v2 << "\n";
      continue;
    }

    actualDistancesFile << distance << "\n";
    auto duration_actual = duration_cast<std::chrono::nanoseconds>(end_actual - start_actual);

    totalRuntime += (long long int)(duration_actual.count());
    totalRounds++;
    actualRuntimesFile << (long long int)(duration_actual.count()) << "\n";
  }

  printf("completed %d inqueries in %lld nanoseconds \n", totalRounds, totalRuntime);
  actualDistancesFile.close();
  actualRuntimesFile.close();
  actualSampleFile.close();
  return 0;
}