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

// exe/filname [BiBFS/L0] [graphname/graphname_seed_{n}] [start_from] [end_at]
// note: graphname cannot contain any instances of _
int main(int argc, char *argv[])
{
  std::string command = argv[1];

  std::string graph_name = argv[2];
  std::string graph_L0_name = argv[2];

  // std::string delimiter = "_";
  // std::size_t d_pos = graph_L0_name.find(delimiter);
  // if (d_pos != std::string::npos)
  //   graph_name = graph_L0_name.substr(0, d_pos);

  checkSetupFor(graph_name);

  std::string str_start = argv[3];
  std::string str_finish = argv[4];

  std::string results_graph_sub_folder = getSubfolderName(command);
  std::string results_folder = RESULTS_FOLDER + graph_name + "/" + results_graph_sub_folder;
  std::string results_prefix = str_start + "_" + str_finish + "_";
  std::string results_suffix = "_" + command;

  std::string results_path = results_folder + results_prefix + graph_L0_name + results_suffix;

  CGraph cg;
  cg.loadGraphFromFile(graph_name);
  BiBFSGraph bg = BiBFSGraph(cg);
  int startFrom = stoi(str_start);
  int ROUNDS = stoi(str_finish);

  // make files
  L0Graph L0 = L0Graph(cg, graph_L0_name);
  L0.checkForBadL0();

  VertexIdx vcount = cg.nVertices;
  bool *nodesTouched = new bool[vcount];
  bool *nodesTouchedOutsideCore = new bool[vcount];
  std::fill_n(nodesTouched, vcount, false);
  std::fill_n(nodesTouchedOutsideCore, vcount, false);

  std::ofstream distancesFile(results_path + "_distances.txt");
  std::ofstream pathsFile(results_path + "_paths.txt");
  std::ofstream runtimesFile(results_path + "_runtimes_nano.txt");
  std::ofstream failedSampleFile(results_path + "_failed.txt");
  std::ofstream queriesFile(results_path + "_queries.txt");
  std::ofstream queriesAccumulatedFile(results_path + "_queries-accumulated.txt");

  std::ofstream leftPathsFile;
  std::ofstream rightPathsFile;
  std::ofstream queriesOutsideCoreFile;
  std::ofstream runtimeToReachCoreFile;
  std::ofstream distToReachCoreFile;
  std::ofstream queriesAccumulatedOutsideCoreFile;

  std::ofstream coreQueriesRandomL0File;
  std::ofstream coreQueriesRandomL0RandomL1File;
  std::ofstream coreQueriesHighestDegL0File;
  std::ofstream coreQueriesAllPairsL0;

  if (command == "L0-BiBFS")
  {
    distToReachCoreFile.open(results_path + "_dist-to-core.txt");
    runtimeToReachCoreFile.open(results_path + "_runtime-to-core_nano.txt");
    queriesOutsideCoreFile.open(results_path + "_queries-outside-core.txt");
    queriesAccumulatedOutsideCoreFile.open(results_path + "_queries-accumulated-outside-core.txt");

    leftPathsFile.open(results_path + "_paths_v1-to-core.txt");
    rightPathsFile.open(results_path + "_paths_core-to-v2.txt");

    coreQueriesRandomL0File.open(results_path + "_core-queries_random-L0.txt");
    coreQueriesRandomL0RandomL1File.open(results_path + "_core-queries_random-L0-random-L1.txt");
    coreQueriesHighestDegL0File.open(results_path + "_core-queries_highest-deg-L0.txt");
    coreQueriesAllPairsL0.open(results_path + "_core-queries_all-pairs-L0.txt");
  }

  std::ifstream inputFile(INPUT_FOLDER + graph_name + "_input.txt");
  printf("starting round\n");
  VertexIdx v1, v2;
  long long int totalRuntime = 0;
  long long int totalTimeSpentOutsideCore = 0;
  int inquryCount = 0;
  for (int round = 0; round < ROUNDS; round++)
  {
    inputFile >> v1 >> v2;
    // std::cout << v1 << " " << v2 << std::endl;
    if (round < startFrom || v1 == v2)
      continue;

    std::chrono::high_resolution_clock::time_point start_clock;
    std::chrono::high_resolution_clock::time_point end_clock;
    int distance;
    bool *visited;

    if (command == "L0-BiBFS")
    {
      start_clock = std::chrono::high_resolution_clock::now();
      distance = L0.estimate_distance(v1, v2, false); // for now don't run with PLL
      end_clock = std::chrono::high_resolution_clock::now();
      visited = L0.getVisited();
    }
    else
    {
      start_clock = std::chrono::high_resolution_clock::now();
      distance = bg.BidirectionalBFS(v1, v2);
      end_clock = std::chrono::high_resolution_clock::now();
      visited = bg.getVisited();
    }

    if (distance == -1)
    {
      distancesFile << -1 << std::endl;
      failedSampleFile << (int64_t)v1 << " " << (int64_t)v2 << " " << round << "\n";
      continue;
    }

    visited[v1] = true;
    visited[v2] = true;

    VertexIdx countVerticesOutsideCore = 0;
    for (VertexIdx i = 0; i < vcount; i++)
    {
      if (visited[i])
      {
        nodesTouched[i] = true;
        if (command == "L0-BiBFS" && !L0.L0[i])
        {
          countVerticesOutsideCore++;
          nodesTouchedOutsideCore[i] = true;
        }
      }
    }

    inquryCount++;
    distancesFile << distance << "\n";
    auto duration = duration_cast<std::chrono::nanoseconds>(end_clock - start_clock);
    totalRuntime += (long long int)(duration.count());
    runtimesFile << (long long int)(duration.count()) << "\n";
    queriesFile << countQueries(visited, bg.nVertices) << "\n";
    queriesAccumulatedFile << countQueries(nodesTouched, vcount) << "\n";

    if (command == "L0-BiBFS")
    {
      auto duration = duration_cast<std::chrono::nanoseconds>(L0.outside_core_end_clock - L0.outside_core_start_clock);
      totalTimeSpentOutsideCore += (long long int)(duration.count());

      runtimeToReachCoreFile << (long long int)(duration.count()) << "\n";
      distToReachCoreFile << L0.acc_distance << "\n";
      queriesOutsideCoreFile << countVerticesOutsideCore << "\n";
      queriesAccumulatedOutsideCoreFile << countQueries(nodesTouchedOutsideCore, vcount) << "\n";

      vector<VertexIdx> path = L0.recoverPath(v1, v2);
      bool reachedCore = false;
      for (const auto &p : path)
      {
        pathsFile << p << " ";
        if (L0.acc_distance != distance)
        {
          if (!L0.L0[p])
            if (!reachedCore)
              leftPathsFile << p << " ";
            else
              rightPathsFile << p << " ";
          else
            reachedCore = true;
        }
      }
      pathsFile << "\n";
      leftPathsFile << "\n";
      rightPathsFile << "\n";

      L0.sampleL0Vertices(v1, v2);
      coreQueriesRandomL0File << L0.randomL0V1 << " " << L0.randomL0V2 << "\n";
      coreQueriesRandomL0RandomL1File << L0.randomL0FromRandomL1V1 << " " << L0.randomL0FromRandomL1V2 << "\n";
      coreQueriesHighestDegL0File << L0.highestDegL0V1 << " " << L0.highestDegL0V2 << "\n";

      if (L0.acc_distance == distance)
        coreQueriesAllPairsL0 << -1 << " " << -1 << "\n";
      else
      {
        for (VertexIdx i = L0.Q_START; i < L0.level1Record[L0.LEVEL_START]; i++)
        {
          VertexIdx v1 = L0.queue1[i];
          for (VertexIdx j = L0.Q_START; j < L0.level2Record[L0.LEVEL_START]; j++)
          {
            VertexIdx v2 = L0.queue2[j];
            coreQueriesAllPairsL0 << v1 << " " << v2 << "\n";
          }
        }
      }
      coreQueriesAllPairsL0 << -2 << " " << -2 << "\n";
    }
    else
    {
      vector<VertexIdx> path = bg.recoverPath(v1, v2);
      for (const auto &p : path)
        pathsFile << p << " ";
      pathsFile << "\n";
    }
  }
  if (command == "L0-BiBFS")
  {
    printf("Wormhole-E ran %d inqueries in %lld nanoseconds \n", inquryCount, totalRuntime);
    printf("Outside the core, Wormhole-E ran %d inqueries in %lld nanoseconds \n", inquryCount, totalTimeSpentOutsideCore);
    runtimeToReachCoreFile.close();
    distToReachCoreFile.close();
    queriesOutsideCoreFile.close();
    queriesAccumulatedOutsideCoreFile.close();

    leftPathsFile.close();
    rightPathsFile.close();

    coreQueriesAllPairsL0.close();
    coreQueriesHighestDegL0File.close();
    coreQueriesRandomL0File.close();
    coreQueriesRandomL0RandomL1File.close();
  }
  else
    printf("BiBFS ran %d inqueries in %lld nanoseconds \n", inquryCount, totalRuntime);

  distancesFile.close();
  pathsFile.close();
  runtimesFile.close();
  failedSampleFile.close();
  queriesAccumulatedFile.close();

  return 0;
}