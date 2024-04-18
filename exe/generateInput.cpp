#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Graph.h"
#include "Escape/L0Graph.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <chrono>
#include <string>
#include <cstring>

using namespace Escape;

// generateInput [graph_name]
int main(int argc, char *argv[])
{
  std::string graph_name = argv[1];
  std::ofstream inputFile(INPUT_FOLDER + graph_name + "_input.txt");

  CGraph cg;
  cg.loadGraphFromFile(graph_name);

  int ROUNDS = 10000;
  int count = 0;
  VertexIdx nVertices = cg.nVertices;

  while (count < ROUNDS)
  {
    auto v1 = rand() % (int64_t)nVertices;
    auto v2 = rand() % (int64_t)nVertices;
    if (v1 == v2)
      continue;
    // SamplingInformation distance = cg.BidirectionalBFS(v1, v2);
    // if (distance.estimate_distance < 5)
    //   continue;
    // if(count%100==0) // printf("%d\n", count);
    count += 1;
    inputFile << (int64_t)v1 << " " << (int64_t)v2 << "\n";
    inputFile.flush();
  }
  inputFile.close();
  return 0;
}