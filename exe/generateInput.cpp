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

//generateInput [graph_name] [input_ind]
int main(int argc, char *argv[])
{
  std::string  graph_name = argv[1];
  std::string input_ind = argv[2];
  std::string graph_filename = graph_name + ".edges";
  std::string graph_folder = "graphs/";
  std::string graph_loc = graph_folder + graph_filename;
  // Graph g;
  // if (loadGraph(graph_loc.c_str(), g, 1, IOFormat::escape))
  // exit(1);
  // printf("Loaded graph\n");
  // CGraph cg = makeCSR(g);
  // cg.sortById();
  // printf("Converted to CSR\n");
  CGraph cg;
  cg.loadGraphFromFile(graph_name);

  int ROUNDS = 100;
  int count = 0;
  VertexIdx nVertices = cg.nVertices;
  // VertexIdx nVertices = 41652230;
  std::string folder = "input/";
  std::ofstream inputFile(folder + graph_name + "_" + input_ind +"_input.txt");
  while(count < ROUNDS){
      auto v1 = rand() % (int64_t) nVertices;
      auto v2 = rand() % (int64_t) nVertices;
      if (v1==v2) continue;
      SamplingInformation distance = cg.BidirectionalBFS(v1, v2);
      if(distance.estimate_distance<5) continue;
      // if(count%100==0) printf("%d\n", count);
      count+=1;
      inputFile << (int64_t) v1 << " " << (int64_t) v2 << "\n";
      inputFile.flush();
  }
  inputFile.close();

  return 0;
  
}