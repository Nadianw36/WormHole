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
#include <fstream>


using namespace Escape;
using namespace std::chrono;

// exe/createL0SampleFile [graph_name] [L0_size (ex "10.0")]
int main(int argc, char *argv[])
{
  std::string  graph_name = argv[1];
  std::string L0_size = argv[2];
  std::string graph_filename = graph_name + ".edges";
  std::string graph_folder = "graphs/";
  std::string graph_loc = graph_folder + graph_filename;
  std::string graph_L0_name = graph_name + "_seed_" + L0_size;
  
  std::string results_folder = "results/";
  std::string results_file = results_folder + "10000_" + graph_L0_name + "_L0_distances.txt";
  CGraph cg;
  cg.loadGraphFromFile(graph_name);

  // initialize varibles for sampling

  //make files
  printf("created variables\n");
  L0Graph L0 = L0Graph(cg, graph_L0_name);;
  VertexIdx nVerticesL0 = L0.nVerticesL0;

  std::string folder = "input/";
  std::ofstream inputFile(folder + graph_L0_name + "_L0_input.txt");
  std::ofstream outputFile(results_file);

  int ROUNDS = 10000;
  int count = 0;

  while(count < ROUNDS){
      auto v1 = rand() % (int64_t) nVerticesL0;
      auto v2 = rand() % (int64_t) nVerticesL0;
      //   printf("vs %ld %ld\n", v1, v2);
      if (v1==v2) continue;
      count+=1;
      VertexIdx L0_v1 = 0;
      VertexIdx L0_v2 = 0;
      while(v1 >= 0 || v2 >= 0){
        // printf("status1 v1 %ld v2 %ld \n", v1, v2);
        // printf("    status2 L0v1 %ld L0v2 %ld \n", L0_v1, L0_v2);
        // printf("    status3 v1_inL0 %d v2_inL0 %d \n", L0.L0[L0_v1], L0.L0[L0_v2]);
        if(v1 >= 0 && L0.L0[L0_v1]) v1--;
        if(v2 >= 0 && L0.L0[L0_v2]) v2--;
        if(v1 >= 0) L0_v1++;
        if(v2 >= 0) L0_v2++;
        // printf("    status4 v1 %ld v2 %ld \n", v1, v2);
        // printf("    status5 L0v1 %ld L0v2 %ld \n", L0_v1, L0_v2);
        // printf("    status6 v1_inL0 %d v2_inL0 %d \n", L0.L0[L0_v1], L0.L0[L0_v2]);
      }
      //   printf("L0 vs %ld %ld\n", L0_v1, L0_v2);
      if(!L0.L0[L0_v1] || !L0.L0[L0_v2]) throw std::logic_error("something wrong");
      inputFile << (int64_t) L0_v1 << " " << (int64_t) L0_v2 << "\n";
      set<VertexIdx> p1({L0_v1});
      set<VertexIdx> p2({L0_v2});
      L0_to_L0 L0_info = L0.distance_L0_to_L0(p1, p2);
      outputFile << L0_info.dist << "\n";
  }
  inputFile.close();
  outputFile.close();

  
  return 0;
  
}