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

// exe/filname [BiBFS/L0] [graphname/graphname_seed_{n}] [start_from] [end_at]
// note: graphname cannot contain any instances of _
int main(int argc, char *argv[])
{
  std::string command = argv[1];

  std::string graph_name = argv[2];
  std::string graph_L0_name = argv[2];

  std::string delimiter = "_";
  std::size_t d_pos = graph_L0_name.find(delimiter);
  if(d_pos != std::string::npos) graph_name = graph_L0_name.substr(0, d_pos);

  std::string str_start = argv[3];
  std::string str_finish = argv[4];

  std::string graph_folder = "graphs/";
  std::string input_folder = "inputs/";
  std::string results_folder = "results/" + graph_name +"/";
  std::string bin_folder = "bin/";

  std::string graph_filename = graph_name + ".edges";
  std::string graph_loc = graph_folder + graph_filename;

  std::string results_prefix = str_start + "_" + str_finish + "_";
  std::string results_suffix = "_" + command;
  std::string results_path = results_folder + results_prefix + graph_L0_name + results_suffix;

  std::string load_bin_prefix = str_start + "_";
  std::string final_save_bin_prefix = str_finish + "_";
  std::string load_bin_path = bin_folder + load_bin_prefix + graph_L0_name + results_suffix;
  std::string final_save_bin_path = bin_folder + final_save_bin_prefix + graph_L0_name + results_suffix;

  std::ofstream graph_file(results_folder + graph_name + results_suffix + ".txt");
  // Graph g;

  


  // if (loadGraph(graph_loc.c_str(), g, 1, IOFormat::escape))
  //   exit(1);
 
  // CGraph cg = makeCSR(g);
  // cg.sortById();
  // printf("Converted to CSR\n");
  // std::cout<<cg.nEdges<<sstd::endl;

  CGraph cg;
  cg.loadGraphFromFile(graph_name);
  int startFrom = stoi(str_start);
  int ROUNDS = stoi(str_finish);

  // initialize varibles for sampling

  //make files
  printf("created variables\n");
  L0Graph L0;

  if(command == "L0"){
    printf("loading L0 and L1\n");
    L0 = L0Graph(cg, graph_L0_name);
    // L0.checkForBadL0();
    // L0.setL0_PLL();
    // printf("setPLL\n");
  }

  std::ofstream nodesTouchedFile;
  std::ofstream queryTouchedFile;
  // std::ofstream nodesTouchedWithoutSamplesFile;

  bool* nodesTouched;
  uint8_t* queryTouched;
  bool* nodesTouchedWithoutSample;

  VertexIdx vcount = cg.nVertices;
  // if(command == "BiBFS"){
  //   nodesTouched = new bool[vcount];
  //   queryTouched = new uint8_t[vcount];
  //   std::fill_n(nodesTouched, vcount, false);
  //   std::fill_n(queryTouched, vcount, 255);
  // } else{
    nodesTouched = new bool[vcount];
    nodesTouchedWithoutSample = new bool[vcount];
    std::fill_n(nodesTouched, vcount, false);
    std::fill_n(nodesTouchedWithoutSample, vcount, false);
  // }

  if(startFrom==0) {
  
    std::fill_n(nodesTouchedWithoutSample, vcount, false);
  } else {
    std::ifstream loadNodesTouchedFile(load_bin_path + "_nodesTouched.bin", std::ios::binary);
    std::ifstream loadNodesTouchedWithoutSamplesFile(load_bin_path + "_nodesTouchedWithoutSamples.bin", std::ios::binary);

    loadNodesTouchedFile.read((char*)nodesTouched, sizeof(bool) * (int64_t) vcount);
    loadNodesTouchedWithoutSamplesFile.read ((char*)nodesTouchedWithoutSample, sizeof(bool) * (int64_t) vcount);

    loadNodesTouchedFile.close();
    loadNodesTouchedWithoutSamplesFile.close();

  }

  std::vector<int> numQueries;
  std::vector<int> queriesAccumulated;
  std::vector<int> queriesAccumulatedWithoutSample;
  std::vector<int> distances;
  std::vector<int> runtimes;
  std::vector<int> L2_to_L0_runtimes;
  std::vector<int> L2_to_L0_dist;
  std::vector<int> L0_hops;

  std::ofstream distancesFile;
  std::ofstream distToReachCoreFile;
  std::ofstream hopsInCoreFile;
  std::ofstream runtimesFile;
  std::ofstream runtimeToReachCoreFile;
  std::ofstream failedSampleFile;
  std::ofstream queriesFile;
  std::ofstream queriesAccumulatedFile;
  std::ofstream queriesAccumulatedWithoutSampleFile;

  distancesFile.open(results_path + "_distances.txt");
  runtimesFile.open(results_path + "_runtimes.txt");
  failedSampleFile.open(results_path + "_failed.txt");
  queriesFile.open(results_path + "_queries.txt");
  
  // if(command=="L0") {
      queriesAccumulatedFile.open(results_path + "_queries_accumulated.txt");
      queriesAccumulatedWithoutSampleFile.open(results_path + "_queries_accumulated_without_sample.txt");
      // queriesListFile.open(results_path + "_queries_list.txt");
  // }
  std::ifstream inputFile;
  inputFile.open("input/" + graph_name + "_input.txt");

  if(command == "L0"){
    hopsInCoreFile.open(results_path + "_hops_in_core.txt");
    distToReachCoreFile.open(results_path + "_distance_to_core.txt");
    runtimeToReachCoreFile.open(results_path + "_runtime_to_core.txt");
  }

  printf("starting rounds\n");
  for(int round = 0; round < ROUNDS; round ++){
      VertexIdx v1;
      VertexIdx v2;
      inputFile >> v1 >> v2;
      // printf("round %i\n", round);
      if(round < startFrom || v1 == v2) continue;

      if (round%1000==0 && round!=startFrom) {
          printf("round %i\n", round);

          for (const auto &e : distances) distancesFile << e << "\n";
          for (const auto &e : runtimes) runtimesFile << e << "\n";
          for (const auto &e : numQueries) queriesFile << e << "\n";
          for (const auto &e : queriesAccumulatedWithoutSample) queriesAccumulatedWithoutSampleFile << e << "\n";
          for (const auto &e : queriesAccumulated) queriesAccumulatedFile << e << "\n";
          if(command == "L0"){
            for (const auto &e : L2_to_L0_runtimes) runtimeToReachCoreFile << e << "\n";
            L2_to_L0_runtimes.clear();
            for (const auto &e : L2_to_L0_dist) distToReachCoreFile << e << "\n";
            L2_to_L0_dist.clear();
            for (const auto &e : L0_hops) hopsInCoreFile << e << "\n";
            L0_hops.clear();
          }

          distances.clear();
          runtimes.clear();
          numQueries.clear();
          queriesAccumulated.clear();
          queriesAccumulatedWithoutSample.clear();

          std::string save_bin_prefix = to_string(round) + "_";
          std::string save_bin_path = bin_folder + save_bin_prefix + graph_L0_name + results_suffix;

          // nodesTouchedFile.open(save_bin_path + "_nodesTouched.bin", std::ios::binary);
          // nodesTouchedWithoutSamplesFile.open(save_bin_path + "_nodesTouchedWithoutSamples.bin", std::ios::binary);

          // nodesTouchedFile.write((char*)nodesTouched, sizeof(bool) * (int64_t) vcount);
          // nodesTouchedWithoutSamplesFile.write((char*)nodesTouchedWithoutSample, sizeof(bool) * (int64_t) vcount);

          // nodesTouchedFile.close();
          // nodesTouchedWithoutSamplesFile.close();

          printf("uploaded through round %i\n", round-1);
      }
      
      if (v1==v2) continue;

      std::chrono::high_resolution_clock::time_point start_clock;
      SamplingInformation sample;
      //int actual = pll_graph.QueryDistance(v1, v2);
      std::chrono::high_resolution_clock::time_point end_clock;
      start_clock = std::chrono::high_resolution_clock::now();
      if(command=="L0")  {
        sample = L0.estimate_distance(v1,v2, false); // for now don't run with PLL
      } else sample = cg.BidirectionalBFS(v1, v2);
      end_clock = std::chrono::high_resolution_clock::now();
      if(sample.estimate_distance==-1) {
        failedSampleFile << (int64_t) v1 << " " << (int64_t) v2 << " " << round <<"\n";
        continue;
      }
      if(command=="L0"){
        L2_to_L0_runtimes.push_back(sample.p1_time_to_L0);
        L2_to_L0_dist.push_back(sample.p1_dist_to_L0);
        L2_to_L0_runtimes.push_back(sample.p2_time_to_L0);
        L2_to_L0_dist.push_back(sample.p2_dist_to_L0);
        L0_hops.push_back(sample.p1_L0_hops);
        L0_hops.push_back(sample.p2_L0_hops);
      }
      distances.push_back((double)(sample.estimate_distance));
      auto duration = duration_cast<std::chrono::microseconds>(end_clock - start_clock);
      runtimes.push_back((double)(duration.count()));

      sample.visited[v1] = false;
      sample.visited[v2] = false;

      numQueries.push_back(sample.countQueries());
      for(VertexIdx i = 0; i< vcount; i++){
          if(sample.visited[i]){
              // if(command=="BiBFS"&& !nodesTouched[i]){
              //     nodesTouched[i] = true;
              //     queryTouched[i] = round - startFrom;
              // } else if (command=="L0"){
                  nodesTouched[i] = true;
                  nodesTouchedWithoutSample[i] = true;
                //queriesListFile << (int64_t) i << " ";
              // }
              
          }
      }
      // if(command=="L0"){
      //   queriesListFile << "\n";
      // }
      nodesTouched[v1] = true;
      nodesTouched[v2] = true;

      queriesAccumulatedWithoutSample.push_back(std::count(nodesTouchedWithoutSample, nodesTouchedWithoutSample + (int64_t) vcount, true));
      queriesAccumulated.push_back(std::count(nodesTouched, nodesTouched + (int64_t) vcount, true));

      delete[] sample.visited;
  }
  printf("ending rounds\n");
  for (const auto &e : distances) distancesFile << e << "\n";
  for (const auto &e : runtimes) runtimesFile << e << "\n";
  for (const auto &e : numQueries) queriesFile << e << "\n";
  for (const auto &e : queriesAccumulatedWithoutSample) queriesAccumulatedWithoutSampleFile << e << "\n";
    for (const auto &e : queriesAccumulated) queriesAccumulatedFile << e << "\n";
  if(command == "L0"){
    for (const auto &e : L2_to_L0_runtimes) runtimeToReachCoreFile << e << "\n";
    runtimeToReachCoreFile.close();
    for (const auto &e : L2_to_L0_dist) distToReachCoreFile << e << "\n";
    distToReachCoreFile.close();
    for (const auto &e : L0_hops) hopsInCoreFile << e << "\n";
    hopsInCoreFile.close();
    //queriesListFile.close();
    queriesAccumulatedFile.close();
    queriesAccumulatedWithoutSampleFile.close();
  }
//   } else{
    
//     nodesTouchedFile.open(final_save_bin_path + "_nodesTouched.bin", std::ios::binary);
//     queryTouchedFile.open(final_save_bin_path + "_queryTouched.bin", std::ios::binary);
//     nodesTouchedFile.write((char*)nodesTouched, sizeof(bool) * (int64_t) vcount);
//     queryTouchedFile.write((char*)queryTouched, sizeof(uint8_t) * (int64_t) vcount);
//     nodesTouchedFile.close();
//     queryTouchedFile.close();
//   }
  //nodesTouchedWithoutSamplesFile.open(final_save_bin_path + "_nodesTouchedWithoutSamples.bin", std::ios::binary);
  //nodesTouchedFile.clear();
  //nodesTouchedWithoutSamplesFile.clear();
  // nodesTouchedFile.seekp(0);
  // nodesTouchedWithoutSamplesFile.seekp(0);
  // nodesTouchedWithoutSamplesFile.write((char*)nodesTouchedWithoutSample, sizeof(bool) * (int64_t) vcount);

  distancesFile.close();
  runtimesFile.close();
  failedSampleFile.close();
  //queriesFile.close();
  // nodesTouchedWithoutSamplesFile.close();
  return 0;
  
}