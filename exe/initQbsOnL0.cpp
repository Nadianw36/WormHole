#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Graph.h"
#include "Escape/L0Graph.h"
#include "Escape/Config.h"
#include "QbS/Preprocess.h"
#include "QbS/Sparsify.h"
#include "QbS/IndexConstruction.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <chrono>
#include <string>
#include <cstring>
#include <fstream>
#include <filesystem>

using namespace Escape;
using namespace QbS;
using namespace std;
namespace fs = std::filesystem;

int k = 2;
// exe/filname [graphname_seed_{n}]
// note: graphname cannot contain any instances of _
int main(int argc, char *argv[])
{
  // grab info from args
  // printf("init\n");
  std::string graph_name = argv[1];
  std::string graph_L0_name = argv[1];

  std::string delimiter = "_";
  std::size_t d_pos = graph_L0_name.find(delimiter);
  if (d_pos != std::string::npos)
    graph_name = graph_L0_name.substr(0, d_pos);
  // printf("check setup\n");

  checkSetupFor(graph_name);

  std::string L0_qbs_results_folder = RESULTS_FOLDER + graph_name + "/" + L0_FOLDER + QBS_FOLDER;
  std::ofstream init_results_file(L0_qbs_results_folder + "init" + ".txt");

  CGraph cg;
  cg.loadGraphFromFile(graph_name);

  // initialize varibles for sampling

  // make files
  // printf("created variables\n");
  L0Graph L0;
  // printf("loading L0 and L1\n");
  L0 = L0Graph(cg, graph_L0_name);
  // printf("generating edgelist\n");
  vector<std::pair<VertexIdx, VertexIdx>> L0_edgelist = L0.L0_edgelist();

  // printf("preprocess L0 with QBS\n");
  cleanData(L0_edgelist, graph_name, L0_qbs_results_folder);

  // printf("sparsify L0 with QBS\n");
  sparsify(20, graph_name, L0_qbs_results_folder);

  // printf("index construction on L0 with QBS\n");
  long constructionTime = indexConstructionEfficient(k, graph_name, L0_qbs_results_folder);

  return 0;
}