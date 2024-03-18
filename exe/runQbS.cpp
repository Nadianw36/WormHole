#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Graph.h"
#include "Escape/L0Graph.h"
#include "Escape/Config.h"
#include "QbS/Preprocess.h"
#include "QbS/Sparsify.h"
#include "QbS/IndexConstruction.h"
#include "QbS/QueryQbS.h"

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

int k = 20;
// exe/filname [graphname_seed_{n}]
// note: graphname cannot contain any instances of _
int main(int argc, char *argv[])
{
  // grab info from args
  std::string graph_name = argv[1];

  std::string graph_loc = GRAPH_FOLDER + graph_name + ".edges";
  checkSetupFor(graph_name);

  std::string qbs_results_folder = RESULTS_FOLDER + graph_name + "/" + QBS_FOLDER;
  std::ofstream init_results_file(qbs_results_folder + "init" + ".txt");

  CGraph cg;
  cg.loadGraphFromFile(graph_name);

  std::vector<std::pair<VertexIdx, VertexIdx>> edgelist;
  // printf("neighbors of a in cg \n");
  for (VertexIdx i = 0; i < cg.nVertices; i++)
  {
    for (EdgeIdx j = cg.offsets[i]; j < cg.offsets[i + 1]; j++)
    {
      VertexIdx nbor = cg.nbors[j];
      std::pair<VertexIdx, VertexIdx> edge(i, nbor);
      edgelist.push_back(edge);
      // if(i == 16807) // printf("%d ", j);
    }
  }
  // printf("\n");

  // printf("preprocess with QbS\n");
  cleanData(edgelist, graph_name, qbs_results_folder);

  // printf("sparsify with QbS\n");
  sparsify(20, graph_name, qbs_results_folder);

  // printf("index construction on L0 with QBS hiii\n");
  long constructionTime = indexConstructionEfficient(k, graph_name, qbs_results_folder);

  // printf("running experiments\n");
  int diameter = 17;
  spQuery test = spQuery(k, diameter, graph_name, qbs_results_folder);
  test.generateTest(10000, INPUT_FOLDER + graph_name + "_input.txt");
  test.seriesTest(10000);

  return 0;
}