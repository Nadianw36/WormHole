#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Graph.h"
#include "Escape/L0Graph.h"
#include "Escape/Config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cctype>

std::vector<std::string> splitString(const std::string &str)
{
  std::istringstream iss(str);
  std::vector<std::string> tokens;

  std::string token;
  while (iss >> token)
  {
    tokens.push_back(token);
  }

  return tokens;
}

using namespace Escape;

// sanitize [graphname] [input_extension]
int main(int argc, char *argv[])
{
  std::string graph_name = argv[1];
  std::string graph_ex = argv[2];
  std::string graph_filename = graph_name + "-new.edges";

  std::string out_graph_loc = GRAPH_FOLDER + graph_filename;
  std::string in_graph_loc = GRAPH_FOLDER + graph_name + graph_ex;
  Graph g;
  if (loadGraph(in_graph_loc.c_str(), g, 1, IOFormat::escape))
    exit(1);
  // printf("Loaded graph\n");
  CGraph cg = makeCSR(g);
  cg.sortById();
  // printf("Converted to CSR\n");
  cg.createSanitizedFile(out_graph_loc);
  //   VertexIdx ind = 0;

  //   std::ifstream inputFile(in_graph_loc);
  //   std::ofstream outputFile(out_graph_loc);
  //   std::unordered_map<std::string, VertexIdx> mapping;
  //   std::unordered_set<std::pair<VertexIdx, VertexIdx>> edges;

  //   std::string eachline;
  //     while (std::getline(inputFile, eachline)) {
  //         // Process the line here
  //         if (eachline[0] == '#') {
  //                 continue; // Ignore lines starting with a hash
  //         }
  //         std::istringstream iss(eachline);
  //         VertexIdx node1;
  //         VertexIdx node2;
  //         inputFile >> node1 >> node2;
  //         if (node1 == node2) {
  //             continue; // Remove self-loop
  //         }
  //         // if (mapping.find(tokens[0]) != mapping.end()) {
  //         //     node1 = mapping[tokens[0]];
  //         // } else {
  //         //     mapping[tokens[0]] = ind;
  //         //     node1 = ind;
  //         //     ind++;
  //         // }
  //         // if (mapping.find(tokens[1]) != mapping.end()) {
  //         //     node2 = mapping[tokens[1]];
  //         // } else {
  //         //     mapping[tokens[1]] = ind;
  //         //     node2 = ind;
  //         //     ind++;
  //         // }
  //         // std::cout << node1 << " " << node2  << std::endl;
  //         std::pair<VertexIdx, VertexIdx> edge(node1, node2);
  //         std::pair<VertexIdx, VertexIdx> reverse_edge(node2, node1);
  //         if(edges.find(reverse_edge) == edges.end()){
  //             edges.insert(edge);
  //         }
  //     }
  //     EdgeIdx m = edges.size();
  //     std::cout << ind << " vertices and " << m << " edges" << std::endl;
  //     outputFile << ind << " " << m << "\n";
  //     for (const auto& edge : edges) {
  //         outputFile  << edge.first << " " << edge.second << "\n";
  //     }

  //     inputFile.close();
  //     outputFile.close();

  return 0;
}