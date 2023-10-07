#ifndef ESCAPE_L0GRAPH_H_
#define ESCAPE_L0GRAPH_H_

#include <cstdlib>
#include <stdexcept>
#include <cstdio>
#include <numeric>
#include <algorithm>
#include <set>
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <chrono>
#include "Escape/Graph.h"
#include "Escape/PrunedLandmark.h"
#include "Escape/IndexedBinaryHeap.h"
#include <queue>

namespace Escape
{

std::string FOLDER = "bin/";

// data structure for storing sampling information from Ln to L2 path
struct Ln_to_L2{
  VertexIdx L2_pt;
  int distance;
};
struct Ln_to_L0{
  int distance;
  bool early_connection;
};
struct L2_to_L2_info{
  VertexIdx p1_L2_pt;
  VertexIdx p2_L2_pt;
  int acc_dist; // set to -1 if not needed
  SamplingInformation estimate;
};
struct L0_to_L0{
  int dist;
  int p1_hops;
  int p2_hops;
};

class L0Graph{
  public:
  CGraph graph;
  bool *L0;  // bitmap of boolean flags
  bool *L1; // bitmap of boolean flags
  int *apsp_L0; // deprecated
  VertexIdx nVerticesL0;
  PrunedLandmarkLabeling<> pll;
  std::ofstream output;
  L0Graph(){}
  /**
   * @brief Construct a new L0Graph object from binary files containing L0 and L1 boolean arrays. 
   * Files should be named and have path bin/{graphName}_L0.bin and bin/{graphName}_L1.bin
   * 
   * @param g 
   * @param graphName 
   */
  L0Graph(CGraph g, std::string graphName){
    graph = g;
    auto vcount = g.nVertices;
    L0 = new bool[vcount];
    L1 = new bool[vcount];

    std::ifstream L0File(FOLDER + graphName + "_L0.bin", std::ios::in | std::ios::binary);
    std::ifstream L1File(FOLDER + graphName + "_L1.bin", std::ios::in | std::ios::binary);
    L0File.read ((char*)L0, sizeof (bool)* (int64_t) vcount);
    L1File.read ((char*)L1, sizeof (bool)* (int64_t) vcount);
    L0File.close();
    L1File.close();

    nVerticesL0 = std::count(L0, L0 + vcount, true);
    printf("L0 size: %ld \n", nVerticesL0);
  }
  /**
   * @brief Construct a new L0Graph object via algorithm with an L0 size of sizeL0
   * 
   * @param g 
   * @param sizeL0 - number of vertices in L0
   * @param optimal - if True, then use high degree vertice as seed vertice, else use random vertice
   */
  L0Graph (CGraph g, VertexIdx sizeL0, bool optimal) {
      output.open("output/output.txt");
      auto vcount = g.nVertices;
      nVerticesL0 = sizeL0;
      L0 = new bool[vcount];
      L1 = new bool[vcount];
      std::fill_n(L0, vcount, false);
      std::fill_n(L1, vcount, false);
      IndexedBinaryHeap L1_map = IndexedBinaryHeap(&g); // maps nodes currently in L1 -> L0 degree

      auto v = rand() % (int16_t) vcount; // picks seed node
      if (optimal){
        VertexIdx *p = new VertexIdx[vcount];
        VertexIdx sampleThresh = (VertexIdx) (0.01*vcount);
        for (VertexIdx v = 0; v < vcount; v++){
          p[v] = v;
        }
        std::sort(p, p+vcount, [&g] (VertexIdx v1, VertexIdx v2){return g.degree(v1) > g.degree(v2);});
        VertexIdx sample = rand() % (int16_t) sampleThresh;
        v = p[sample];
        delete p;
      }

      L0[v] = true;
      for (EdgeIdx i = g.offsets[v]; i < g.offsets[v+1]; i++){
        VertexIdx nbor = g.nbors[i];
        L1_map.increment(nbor);
      }
      for (int i = 1; i < sizeL0; i++){
        auto max_v = L1_map.popMax();
        L0[max_v] = true; 
        for (EdgeIdx j = g.offsets[max_v]; j < g.offsets[max_v+1]; j++){
          VertexIdx nbor = g.nbors[j];
          if(!L0[nbor]){ 
            L1_map.increment(nbor);
          }
          }
      }
      for (auto it = L1_map.keyToPosition.cbegin(); it !=L1_map.keyToPosition.cend(); ++it )
      {
        L1[it->first] = true;
      }
      graph = g;
  }
  /**
   * @brief Assert that:
   *            1) for every vertex in L0, all of its neighbors are in L1
   *            2) for every vertex in L1, at least one neighbor is in L0
   */
  void checkForBadL0() const {
    for(VertexIdx v = 0; v < graph.nVertices; v++){
      if(L0[v]){
        for(EdgeIdx e = graph.offsets[v]; e < graph.offsets[v+1]; e++){
          VertexIdx nbor = graph.nbors[e];
          if(!(L0[nbor] || L1[nbor])) throw std::invalid_argument( "all neibors of L0 should be in L0 or L1" );
        }

      } else if(L1[v]){
        int countL0 = 0;
        for(EdgeIdx e = graph.offsets[v]; e < graph.offsets[v+1]; e++){
          VertexIdx nbor = graph.nbors[e];
          if(L0[nbor]) countL0++;
        }
        if(!(countL0 > 0)) throw std::invalid_argument("L1 vertices are all connected to L0");
      }
    }
  }
  /**
   * @brief Create L0 from the {sizeL0} highest degree vertices in graph {cg}. 
   * If resulting L0 is disconnected, then make L0 the largest component 
   * (meaning that the resulting L0 might have less vertices than {sizeL0})
   * 
   * @param cg 
   * @param sizeL0 - number of vertices in L0
   * @param graphName 
   */
  void generateL0FromPriority(CGraph cg, VertexIdx sizeL0, std::string graphName){
    VertexIdx vcount = cg.nVertices;
    nVerticesL0 = sizeL0;
    graph = cg;
    std::priority_queue<std::pair<EdgeIdx, VertexIdx>> p;
    L0 = new bool[vcount];
    L1 = new bool[vcount];
    std::fill_n(L0, vcount, false);
    std::fill_n(L1, vcount, false);
    for (VertexIdx v = 0; v < vcount; v++){
      EdgeIdx degrees = cg.degree(v);
      p.push(std::make_pair(degrees, v));
    }
    printf("generated pq\n");
    for(VertexIdx i = 0; i < sizeL0; i++){
      VertexIdx vl0 = p.top().second;
      p.pop();
      L0[vl0] = true;
      L1[vl0] = false;
      for(EdgeIdx e = cg.offsets[vl0]; e < cg.offsets[vl0+1]; e++){
        VertexIdx nbor = cg.nbors[e];
        if(!L0[nbor]){
          L1[nbor] = true;
        }
      }
    }
  }
  /**
   * @brief Prints the number of connected components, the size of each component, and
   * the size of the largest component. Also, on default, prunes everything but the largest component. 
   * 
   * @param prune - if True, then prune everything but the largest component; else don't
   * @return VertexIdx - number of vertices in L0
   */
  VertexIdx connectComponentsL0(bool prune = true){
    int count = 0;
    VertexIdx num_nodes  = 0;
    VertexIdx largest_component = 0;
    VertexIdx* component_size = new VertexIdx[graph.nVertices];
    std::fill_n(component_size, graph.nVertices, 0);

    bool *visited = new bool[graph.nVertices];
    std::fill_n(visited, graph.nVertices, false);
    // find components
    for(VertexIdx i = 0; i < graph.nVertices; i++){
      if(!L0[i] || visited[i]) continue;
      count++;
      VertexIdx seed  = i;
      std::queue<VertexIdx> q;
      q.push(seed);
      visited[seed] = true;
      int size = 0;
      while(!q.empty()){
        num_nodes++;
        VertexIdx v = q.front();
        size++;
        q.pop();
        for(EdgeIdx e = graph.offsets[v]; e < graph.offsets[v+1]; e++){
          VertexIdx nbor = graph.nbors[e];
          if(L0[nbor] && !visited[nbor]){
            visited[nbor] = true;
            q.push(nbor);
          }
        }
      }
      component_size[i] = size;
      if (size>largest_component) largest_component = size;
      printf("size of connected component: %d\n", size);
      //printf("L0_size: %d\n", num_nodes);
     }
     printf("number of components: %d\n", count);
     //printf("L0_size: %d\n", num_nodes);
     printf("largest component: %d\n", largest_component);
     // prune
     for(VertexIdx i = 0; i < graph.nVertices; i++){
      if (!L0[i] || component_size[i] == largest_component || component_size[i] == 0) continue;
      count--;
      VertexIdx seed = i;
      std::queue<VertexIdx> q;
      q.push(seed);
      L0[seed] = false; 
      while(!q.empty()){
        VertexIdx v = q.front();
        num_nodes --;
        q.pop();
        for(EdgeIdx e = graph.offsets[v]; e < graph.offsets[v+1]; e++){
          VertexIdx nbor = graph.nbors[e];
          if(L0[nbor]){
            q.push(nbor);
            L0[nbor] = false; 
          }else if(L1[nbor]) L1[nbor] = false;
        }
      }
     }
     printf("number of components: %d\n", count);
    return num_nodes;
  }
  std::vector<std::pair<int, int>> L0_edgelist() {
    std::vector<std::pair<int, int>> edgelist;
    for (VertexIdx i = 0; i<graph.nVertices; i++){
      if(!L0[i]) continue;
      for (EdgeIdx j = graph.offsets[i]; j < graph.offsets[i+1]; j++){
        VertexIdx nbor = graph.nbors[j];
        if(L0[nbor]) {
          std::pair<int, int> edge((int64_t) i,(int64_t) nbor);
          edgelist.push_back(edge);
        }
      }
    }
    return edgelist;
  }
  int L0_APSP_Index(VertexIdx v) const {
    int count = 0;
    for(VertexIdx i = 0; i<graph.nVertices; i++){
      if (v==i) return count;
      if (L0[i]) count++;
    }
  }
  /**
   * @brief Set the L0 PLL object
   * Constructs a Pruned Landmark Labelling Graph in L0 to speed up distance computing 
   */
  void setL0_PLL() {
    std::vector<std::pair<int, int>> edgelist;
    for (VertexIdx i = 0; i<graph.nVertices; i++){
        if(L0[i]) {
          for (EdgeIdx j = graph.offsets[i]; j < graph.offsets[i+1]; j++){
            VertexIdx nbor = graph.nbors[j];
            if(L0[nbor]) {
              std::pair<int, int> edge((int64_t) i,(int64_t) nbor);
              edgelist.push_back(edge);
            }
          }
        }
        
      }
    pll.ConstructIndex(edgelist);
  }
  
  /**
   * @brief Save L0 and L1 boolean arrays as binary files named bin/{graphName}_L0.bin and bin/{graphName}_L1.bin.
   * 
   * @param graphName 
   */
  void writeGraphToFile(std::string graphName){
    std::ofstream L0File(FOLDER + graphName + "_L0.bin", std::ios::out | std::ios::binary);
    std::ofstream L1File(FOLDER + graphName + "_L1.bin", std::ios::out | std::ios::binary);
    L0File.write ((char*)L0, sizeof(bool) * (int64_t) graph.nVertices);
    L1File.write ((char*)L1, sizeof(bool) * (int64_t) graph.nVertices);
    L0File.close();
    L1File.close();
  }
  /**
   * @brief Print L0 -> L0, L0 -> L1, L1 -> L0, and L1 -> L1 degrees, as well as the total degree file
   * in results/ folder. 
   * 
   * @param graph_name 
   */
  void print_size(std::ostream& graph_file, std::string graph_name) const{
    std::string results_folder = "results/";
    std::ofstream degreeFile(results_folder+ graph_name + "_degrees.txt");
    std::ofstream L0DegreeFile(results_folder + graph_name + "_degreesL0.txt");
    std::ofstream L1DegreeFile(results_folder + graph_name + "_degreesL1.txt");
    std::ofstream L0L1DegreeFile(results_folder + graph_name + "_degreesL0L1.txt");
    std::ofstream L1L0DegreeFile(results_folder + graph_name + "_degreesL1L0.txt");
    std::ofstream L2L1DegreeFile(results_folder + graph_name + "_degreesL2L1.txt");

    std::vector<int> degrees;
    std::vector<int> L0_degrees;
    std::vector<int> L1_degrees;
    std::vector<int> L1_L0_degrees;
    std::vector<int> L0_L1_degrees;
    std::vector<int> L2_L1_degrees;
    int mistakes = 0;
    for (VertexIdx i = 0; i<graph.nVertices; i++){
      if(L0[i]&&L1[i]) mistakes++;
      degrees.push_back(graph.offsets[i+1] - graph.offsets[i]);

      int L0_degree = 0;
      int L0_L1_degree = 0;
      int L1_degree = 0;
      int L1_L0_degree = 0;
      int L2_L1_degree = 0;

      for (EdgeIdx j = graph.offsets[i]; j < graph.offsets[i+1]; j++){
          VertexIdx nbor = graph.nbors[j];
          if(L0[i] && L0[nbor]) {L0_degree++;} // L0 -> L0
          else if (L0[i] && L1[nbor]) { L0_L1_degree++; } // L0 -> L1
          else if(L1[i] && L1[nbor]) {L1_degree++;} // L1 -> L1
          else if ( L1[i] && L0[nbor]) { L1_L0_degree++; } // L1 -> L0
          else if (!L0[i] && !L1[i] && L1[nbor]){L2_L1_degree++;} // L2 -> L1
        }

      if(L0[i]) {
        L0_degrees.push_back(L0_degree);
        L0_L1_degrees.push_back(L0_L1_degree);
      } else if(L1[i]) {
        L1_degrees.push_back(L1_degree);
        L1_L0_degrees.push_back(L1_L0_degree);
      } else {
        L2_L1_degrees.push_back(L2_L1_degree);
      }
    }
    for (const auto &e : degrees) degreeFile << e << "\n";
    for (const auto &e : L0_degrees) L0DegreeFile << e << "\n";
    for (const auto &e : L1_degrees) L1DegreeFile << e << "\n";
    for (const auto &e : L0_L1_degrees) L0L1DegreeFile << e << "\n";
    for (const auto &e : L1_L0_degrees) L1L0DegreeFile << e << "\n";
    for (const auto &e : L2_L1_degrees) L2L1DegreeFile << e << "\n";

    graph_file << graph_name << " size L0: " << L0_degrees.size() << "\n";
    graph_file << graph_name << " size L1: " << L1_degrees.size() << "\n";
  };
  /**
   * @brief Computer actual distance between p1 and p2
   * 
   * @param p1 - source Vertex
   * @param p2 - destination Vertex
   * @return SamplingInformation which contains true distance as well as visited vertices
   */
  SamplingInformation actualDistance(VertexIdx p1, VertexIdx p2) const {
    return graph.BidirectionalBFS(p1, p2);
  }
  /** MAIN ALGORITHM 
   * @brief Run L0 algorithm to estimate distance between p1 and p2
   * 
   * @param p1 
   * @param p2 
   * @param optimized - If true, use PLL in L0; else use BiBFS
   * @return SamplingInformation which contains estimated distance as well as visited vertices (excluding visited L0 vertices)
   */
  SamplingInformation estimate_distance(VertexIdx p1, VertexIdx p2, bool optimized) { 
    bool* visited;
    SamplingInformation result = SamplingInformation(visited, graph.nVertices);
    int distance = 0;
    // printf("L0: %d, %d\n", L0[p1], L0[p2]);
    // printf("L1: %d, %d\n", L1[p1], L1[p2]);
    if(p1==p2) return SamplingInformation(visited, 0, graph.nVertices);
    
    VertexIdx L2_p1 = -1;
    VertexIdx L2_p2 = -1;
    if(!L0[p1] && !L1[p1] && !L0[p2] && !L1[p2]){ // BiBFS within L2 neighborhood(s)
      chrono::high_resolution_clock::time_point start_clock = chrono::high_resolution_clock::now();
      L2_to_L2_info p2_distance = L2_to_L2(p1, p2);
      chrono::high_resolution_clock::time_point end_clock = chrono::high_resolution_clock::now();
      // printf("L2->L2 est_dist: %d\n", p2_distance.estimate.estimate_distance);
      // printf("L2->L2 acc_dist: %d\n", p2_distance.acc_dist);

      if (p2_distance.estimate.estimate_distance>0 || p2_distance.acc_dist < 0) return p2_distance.estimate; // p2_distance.acc_dist < 0 means that at least one of the neighborhoods is disconnected

      L2_p1 = p2_distance.p1_L2_pt; // while traversing L>=2 neighborhoods, we collect the first L2 point we hit so we don't need to traverse a second time
      L2_p2 = p2_distance.p2_L2_pt;

      distance = p2_distance.acc_dist;
      result.visited = p2_distance.estimate.visited;

      auto duration = chrono::duration_cast<chrono::microseconds>(end_clock - start_clock);
      int time = (double) duration.count();
      result.p1_time_to_L0 = time / 2;
      result.p2_time_to_L0 = time / 2;
    } else{
      result.visited = new bool[graph.nVertices];
      std::fill_n(result.visited, graph.nVertices, false);
    }
    std::set<VertexIdx> L0_p1s; // L0 points connected to p1
    std::set<VertexIdx> L0_p2s; // L0 points connected to p2
    // printf("visited status: %d\n", visited[17619650]);
    // printf("about to traverse Ln to L0\n");
    chrono::high_resolution_clock::time_point p1_start_clock = chrono::high_resolution_clock::now();
    Ln_to_L0 p1_distance = distance_Ln_to_L0(p1, L0_p1s, L2_p1, result.visited);
    // printf("traversed from p1\n");
    chrono::high_resolution_clock::time_point p1_end_clock = chrono::high_resolution_clock::now();
    chrono::high_resolution_clock::time_point p2_start_clock = chrono::high_resolution_clock::now();
    Ln_to_L0 p2_distance = distance_Ln_to_L0(p2, L0_p2s, L2_p2, result.visited);
    // printf("traversed from p2\n");
    chrono::high_resolution_clock::time_point p2_end_clock = chrono::high_resolution_clock::now();

    auto p1_duration = chrono::duration_cast<chrono::microseconds>(p1_end_clock - p1_start_clock);
    auto p2_duration = chrono::duration_cast<chrono::microseconds>(p2_end_clock - p2_start_clock);
    if(!L0[p1] && !L1[p1] ) {
      if(result.p1_time_to_L0== -1) result.p1_time_to_L0 = (double) p1_duration.count();
      else result.p1_time_to_L0 += (double) p1_duration.count();
      result.p1_dist_to_L0 = p1_distance.distance;
    }
    if(!L0[p2] && !L1[p2] ) {
      if(result.p2_time_to_L0== -1) result.p2_time_to_L0 = (double) p2_duration.count();
      else result.p2_time_to_L0 += (double) p2_duration.count();
      result.p2_dist_to_L0 = p2_distance.distance;
    }
    // printf("p1 L2: %d\n", p1_distance.distance);
    // printf("p2 L2: %d\n", p2_distance.distance);
    if(p1_distance.distance < 0 || p2_distance.distance < 0)
      return result;
    
    distance = distance + p1_distance.distance + p2_distance.distance;
    // check if P1 and P2 connect in L1; otherwise, continue to L0
    if(p1_distance.early_connection) throw std::logic_error("found early connection with P1\n"); // this should not be possible
    if(p2_distance.early_connection){ 
      // printf("early connection\n");
      distance = distance  - 2; // early connection can only happen in L1 so need to remove traversal steps to L0
      result.estimate_distance = distance;
      return result;
    }

    result.visited[p1] = false;
    result.visited[p2] = false;

    int L0_distance = 0;
    if(optimized) L0_distance = distance_L0_to_L0_Optimized(L0_p1s, L0_p2s);
    else {
      
      L0_to_L0 L0_info = distance_L0_to_L0(L0_p1s, L0_p2s);
      L0_distance = L0_info.dist;
      result.p1_L0_hops = L0_info.p1_hops;
      result.p2_L0_hops = L0_info.p2_hops;
    }
    // printf("L0 distance: %d\n", L0_distance);
    if(L0_distance < 0) throw std::logic_error("unable to find a path through L0\n");
    distance = distance + L0_distance;
    result.estimate_distance = distance;
    return result;
  }
  /**
   * @brief Computes vertex {p}'s distance to L0 and tracks newly visited vertices in {visited}
   * 
   * @param p 
   * @param L0_pts - L0 points p is connected to
   * @param L2_p - if we previously traversed {p}'s L>=2 neighborhood, then this is the L2 vertex p is connected to 
   * @param visited - boolean array where true means we've seen that vertex; this is updated in this algorithm
   * @return Ln_to_L0 {distance, early connection through L1}
   */
  Ln_to_L0 distance_Ln_to_L0(VertexIdx p, std::set<VertexIdx> &L0_pts, VertexIdx L2_p, bool* visited) const{
    // p in L0
    if(L0[p]){
      L0_pts.insert(p);
      return {0, false};
    // p in L1
    }else if (L1[p]){
      std::set<VertexIdx> L1_pts;
      L1_pts.insert(p);
      if(visited[p]) return {1, true}; // this means there's a connection at L1
      visited[p] = true;
      L1_to_L0(&L1_pts, L0_pts);
      return {1, false};
    // p in L2+
    } else{
      int distance = 0;
      if(L2_p < 0){ // means we have not seen p's L2+ neighborhood
        Ln_to_L2 l1_connection = shortest_path_to_L2(p, visited);
        if (l1_connection.distance < 0 ) return {-1, false}; // neighborhood is disconnected
        distance = l1_connection.distance;
        L2_p = l1_connection.L2_pt;
      }
      // Collect all L1 points connected to L2 endpoint
      std::set<VertexIdx> L1_pts;
      VertexIdx max_L1 = 0;
      EdgeIdx d_max_L1 = 0;
      for (EdgeIdx j = graph.offsets[L2_p]; j < graph.offsets[L2_p+1]; j++){
        VertexIdx nbor = graph.nbors[j];
        if(L1[nbor]) {
          L1_pts.insert(nbor);
           // EdgeIdx d = graph.degree(nbor);
           // if(d>d_max_L1){
           //  max_L1 = nbor;
           // d_max_L1 = d;
          // }
          if(visited[nbor]) {
            // printf("found early connection with %ld\n", nbor);
            return {distance + 2, true}; // we correct for overcounted distance in estimate_distance
          }
          visited[nbor] = true;
        }
      }
      //L1_pts.insert(max_L1);
      L1_to_L0(&L1_pts, L0_pts);
      distance = distance + 2; // plus 2 for L2 -> L1 -> L0
      //printf("      P2 -> P1 for %ld: %i\n", p, l1_connection.estimate.estimate_distance);
      return {distance, false};
    }
  }
  /**
   * @brief 
   * 
   * @param p1s 
   * @param p2s 
   * @return int 
   */
  L0_to_L0 distance_L0_to_L0(std::set<VertexIdx> &p1s, std::set<VertexIdx> &p2s) const {
    // Mark all the vertices as not visited
    bool* visited = new bool[graph.nVertices];
    std::fill_n(visited, graph.nVertices, false);
  
    // Create a queue for BFS
    std::queue<VertexIdx> q1;
    std::queue<VertexIdx> q2;
    std::unordered_map<VertexIdx, int> distance1; // map to retrace BFS (node -> parent)
    std::unordered_map<VertexIdx, int> distance2;
 
    // Mark starting nodes as visited and enqueue it
    for(auto p1 : p1s) {
      // printf("p1: %ld\n", p1);
      visited[p1] = true;
      q1.push(p1);
      distance1[p1] = 0;
    }

    for(auto p2 : p2s) {
      // printf("p2: %ld\n", p2);
      if(visited[p2]) return {0, 0, 0};
      visited[p2] = true;
      q2.push(p2);
      distance2[p2] = 0;
    }
    int p1_hops = 0;
    int p2_hops = 0;
    int mu = nVerticesL0*2;
    while(!q1.empty() && !q2.empty()) {
        // Dequeue a vertex from queue and print it
        VertexIdx s1 = q1.front();
        VertexIdx s2 = q2.front();
        if(distance1[s1] + distance2[s2]>mu){
          return {mu, p1_hops, p2_hops};
        }
        if(distance1[s1] < distance2[s2]){
           q1.pop();
           for (EdgeIdx j = graph.offsets[s1]; j < graph.offsets[s1+1]; j++){
            VertexIdx nbor = graph.nbors[j];
            if(L0[nbor]){ //only check nbors in L0
              if(!visited[nbor]){
                  visited[nbor] = true;
                  q1.push(nbor);
                  distance1[nbor] = distance1[s1]+1;
              } else if(distance1.find(nbor) == distance1.end()){
                  distance1[nbor] = distance1[s1]+1;
                  q1.push(nbor);
                  if( distance1[nbor] + distance2[nbor]< mu){
                      mu = distance1[nbor] + distance2[nbor];
                      p1_hops = distance1[nbor];
                      p2_hops = distance2[nbor];
                      if(distance2[nbor]==0) return {mu, p1_hops, p2_hops};
                  }
              }
            }
          }
        } else {
          q2.pop();
           for (EdgeIdx j = graph.offsets[s2]; j < graph.offsets[s2+1]; j++){
            VertexIdx nbor = graph.nbors[j];
            if(L0[nbor]){ //only check nbors in L0
              if(!visited[nbor]){
                  visited[nbor] = true;
                  q2.push(nbor);
                  distance2[nbor] = distance2[s2]+1;
              }else if(distance2.find(nbor) == distance2.end()){
                  distance2[nbor] = distance2[s2]+1;
                  q2.push(nbor);
                  if(distance1[nbor] + distance2[nbor]< mu){
                      mu = distance1[nbor] + distance2[nbor];
                      p1_hops = distance1[nbor];
                      p2_hops = distance2[nbor];
                      if(distance1[nbor]==0) return {mu, p1_hops, p2_hops};
                  }
              }
            }
          }
        }
    }
    if (mu!= 2*nVerticesL0) return {mu, p1_hops, p2_hops};
    return {-1, p1_hops, p2_hops};
  }
  /**
   * @brief Using PLL, xompute distance between set of p1 vertices and set of p2 vertices within L0
   * 
   * @param p1s - set of source vertices in L0
   * @param p2s - set of destination vertices in L0
   * @return int - computer distance
   */
  int distance_L0_to_L0_Optimized(std::set<VertexIdx> &p1s, std::set<VertexIdx> &p2s) {
    int shortest = graph.nVertices;

    for(auto p1 : p1s) {
      for(auto p2 : p2s){
        if(p1==p2) return 0;
        const int distance = pll.QueryDistance(p1, p2);
        if (distance<shortest && distance!=-1) return shortest=distance;
      }
    }
    return shortest;
  }
  /**
   * @brief Collects all L0 neighbors of each vertice in {L1_pts} and adds them to {L0_pts}
   * 
   * @param L1_pts 
   * @param L0_pts 
   * @return std::set<VertexIdx> 
   */
  void L1_to_L0(std::set<VertexIdx> *L1_pts, std::set<VertexIdx> &L0_pts) const {
    // VertexIdx max_L0 = 0;
    // EdgeIdx d_max_L0 = 0;
    for(auto p1: *L1_pts){
      for (EdgeIdx j = graph.offsets[p1]; j < graph.offsets[p1+1]; j++){
        VertexIdx nbor = graph.nbors[j];
        if(L0[nbor]){
          L0_pts.insert(nbor);
          // EdgeIdx d = graph.degree(nbor);
          // if(d > d_max_L0){
          //   max_L0 = nbor;
          //   d_max_L0 = d;
          // }
        }
      }
    }
    // L0_pts.insert(max_L0);
  }
  /**
   * @brief Traverses {p1}'s and {p2's} L>=2 neighborhood(s). If p1 and p2 are in the same neighborhood,
   * then returns the distance. Else, returns the L2 points they are both connected to, 
   * as well as their accumulated distances to L1. Also, returns array fo visited vertices. 
   * 
   * @param p1 
   * @param p2 
   * @return L2_to_L2_info 
   */
  L2_to_L2_info L2_to_L2(VertexIdx p1, VertexIdx p2) const {
    bool* visited = new bool[graph.nVertices];
    std::fill_n(visited, graph.nVertices, false);

    // Create a queue for BFS
    std::queue<VertexIdx> q1;
    q1.push(p1);

    std::queue<VertexIdx> q2;
    q2.push(p2);
 
    // Mark the current node as visited and enqueue it
    visited[p1] = true;
    visited[p2] = true;

    std::unordered_map<VertexIdx, int> distance1; // map to retrace BFS (node -> parent)
    std::unordered_map<VertexIdx, int> distance2;
    distance1[p1] = 0;
    distance2[p2] = 0;

    int mu = 2*graph.nVertices;

    VertexIdx p1_L2_pt = -1;
    VertexIdx p2_L2_pt = -1;

    int hop_limit_p1 = 2;
    int hop_limit_p2 = 2;

    while(!q1.empty() && !q2.empty() ){ // checks if p1 and p2 are in the same neighborhood
      VertexIdx s1 = q1.front();
      VertexIdx s2 = q2.front();
      if(distance1[s1]>hop_limit_p1 && distance2[s2] >hop_limit_p2) break;
      if(distance1[s1]+distance2[s2]>mu){
          return {-1, -1, -1, SamplingInformation(visited, mu, graph.nVertices)};
      }
      if(distance1[s1] < distance2[s2] && distance1[s1] <= hop_limit_p1){
        q1.pop();
        for (EdgeIdx j = graph.offsets[s1]; j < graph.offsets[s1+1]; j++){
          VertexIdx nbor = graph.nbors[j];
          if(L1[nbor]){ // set L2 pt if one hasn't alread been found
            if(p1_L2_pt < 0) p1_L2_pt = s1;
            continue;
          }
          if(!visited[nbor]){
            visited[nbor] = true;
            q1.push(nbor);
            distance1[nbor] = distance1[s1] + 1;
          } else if(distance1.find(nbor) == distance1.end()){
            distance1[nbor] = distance1[s1]+1;
            q1.push(nbor);
            if( distance1[nbor] + distance2[nbor]< mu){
                mu = distance1[nbor] + distance2[nbor];
                if (distance2[nbor]==0) return {-1, -1, -1, SamplingInformation(visited, mu, graph.nVertices)};
            }
          }
        }
      } else if (distance2[s2] <= hop_limit_p2) {
        q2.pop();
        for (EdgeIdx j = graph.offsets[s2]; j < graph.offsets[s2+1]; j++){
            VertexIdx nbor = graph.nbors[j];
            if(L1[nbor]){ // set L2 pt if one hasn't alread been found
              if(p2_L2_pt < 0) p2_L2_pt = s2;
              continue;
            }
            if(!visited[nbor]){
              visited[nbor] = true;
              q2.push(nbor);
              distance2[nbor] = distance2[s2]+1;
            } else if(distance2.find(nbor) == distance2.end()){
              distance2[nbor] = distance2[s2]+1;
              q2.push(nbor);
              if(distance1[nbor] + distance2[nbor]< mu){
                  mu = distance1[nbor] + distance2[nbor];
                  if (distance1[nbor]==0) return {-1, -1, -1, SamplingInformation(visited, mu, graph.nVertices)};
              }
            } 
        }
      }
    }
    if(mu != 2*graph.nVertices) return {-1, -1, -1, SamplingInformation(visited, mu, graph.nVertices)};
    if(p1_L2_pt < 0 && !q1.empty()){ // continue traversing p1's neighborhood if we haven't hit L2 yet
      p1_L2_pt = shortest_path_to_L2(&distance1, visited, &q1);
    }
    if (p2_L2_pt < 0 && !q2.empty()){ // continue traversing p2's neighborhood if we haven't hit L2 yet
      p2_L2_pt = shortest_path_to_L2(&distance2, visited, &q2);
    }
    if(p1_L2_pt < 0 && p2_L2_pt < 0) return {-1, -1, -1, SamplingInformation(visited, graph.nVertices)}; // means path is not possible
    return {p1_L2_pt, p2_L2_pt, distance1[p1_L2_pt] + distance2[p2_L2_pt], SamplingInformation(visited, graph.nVertices)};
  }
  /**
   * @brief Computes the distance between p1 and L2. Also tracks the visited vertices. 
   * 
   * @param p1 
   * @param visited 
   * @return Ln_to_L2 storing {L2 endpoint, distance} or {-1, -1} if no path was found
   */
  Ln_to_L2 shortest_path_to_L2(VertexIdx p1, bool* visited) const {
    std::queue<VertexIdx> q;
    visited[p1] = true;
    q.push(p1);
    std::unordered_map<VertexIdx, int> distance;
    distance[p1] = 0;

    VertexIdx L2_pt = shortest_path_to_L2(&distance, visited, &q);

    int d = -1;
    if(L2_pt != -1) d = distance[L2_pt];
    return {L2_pt, d};
  }

  /**
   * @brief Computes the distance between p1 and L2, contiuing the BiBFS from {prev_distance}, {prev_q}, and {prev_visited}. 
   * Also tracks the visited vertices. 
   * 
   * @param prev_distance
   * @param prev_visited
   * @param prev_q 
   * @return VertexIdx L2 point first reached
   */
  VertexIdx shortest_path_to_L2(std::unordered_map<VertexIdx, int> *prev_distance, bool* prev_visited, std::queue<VertexIdx> *prev_q) const {
    bool* visited = prev_visited;

    std::queue<VertexIdx> *q = prev_q;
    std::unordered_map<VertexIdx, int> *distance = prev_distance; // map to retrace BFS (node -> parent)
    while(!q->empty()){
        VertexIdx s = q->front();
        q->pop();
        for (EdgeIdx j = graph.offsets[s]; j < graph.offsets[s+1]; j++){
            VertexIdx nbor = graph.nbors[j];
            if (L1[nbor]){
                return s;
            }
            if(!visited[nbor]){
                visited[nbor] = true;
                q->push(nbor);
                (*distance)[nbor] = (*distance)[s] + 1;
            }
            // Check if L1 has been reached
        }
    }
    return -1;
  }
  
  // void L2_stats(std::string graph_name) const {
  //   bool *visited_L2 = new bool[graph.nVertices];
  //   for(VertexIdx i = 0; i < graph.nVertices; i++){
  //     if(L0[i] || L1[i] || visited_L2[i]) continue;
  //     std::queue<VertexIdx> q;
  //     visited_L2[i] = true;
  //     q.push(i);
  //     std::unordered_map<VertexIdx, int> distance;
  //     distance[i] = 0;
  //     shortest_path_to_L2(&distance, visited_L2, &q, true);
  //     for(VertexIdx e; )
  //   }
  // }
  /**
   * @brief Print L0 and L1 boolean arrays to {f}. Default is to stdout. 
   * 
   * @param f 
   */
  void print(FILE* f = stdout) const {
    graph.print(f);
    fprintf(f, "L0 Boolean Array: \n");
    for (int i = 0; i < graph.nVertices; ++i)
    {
        fprintf(f, L0[i]? "0 " : "1 ");
    }
    fprintf(f, "\n");
    fprintf(f, "L1 Boolean Array: \n");
    for (int i = 0; i < graph.nVertices; ++i)
    {
        fprintf(f, L1[i]? "0 " : "1 ");
    }
    fprintf(f, "\n");


  }
  // depracated 
  void generate_histogram(std::vector<double> degrees, std::ostream& DegreeFile, std::string name)const {
    DegreeFile << "PLOT "<< name << std::endl;
  int max_value = *max_element(degrees.begin(), degrees.end());
  const int num_bins = 6;
  int bin_size = max_value / num_bins;
  int bins[num_bins] = {0};
  DegreeFile << "DEGREES ";
  for(int i=0; i < degrees.size(); i++){
    DegreeFile << degrees[i] <<",";
    int bin_number = int(degrees[i]/bin_size);
    if(bin_number==num_bins){
        bins[bin_number-1]++;
    } else{
      bins[bin_number]++;
    }
  }
  DegreeFile << std::endl;
  for(int i = 0; i< num_bins; i++){
    DegreeFile << " bin " << i << "[ " << i*bin_size <<", "<<(i+1)*bin_size<<" ): " << bins[i]<<std::endl;
  }
}
};
}

#endif