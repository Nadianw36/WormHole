#ifndef ESCAPE_L0GRAPH_H_
#define ESCAPE_L0GRAPH_H_

#include <cstdlib>
#include <stdexcept>
#include <cstdio>
#include <numeric>
#include <algorithm>
#include <set>
#include <vector>
#include <utility>
#include <unordered_map>
#include <string>
#include <iostream>
#include <chrono>
#include "Escape/Graph.h"
#include "Escape/BiBFSGraph.h"
#include "Escape/Config.h"
#include "Escape/PrunedLandmark.h"
#include "Escape/IndexedBinaryHeap.h"
#include <queue>

using namespace std;

namespace Escape
{

  // data structure for storing sampling information from Ln to L2 path
  struct Ln_to_L2
  {
    VertexIdx L2_pt;
    int distance;
  };
  struct Ln_to_L0
  {
    int distance;
    bool early_connection;
  };

  std::string extractGraphNameFromL0Name(std::string graphL0Name)
  {
    std::string delimiter = "_";
    std::size_t d_pos = graphL0Name.find(delimiter);
    std::string graphName = graphL0Name.substr(0, d_pos);
    return graphName;
  }

  class L0Graph : public BiBFSGraph
  {
  public:
    bool *L0; // bitmap of boolean flags
    bool *L1; // bitmap of boolean flags
    VertexIdx nVerticesL0;
    PrunedLandmarkLabeling<> pll;
    ofstream output;
    chrono::high_resolution_clock::time_point outside_core_start_clock, outside_core_end_clock;
    VertexIdx randomL0V1, randomL0V2, randomL0FromRandomL1V1, randomL0FromRandomL1V2;

    // shared variables
    int acc_distance, L0_distance;
    VertexIdx p1_L2_pt, p2_L2_pt;
    bool earlyConnection;
    bool L0traversal; // true means BiBFS through L0 else BiBFS through L2

    L0Graph() {}
    /**
     * @brief Construct a new L0Graph object from binary files containing L0 and L1 boolean arrays.
     * Files should be named and have path bin/{graphName}_L0.bin and bin/{graphName}_L1.bin
     *
     * @param g
     * @param graphName
     */
    L0Graph(const BiBFSGraph &graph, string graphName) : BiBFSGraph(graph)
    {

      string graphNameStem = graphName;
      string delimiter = "_";
      size_t d_pos = graphName.find(delimiter);
      if (d_pos != string::npos)
        graphNameStem = graphName.substr(0, d_pos);
      auto vcount = nVertices;
      // printf("vcount: %ld \n", vcount);
      L0 = new bool[vcount];
      L1 = new bool[vcount];

      string graph_L0_folder = BIN_FOLDER + graphNameStem + "/" + L0_FOLDER;

      ifstream L0File(graph_L0_folder + graphName + "_L0.bin", ios::in | ios::binary);
      ifstream L1File(graph_L0_folder + graphName + "_L1.bin", ios::in | ios::binary);
      L0File.read((char *)L0, sizeof(bool) * (int64_t)vcount);
      L1File.read((char *)L1, sizeof(bool) * (int64_t)vcount);
      L0File.close();
      L1File.close();

      nVerticesL0 = count(L0, L0 + vcount, true);
      // printf("L0 size: %ld \n", nVerticesL0);

      VertexIdx nVerticesL1 = count(L1, L1 + vcount, true);
      // printf("L1 size: %ld \n", nVerticesL1);
    }
    /**
     * @brief Construct a new L0Graph object via algorithm with an L0 size of sizeL0
     *
     * @param g
     * @param sizeL0 - number of vertices in L0
     * @param optimal - if True, then use high degree vertice as seed vertice, else use random vertice
     */
    L0Graph(const CGraph &graph, VertexIdx sizeL0, bool optimal = false, bool usePriorityQueue = false) : BiBFSGraph(graph)
    {
      output.open("output/output.txt");
      auto vcount = nVertices;
      nVerticesL0 = sizeL0;

      L0 = new bool[vcount];
      L1 = new bool[vcount];

      fill_n(L0, vcount, false);
      fill_n(L1, vcount, false);

      if (usePriorityQueue)
        generateL0FromPriority();
      else
        generateL0FromGreedy(optimal);
    }
    /**
     * @brief Assert that:
     *            1) for every vertex in L0, all of its neighbors are in L1
     *            2) for every vertex in L1, at least one neighbor is in L0
     */
    void checkForBadL0() const
    {
      for (VertexIdx v = 0; v < nVertices; v++)
      {
        if (L0[v])
        {
          for (EdgeIdx e = offsets[v]; e < offsets[v + 1]; e++)
          {
            VertexIdx nbor = nbors[e];
            if (!(L0[nbor] || L1[nbor]))
            {
              printf("%ld %ld \n", v, nbor);
              throw invalid_argument("all neighbors of L0 should be in L0 or L1");
            }
          }
        }
        else if (L1[v])
        {
          int countL0 = 0;
          for (EdgeIdx e = offsets[v]; e < offsets[v + 1]; e++)
          {
            VertexIdx nbor = nbors[e];
            if (L0[nbor])
              countL0++;
          }
          if (!(countL0 > 0))
          {
            for (EdgeIdx e = offsets[v]; e < offsets[v + 1]; e++)
            {
              VertexIdx nbor = nbors[e];
              printf("nbor %ld \n", nbor);
            }
            printf(" %ld %ld %d %d\n", v, countL0, L1[v], L0[v]);
            throw invalid_argument("L1 vertices are all connected to L0 %ld ");
          }
        }
      }
    }

    /**
     * @brief Construct a new L0Graph object via algorithm with an L0 size of {nVerticesL0}
     */
    void generateL0FromGreedy(bool optimal)
    {
      IndexedBinaryHeap L1_map = IndexedBinaryHeap(nVertices); // maps nodes currently in L1 -> L0 degree

      auto v = rand() % (VertexIdx)nVertices; // picks seed node
      if (optimal)
      {
        VertexIdx *p = new VertexIdx[nVertices];
        VertexIdx sampleThresh = (VertexIdx)(0.01 * nVertices);
        for (VertexIdx v = 0; v < nVertices; v++)
        {
          p[v] = v;
        }
        sort(p, p + nVertices, [this](VertexIdx v1, VertexIdx v2)
             { return degree(v1) > degree(v2); });
        VertexIdx sample = rand() % (int16_t)sampleThresh;
        v = p[sample];
        delete p;
      }

      L0[v] = true;

      for (EdgeIdx i = offsets[v]; i < offsets[v + 1]; i++)
      {
        VertexIdx nbor = nbors[i];
        L1_map.increment(nbor);
      }
      for (int i = 1; i < nVerticesL0; i++)
      {

        auto max_v = L1_map.popMax();
        L0[max_v] = true;
        for (EdgeIdx j = offsets[max_v]; j < offsets[max_v + 1]; j++)
        {
          VertexIdx nbor = nbors[j];
          if (!L0[nbor])
          {
            L1_map.increment(nbor);
          }
        }
      }
      for (VertexIdx i = 0; i < nVertices; i++)
      {
        if (L1_map.keyToPosition[i] > 0)
          L1[i] = true;
      }
    }

    /**
     * @brief Create L0 from the {nVerticesL0} highest degree vertices in graph.
     * If resulting L0 is disconnected, then make L0 the largest component
     * (meaning that the resulting L0 might have less vertices than {nVerticesL0})
     */
    void generateL0FromPriority()
    {
      priority_queue<pair<EdgeIdx, VertexIdx>> p;
      for (VertexIdx v = 0; v < nVertices; v++)
      {
        EdgeIdx degrees = degree(v);
        p.push(make_pair(degrees, v));
      }
      // printf("generated pq\n");
      for (VertexIdx i = 0; i < nVerticesL0; i++)
      {
        VertexIdx vl0 = p.top().second;
        p.pop();
        L0[vl0] = true;
        L1[vl0] = false;
        for (EdgeIdx e = offsets[vl0]; e < offsets[vl0 + 1]; e++)
        {
          VertexIdx nbor = nbors[e];
          if (!L0[nbor])
          {
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
    VertexIdx connectComponentsL0(bool prune = true)
    {
      int count = 0;
      VertexIdx num_nodes = 0;
      VertexIdx largest_component = 0;
      VertexIdx *component_size = new VertexIdx[nVertices];
      fill_n(component_size, nVertices, 0);

      bool *visited = new bool[nVertices];
      fill_n(visited, nVertices, false);
      // find components
      for (VertexIdx i = 0; i < nVertices; i++)
      {
        if (!L0[i] || visited[i])
          continue;
        count++;
        VertexIdx seed = i;
        queue<VertexIdx> q;
        q.push(seed);
        visited[seed] = true;
        int size = 0;
        while (!q.empty())
        {
          num_nodes++;
          VertexIdx v = q.front();
          size++;
          q.pop();
          for (EdgeIdx e = offsets[v]; e < offsets[v + 1]; e++)
          {
            VertexIdx nbor = nbors[e];
            if (L0[nbor] && !visited[nbor])
            {
              visited[nbor] = true;
              q.push(nbor);
            }
          }
        }
        component_size[i] = size;
        if (size > largest_component)
          largest_component = size;
        // printf("size of connected component: %d\n", size);
        // printf("L0_size: %d\n", num_nodes);
      }
      // printf("number of components: %d\n", count);
      // printf("L0_size: %d\n", num_nodes);
      // printf("largest component: %d\n", largest_component);
      // prune
      for (VertexIdx i = 0; i < nVertices; i++)
      {
        if (!L0[i] || component_size[i] == largest_component || component_size[i] == 0)
          continue;
        count--;
        VertexIdx seed = i;
        queue<VertexIdx> q;
        q.push(seed);
        L0[seed] = false;
        while (!q.empty())
        {
          VertexIdx v = q.front();
          num_nodes--;
          q.pop();
          for (EdgeIdx e = offsets[v]; e < offsets[v + 1]; e++)
          {
            VertexIdx nbor = nbors[e];
            if (L0[nbor])
            {
              q.push(nbor);
              L0[nbor] = false;
            }
            else if (L1[nbor])
              L1[nbor] = false;
          }
        }
      }
      // printf("number of components: %d\n", count);
      return num_nodes;
    }
    vector<pair<VertexIdx, VertexIdx>> L0_edgelist()
    {
      vector<pair<VertexIdx, VertexIdx>> edgelist;
      for (VertexIdx i = 0; i < nVertices; i++)
      {
        if (!L0[i])
          continue;
        for (EdgeIdx j = offsets[i]; j < offsets[i + 1]; j++)
        {
          VertexIdx nbor = nbors[j];
          if (L0[nbor])
          {
            pair<VertexIdx, VertexIdx> edge(i, nbor);
            edgelist.push_back(edge);
          }
        }
      }
      return edgelist;
    }

    /**
     * @brief Set the L0 PLL object
     * Constructs a Pruned Landmark Labelling Graph in L0 to speed up distance computing
     */
    void setL0_PLL()
    {
      vector<pair<int, int>> edgelist;
      for (VertexIdx i = 0; i < nVertices; i++)
      {
        if (L0[i])
        {
          for (EdgeIdx j = offsets[i]; j < offsets[i + 1]; j++)
          {
            VertexIdx nbor = nbors[j];
            if (L0[nbor])
            {
              pair<int, int> edge((int64_t)i, (int64_t)nbor);
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
    void writeGraphToFile(string graphName)
    {
      std::replace(graphName.begin(), graphName.end(), '.', '-');

      std::string delimiter = "_";
      std::size_t d_pos = graphName.find(delimiter);
      std::string graphFolder = graphName.substr(0, d_pos) + "/";

      string graph_L0_folder = BIN_FOLDER + graphFolder + L0_FOLDER;
      ofstream L0File(graph_L0_folder + graphName + "_L0.bin", ios::out | ios::binary);
      ofstream L1File(graph_L0_folder + graphName + "_L1.bin", ios::out | ios::binary);
      L0File.write((char *)L0, sizeof(bool) * (int64_t)nVertices);
      L1File.write((char *)L1, sizeof(bool) * (int64_t)nVertices);
      L0File.close();
      L1File.close();
    }
    /**
     * @brief Print L0 -> L0, L0 -> L1, L1 -> L0, and L1 -> L1 degrees, as well as the total degree file
     * in results/ folder.
     *
     * @param graph_name
     */
    void print_size(ostream &graph_file, string graphL0Name) const
    {
      std::replace(graphL0Name.begin(), graphL0Name.end(), '.', '-');

      std::string graphFolder = extractGraphNameFromL0Name(graphL0Name);
      std::string fileLocationAndPrefix = RESULTS_FOLDER + graphFolder + "/" + L0_FOLDER + DEGREES_FOLDER + graphL0Name;
      ofstream degreeFile(fileLocationAndPrefix + "_degrees.txt");
      ofstream L0DegreeFile(fileLocationAndPrefix + "_degreesL0.txt");
      ofstream L1DegreeFile(fileLocationAndPrefix + "_degreesL1.txt");
      ofstream L0L1DegreeFile(fileLocationAndPrefix + "_degreesL0L1.txt");
      ofstream L1L0DegreeFile(fileLocationAndPrefix + "_degreesL1L0.txt");
      ofstream L2L1DegreeFile(fileLocationAndPrefix + "_degreesL2L1.txt");

      vector<int> degrees;
      vector<int> L0_degrees;
      vector<int> L1_degrees;
      vector<int> L1_L0_degrees;
      vector<int> L0_L1_degrees;
      vector<int> L2_L1_degrees;
      int mistakes = 0;
      for (VertexIdx i = 0; i < nVertices; i++)
      {
        if (L0[i] && L1[i])
          mistakes++;
        degrees.push_back(offsets[i + 1] - offsets[i]);

        int L0_degree = 0;
        int L0_L1_degree = 0;
        int L1_degree = 0;
        int L1_L0_degree = 0;
        int L2_L1_degree = 0;

        for (EdgeIdx j = offsets[i]; j < offsets[i + 1]; j++)
        {
          VertexIdx nbor = nbors[j];
          if (L0[i] && L0[nbor])
          {
            L0_degree++;
          } // L0 -> L0
          else if (L0[i] && L1[nbor])
          {
            L0_L1_degree++;
          } // L0 -> L1
          else if (L1[i] && L1[nbor])
          {
            L1_degree++;
          } // L1 -> L1
          else if (L1[i] && L0[nbor])
          {
            L1_L0_degree++;
          } // L1 -> L0
          else if (!L0[i] && !L1[i] && L1[nbor])
          {
            L2_L1_degree++;
          } // L2 -> L1
        }

        if (L0[i])
        {
          L0_degrees.push_back(L0_degree);
          L0_L1_degrees.push_back(L0_L1_degree);
        }
        else if (L1[i])
        {
          L1_degrees.push_back(L1_degree);
          L1_L0_degrees.push_back(L1_L0_degree);
        }
        else
        {
          L2_L1_degrees.push_back(L2_L1_degree);
        }
      }
      for (const auto &e : degrees)
        degreeFile << e << "\n";
      for (const auto &e : L0_degrees)
        L0DegreeFile << e << "\n";
      for (const auto &e : L1_degrees)
        L1DegreeFile << e << "\n";
      for (const auto &e : L0_L1_degrees)
        L0L1DegreeFile << e << "\n";
      for (const auto &e : L1_L0_degrees)
        L1L0DegreeFile << e << "\n";
      for (const auto &e : L2_L1_degrees)
        L2L1DegreeFile << e << "\n";

      graph_file << graphL0Name << " size L0: " << L0_degrees.size() << "\n";
      graph_file << graphL0Name << " size L1: " << L1_degrees.size() << "\n";
    };

    /** MAIN ALGORITHM
     * @brief Run L0 algorithm to estimate distance between p1 and p2
     *
     * @param p1
     * @param p2
     * @param optimized - If true, use PLL in L0; else use BiBFS
     * @return SamplingInformation which contains estimated distance as well as visited vertices (excluding visited L0 vertices)
     */
    int estimate_distance(VertexIdx p1, VertexIdx p2, bool optimized)
    {
      reset();
      p1_L2_pt = -1;
      p2_L2_pt = -1;
      int pathDistance = 0;
      earlyConnection = false;

      // printf("vertices: %d, %d\n", p1, p2);
      // printf("L0: %d, %d\n", L0[p1], L0[p2]);
      // printf("L1: %d, %d\n", L1[p1], L1[p2]);
      if (p1 == p2)
        return 0;

      outside_core_start_clock = chrono::high_resolution_clock::now();

      if (!L0[p1] && !L1[p1] && !L0[p2] && !L1[p2])
      { // BiBFS within L2 neighborhood(s)
        int L2_distance = L2_to_L2(p1, p2);
        // printf("L2->L2: %d, p1 L2 %ld %d, p2 L2 %ld %d\n", distance, p1_L2_pt, distance[p1_L2_pt], p2_L2_pt, distance[p2_L2_pt]);

        if (L2_distance > 0 || (p1_L2_pt < 0 || p2_L2_pt < 0))
        {
          earlyConnection = true;
          return L2_distance; // found path or path not possible
        }
      }

      Ln_to_L0 p1_distance = distance_Ln_to_L0(p1, p1_L2_pt, queue1, 1);
      // printf("p1 dist %d ec %d\n", p1_distance.distance, p1_distance.early_connection);

      Ln_to_L0 p2_distance = distance_Ln_to_L0(p2, p2_L2_pt, queue2, -1);
      // printf("p2 dist %d ec %d\n", p2_distance.distance, p2_distance.early_connection);

      outside_core_end_clock = chrono::high_resolution_clock::now();

      if (p1_distance.distance < 0 || (!p2_distance.early_connection && p2_distance.distance < 0))
        return -1;

      pathDistance += p1_distance.distance + p2_distance.distance;

      if (p2_distance.early_connection)
      {
        earlyConnection = true;
        return pathDistance;
      }

      if (optimized)
        pathDistance += distance_L0_to_L0_Optimized();
      else
      {
        int L0dist = distance_L0_to_L0();
        // printf("l0 dist %ld \n", L0dist);
        pathDistance += L0dist;
      }
      return pathDistance;
    }

    /**
     * @brief Computes vertex {p}'s distance to L0 and tracks newly visited vertices in {visited}
     *
     * @param p
     * @param L0_pts - L0 points p is connected to
     * @param L2_p - if we previously traversed {p}'s L>=2 neighborhood, then this is the L2 vertex p is connected to
     * @param queue - queue to load L0 vertices on
     * @param dist - should be 1 if from side 1; otherwise -1 (from side 2)
     * @return Ln_to_L0 {distance, early connection through L1}
     */
    Ln_to_L0 distance_Ln_to_L0(VertexIdx p, VertexIdx &L2_p, VertexIdx *queue, int dist)
    {
      resetQueue(queue);

      if (L0[p])
      {
        if (visited[p])
        {
          mergeVertex = p;
          secondParent = p;
          return {0, true};
        }

        queueVertex(p, -1, queue, dist);
        return {0, false};
      }
      else if (L1[p])
      {
        if (visited[p])
        {
          mergeVertex = p;
          secondParent = p;
          return {-1, true}; // can only happen during second point traversal and we need to undo steps to L0
        }

        visited[p] = true;
        parent[p] = -1;

        VertexIdx nbor;
        for (EdgeIdx j = offsets[p]; j < offsets[p + 1]; j++)
        {
          nbor = nbors[j];

          if (L0[nbor])
          {
            if (visited[nbor])
            {
              mergeVertex = nbor;
              secondParent = p;
              return {1, true};
            }
            queueVertex(nbor, p, queue, dist);
          }
        }

        return {1, false};
      }
      else
      {
        int pathDistance = 0;
        if (L2_p < 0)
        { // means we have not seen p's L2+ neighborhood
          L2_p = shortest_path_to_L2(p, queue, dist);
          if (L2_p < 0)
            return {-1, false}; // neighborhood is disconnected
        }
        pathDistance = (dist * distance[L2_p]) - 1;
        // printf("L2_p %ld distance %d %d %d\n", L2_p, pathDistance, distance[L2_p], distance[p]);

        resetQueue(queue);

        VertexIdx nbor;
        VertexIdx L0nbor;

        for (EdgeIdx j = offsets[L2_p]; j < offsets[L2_p + 1]; j++)
        {
          nbor = nbors[j];
          if (L1[nbor])
          {
            // printf("   neighbor in L1 %ld %d\n", nbor, visited[nbor]);
            if (visited[nbor])
            {
              mergeVertex = nbor;
              secondParent = p;
              return {pathDistance, true}; // see line 595
            }
            visited[nbor] = true;
            parent[nbor] = L2_p;
            for (EdgeIdx i = offsets[nbor]; i < offsets[nbor + 1]; i++)
            {
              L0nbor = nbors[i];
              if (L0[L0nbor])
              {
                // printf("      neighbor in L0 %ld %d\n", L0nbor, visited[L0nbor]);
                if (visited[L0nbor] && distance[L0nbor] != dist)
                {
                  mergeVertex = L0nbor;
                  secondParent = nbor;
                  return {pathDistance + 2, true};
                }

                queueVertex(L0nbor, nbor, queue, dist);
              }
            }
          }
        }
        pathDistance = pathDistance + 2; // plus 2 for L2 -> L1 -> L0
        return {pathDistance, false};
      }
    }
    bool continueSearch()
    {
      return !collision && !Q1empty() && !Q2empty() &&
             (L0traversal || (level1Record[level1Record[LEVEL_IDX]] <= 3 && level2Record[level2Record[LEVEL_IDX]] <= 3));
    }

    bool exploreVertexQ1(VertexIdx nbor, VertexIdx parent)
    {
      if (L0traversal)
        return L0[nbor];
      else
      {
        if (L1[nbor])
        {
          if (p1_L2_pt < 0)
            p1_L2_pt = parent;
          return false;
        }
        return true;
      }
    }

    bool exploreVertexQ2(VertexIdx nbor, VertexIdx parent)
    {
      if (L0traversal)
        return L0[nbor];
      else
      {
        if (L1[nbor])
        {
          if (p2_L2_pt < 0)
            p2_L2_pt = parent;
          return false;
        }
        return true;
      }
    }

    int distance_L0_to_L0(VertexIdx v1, VertexIdx v2)
    {
      resetQueue(queue1);
      queueVertex(v1, -1, queue1, 1);

      resetQueue(queue2);
      queueVertex(v2, -1, queue2, -1);

      return distance_L0_to_L0();
    }

    int distance_L0_to_L0()
    {

      level1Record[LEVEL_START] = queue1[Q_END_IDX];
      level1Record[LEVEL_IDX] = LEVEL_START;

      level2Record[LEVEL_START] = queue2[Q_END_IDX];
      level2Record[LEVEL_IDX] = LEVEL_START;

      L0traversal = true;
      // printf("check L0 state \n");
      // checkBiBFSState();

      return bidirectionalBFS();
    }

    /**
     * @brief Using PLL, xompute distance between set of p1 vertices and set of p2 vertices within L0
     *
     * @param p1s - set of source vertices in L0
     * @param p2s - set of destination vertices in L0
     * @return int - computer distance
     */
    int distance_L0_to_L0_Optimized()
    {
      int shortest = nVertices;
      VertexIdx p1, p2;

      for (VertexIdx i = queue1[queue1[Q_START_IDX]]; i < queue1[queue1[Q_END_IDX]]; i++)
      {
        p1 = queue1[i];
        for (VertexIdx j = queue2[queue2[Q_START_IDX]]; j < queue2[queue2[Q_END_IDX]]; j++)
        {
          p2 = queue2[i];
          if (p1 == p2)
            return 0;
          const int distance = pll.QueryDistance(p1, p2);
          if (distance < shortest && distance != -1)
            return shortest = distance;
        }
      }
      return shortest;
    }

    /**
     * @brief Traverses {p1}'s and {p2's} L>=2 neighborhood(s). If p1 and p2 are in the same neighborhood,
     * then returns the distance. Else, returns the L2 points they are both connected to,
     * as well as their accumulated distances to L1. Also, returns array fo visited vertices.
     *
     * @param p1
     * @param p2
     * @return int distance
     */
    int L2_to_L2(VertexIdx p1, VertexIdx p2)
    {
      L0traversal = false;

      queue1[Q_START_IDX] = Q_START;
      queue1[Q_END_IDX] = Q_START;
      queueVertex(p1, -1, queue1, 1);

      queue2[Q_START_IDX] = Q_START;
      queue2[Q_END_IDX] = Q_START;
      queueVertex(p2, -1, queue2, -1);

      level1Record[LEVEL_START] = queue1[Q_END_IDX];
      level1Record[LEVEL_IDX] = LEVEL_START;

      level2Record[LEVEL_START] = queue2[Q_END_IDX];
      level2Record[LEVEL_IDX] = LEVEL_START;

      p1_L2_pt = -1;
      p2_L2_pt = -1;

      int pathDistance = bidirectionalBFS();

      if (pathDistance > 0)
        return pathDistance;
      if (p1_L2_pt < 0 && !Q1empty())
      { // continue traversing p1's neighborhood if we haven't hit L2 yet
        p1_L2_pt = shortest_path_to_L2(queue1, 1);
      }
      if (p2_L2_pt < 0 && !Q2empty())
      { // continue traversing p2's neighborhood if we haven't hit L2 yet
        p2_L2_pt = shortest_path_to_L2(queue2, -1);
      }
      return -1; // means path is not possible
    }

    /**
     * @brief Computes the distance between p1 and L2. Also tracks the visited vertices.
     *
     * @param p1
     * @param visited
     * @return VertexIdx L2 point first reached
     */
    VertexIdx shortest_path_to_L2(VertexIdx p, VertexIdx *queue, int distanceIncrement)
    {
      queue[Q_START_IDX] = Q_START;
      queue[Q_END_IDX] = Q_START;
      queueVertex(p, -1, queue, distanceIncrement);
      // printf("q start + end %ld %ld %ld\n", queue[Q_START_IDX], queue[Q_END_IDX], queue[queue[Q_START_IDX]]);

      return shortest_path_to_L2(queue, distanceIncrement);
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
    VertexIdx shortest_path_to_L2(VertexIdx *queue, int distanceIncrement) // distances are negative from one side
    {
      VertexIdx v;
      while (queue[Q_START_IDX] != queue[Q_END_IDX])
      {
        v = queue[queue[Q_START_IDX]++];
        for (EdgeIdx j = offsets[v]; j < offsets[v + 1]; j++)
        {
          VertexIdx nbor = nbors[j];
          if (L1[nbor])
          {
            return v;
          }
          if (!visited[nbor])
          {
            queueVertex(nbor, v, queue, distance[v] + distanceIncrement);
          }
        }
      }
      return -1;
    }

    /**
     * @brief Print L0 and L1 boolean arrays to {f}. Default is to stdout.
     *
     * @param f
     */
    void print(FILE *f = stdout) const
    {
      print(f);
      fprintf(f, "L0 Boolean Array: \n");
      for (int i = 0; i < nVertices; ++i)
      {
        fprintf(f, L0[i] ? "0 " : "1 ");
      }
      fprintf(f, "\n");
      fprintf(f, "L1 Boolean Array: \n");
      for (int i = 0; i < nVertices; ++i)
      {
        fprintf(f, L1[i] ? "0 " : "1 ");
      }
      fprintf(f, "\n");
    }

    void checkBiBFSState()
    {
      for (VertexIdx i = queue1[Q_START_IDX]; i < queue1[Q_END_IDX]; i++)
      {
        // printf("%ld \n", queue1[i]);
      }

      for (VertexIdx i = queue2[Q_START_IDX]; i < queue2[Q_END_IDX]; i++)
      {
        // printf("%ld \n", queue2[i]);
      }

      for (VertexIdx i = 0; i < nVertices; i++)
      {
        if (L0[i])
        {
          if (distance[i] != 1 && distance[i] != -1 && distance[i] != INFI)
          {
            // printf("L0 vertex has distance %ld %ld\n", i, visited);
          }
          else if (distance[i] == INFI && visited[i])
          {
            // printf("L0 vertex has been visited %ld %ld\n", i, visited);
          }
        }
      }
    }
    /**
     * sample random L0 vertices --> assumes that L0 has already been loaded on queue
     * so must call distanceToCore first
     */
    void sampleL0Vertices(VertexIdx v1, VertexIdx v2)
    {
      srand(time(NULL));

      randomL0FromRandomL1V1 = -1;
      randomL0FromRandomL1V2 = -1;

      randomL0V1 = -1;
      randomL0V2 = -1;

      if (earlyConnection)
        return;

      if (!L0[v1])
      {
        // printf("looking for randomL0V1\n");
        VertexIdx sample1 = (rand() % (level1Record[LEVEL_START] - Q_START)) + Q_START;
        randomL0V1 = queue1[sample1];
        randomL0FromRandomL1V1 = randomL0V1;
        // printf(" sampling for random L0 v %ld %ld %ld\n", (level1Record[LEVEL_START] - Q_START), sample1, queue1[sample1]);
      }
      else
      {
        randomL0V1 = v1;
        randomL0FromRandomL1V1 = v1;
      }

      if (!L0[v2])
      {
        // printf("looking for randomL0V2\n");
        VertexIdx sample2 = (rand() % (level2Record[LEVEL_START] - Q_START)) + Q_START;
        randomL0V2 = queue2[sample2];
        randomL0FromRandomL1V2 = randomL0V2;
        // printf("  sampling for random L0 v %ld %ld %ld\n", (level2Record[LEVEL_START] - Q_START), sample2, queue2[sample2]);
      }
      else
      {
        randomL0V2 = v2;
        randomL0FromRandomL1V2 = v2;
      }
      VertexIdx tempRandomL0FromRandomL1V1 = sampleL0VertexFromL1(p1_L2_pt);
      VertexIdx tempRandomL0FromRandomL1V2 = sampleL0VertexFromL1(p2_L2_pt);
      // printf("%ld %ld\n", tempRandomL0FromRandomL1V1, tempRandomL0FromRandomL1V2);

      if (tempRandomL0FromRandomL1V1 > 0)
        randomL0FromRandomL1V1 = tempRandomL0FromRandomL1V1;
      if (tempRandomL0FromRandomL1V2 > 0)
        randomL0FromRandomL1V2 = tempRandomL0FromRandomL1V2;

      if ((randomL0FromRandomL1V1 > 0 && !L0[randomL0FromRandomL1V1]) ||
          (randomL0FromRandomL1V2 > 0 && !L0[randomL0FromRandomL1V2]) ||
          (randomL0V1 > 0 && !L0[randomL0V1]) ||
          (randomL0V2 > 0 && !L0[randomL0V2]))
        printf("fucked %ld %ld %d\n", v1, v2, earlyConnection);
    }

    VertexIdx sampleL0VertexFromL1(VertexIdx L2_p)
    {

      if (L2_p > 0)
      {
        // printf(" sampling for L1 v %ld\n", L2_p);
        VertexIdx *L1Vertices = new VertexIdx[degree(L2_p)];
        VertexIdx endL1Vertices = 0;
        for (EdgeIdx j = offsets[L2_p]; j < offsets[L2_p + 1]; j++)
        {
          VertexIdx nbor = nbors[j];
          if (L1[nbor])
            L1Vertices[endL1Vertices++] = nbor;
        }

        VertexIdx sample = rand() % endL1Vertices;
        // printf(" sampling for L1 v %ld %ld %ld\n", endL1Vertices, sample, L1Vertices[sample]);
        VertexIdx L1VertexSample = L1Vertices[sample];

        VertexIdx *L0Vertices = new VertexIdx[degree(L1VertexSample)];
        VertexIdx endL0Vertices = 0;

        for (EdgeIdx j = offsets[L1VertexSample]; j < offsets[L1VertexSample + 1]; j++)
        {
          VertexIdx nbor = nbors[j];
          if (L0[nbor])
            L0Vertices[endL0Vertices++] = nbor;
        }

        sample = rand() % endL0Vertices;
        // printf(" sampling for L0 v %ld %ld %ld\n", endL0Vertices, sample, L0Vertices[sample]);

        return L0Vertices[sample];
      }
      return -1;
    }

    EdgeIdx getCoreEdgeCount()
    {
      EdgeIdx nEdgesCore = 0;
      for (VertexIdx i = 0; i < nVertices; i++)
      {
        if (!L0[i])
          continue;
        for (EdgeIdx j = offsets[i]; j < offsets[i + 1]; j++)
        {
          VertexIdx nbor = nbors[j];
          if (L0[nbor])
            nEdgesCore++;
        }
      }
    }

    void writeCoreCOO(std::string graphName)
    {
      ofstream L0File(GRAPH_FOLDER + graphName + "_core-COO.txt", ios::out | ios::binary);

      EdgeIdx nEdgesCore = getCoreEdgeCount();

      L0File << nVerticesL0 << "  " << nEdgesCore << std::endl;
      for (VertexIdx i = 0; i < nVertices; i++)
      {
        if (!L0[i])
          continue;
        for (EdgeIdx j = offsets[i]; j < offsets[i + 1]; j++)
        {
          VertexIdx nbor = nbors[j];
          if (L0[nbor] && nbor > i)
            L0File << i << "  " << nbor << std::endl;
        }
      }
      L0File.close();
    }
  };

}

#endif