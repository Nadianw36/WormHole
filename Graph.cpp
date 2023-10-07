#include "Escape/Graph.h"
#include "Escape/JointSort.h"
#include <algorithm>
#include <string>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

using namespace Escape;

//Basic binary search procedure
// Input: pointer array, index of last entry end, and val to search for
// Output: index if val is found, and -1 otherwise
VertexIdx Escape::binarySearch(EdgeIdx* array, VertexIdx end, EdgeIdx val)
{
    VertexIdx low = 0;
    VertexIdx high = end-1;
    VertexIdx mid;

    while (low <= high)
    {
        mid = (low+high)/2;
        if (array[mid] == val)
            return mid;
        if (array[mid] > val)
            high = mid-1;
        if (array[mid] < val)
            low = mid+1;
    }
    return -1;
}

// comparator that only compares the first in pair
bool Escape::pairCompareFirst(Pair firstPair, Pair nextPair)
{
    return firstPair.first < nextPair.first;
}

// comparator that compares the second in pair. If they are equal, then compare first in pair
bool Escape::pairCompareSecond(Pair firstPair, Pair nextPair)
{
    if (firstPair.second != nextPair.second)
        return firstPair.second < nextPair.second;
    return firstPair.first < nextPair.first;
}


Graph Escape::newGraph(VertexIdx nVertices, EdgeIdx nEdges)
{
  return {nVertices, nEdges, new VertexIdx[nEdges], new VertexIdx[nEdges]};
}


Graph Graph::copy() const
{
  Graph ret = newGraph(nVertices, nEdges);
  std::copy(srcs, srcs + nEdges, ret.srcs);
  std::copy(dsts, dsts + nEdges, ret.dsts);
  return ret;
}


void Graph::print(FILE* f) const
{
  fprintf(f, "nVertices=%ld nEdges=%ld\n", (int64_t) nVertices, (int64_t) nEdges);
  for (EdgeIdx i = 0; i < nEdges; ++i)
    fprintf(f, "  %ld   %ld -> %ld\n", (int64_t) i, (int64_t) srcs[i], (int64_t) dsts[i]);
}


void Escape::delGraph(Graph g)
{
  delete[] g.srcs;
  delete[] g.dsts;
}


CGraph Escape::newCGraph(VertexIdx nVertices, EdgeIdx nEdges)
{
  return {nVertices, nEdges, new EdgeIdx[nVertices + 1], new VertexIdx[nEdges]};
}


CGraph CGraph::copy() const
{
  CGraph ret = newCGraph(nVertices, nEdges);
  std::copy(offsets, offsets + nVertices + 1, ret.offsets);
  std::copy(nbors, nbors + nEdges, ret.nbors);
  return ret;
}


void CGraph::print(FILE *f) const
{
  fprintf(f, "nVertices=%ld nEdges=%ld\n", (int64_t) nVertices, (int64_t) nEdges);
  for (VertexIdx i = 0; i < nVertices; ++i)
  {
    fprintf(f, "%ld: ", (int64_t) i);
    for (EdgeIdx j = offsets[i]; j < offsets[i + 1]; ++j)
      fprintf(f, "%ld ", (int64_t) nbors[j]);
    fprintf(f, "\n");
  }
}


void Escape::delCGraph(CGraph g)
{
  delete[] g.offsets;
  delete[] g.nbors;
}


// Checks if edge (v1, v2) is present in CGraph
// If edge is present: Return index (in nbors) of v2 in neighbors of v1
// If edge is not present: Return -1
int CGraph::isEdge(VertexIdx v1, VertexIdx v2)
{
    if (v1 >= nVertices)
        return -1;
    for (EdgeIdx i=offsets[v1]; i < offsets[v1+1]; ++i)
        if (nbors[i] == v2)
            return i;
    return -1;
}

//Checks if edge (v1, v2) is present using binary search
//Assumes CGraph is sorted by ID
//If edge is present: Return index (in nbors) of v2 in neighbors of v2
//If edge is not present: return -1
EdgeIdx CGraph::getEdgeBinary(VertexIdx v1, VertexIdx v2) const
{
    if (v1 >= nVertices)
        return -1;
    EdgeIdx low = offsets[v1];
    EdgeIdx high = offsets[v1+1]-1;
    EdgeIdx mid;

    while(low <= high)
    {
        mid = (low+high)/2;
        if (nbors[mid] == v2)
            return mid;
        if (nbors[mid] > v2)
            high = mid-1;
        else
            low = mid+1;
    }
    return -1;
}

//Checks if edge (v1, v2) is present using binary search
//Assumes CGraph is sorted by ID
//If edge is present: return true
//If edge is not present: return false
bool CGraph::isEdgeBinary(VertexIdx v1, VertexIdx v2) const
{
//     VertexIdx deg1 = offsets[v1+1] - offsets[v1];
//     VertexIdx deg2 = offsets[v2+1] - offsets[v2];

//     if(deg2 < deg1)
//     {
//         VertexIdx swp = v1;
//         v1 = v2;
//         v2 = swp;
//     }
//
    EdgeIdx low = offsets[v1];
    EdgeIdx high = offsets[v1+1]-1;
    EdgeIdx mid;

    while(low <= high)
    {
        mid = (low+high)/2;

        if (nbors[mid] == v2)
            return true;
        if (nbors[mid] > v2)
            high = mid-1;
        else
            low = mid+1;
    }
    return false;
}

//Same as above, changing name since we retrieve the edge index as well.
//There is a good reason to have both isEdge and getEdge, since it's easier
//to answer the first question.  Recommend changing isEdge to return bool - vv
EdgeIdx CGraph::getEdge(VertexIdx v1, VertexIdx v2) const
{
    if (v1 >= nVertices)
        return -1;
    for (EdgeIdx i=offsets[v1]; i < offsets[v1+1]; ++i)
        if (nbors[i] == v2)
            return i;
    return -1;
}


// This sorts each individual adjacency list by vertex ID. This is useful for
// doing a binary search, or for merging neighbor lists to find common neighbors.
void CGraph::sortById() const
{
    for (VertexIdx i=0; i < nVertices; i++)
        std::sort(nbors+offsets[i],nbors+offsets[i+1]);
}

// This outputs a new, isomorphic CGraph where vertex labels are in increasing order corresponding to degree.
// Thus, (after the relabeling), for all i < j, the degree of i is less than that of j.

CGraph CGraph::renameByDegreeOrder() const
{
    CGraph ret = newCGraph(nVertices, nEdges);
    Pair *deg_info = new Pair[nVertices];

    VertexIdx *mapping = new VertexIdx[nVertices];
    VertexIdx *inverse = new VertexIdx[nVertices];


    // Construct array of pairs, storing old vertex label and degree
    for (VertexIdx i=0; i < nVertices; i++)
    {
        deg_info[i].first = i;
        deg_info[i].second = offsets[i+1]-offsets[i];
    }

    // sort the pairs by degree (if degree is same, sort by old vertex label)
    std::sort(deg_info,deg_info+nVertices,Escape::pairCompareSecond);

    // Construct the mapping of old vertex label to new vertex label
    // So mapping[i] is what i is mapped to
    // And inverse[i] is what maps to i
    for (VertexIdx i=0; i < nVertices; i++)
    {
        mapping[deg_info[i].first] = i;
        inverse[i] = deg_info[i].first;
    }

    // Initialize offsets of output CGraph
    ret.offsets[0] = 0;
    EdgeIdx current = 0;


    // Loop over new vertices
    for (VertexIdx new_label=0; new_label < nVertices; new_label++)
    {
        VertexIdx old_label = inverse[new_label]; // corresponding old label for new vertices
        for (EdgeIdx pos = offsets[old_label]; pos < offsets[old_label+1]; pos++) // loop over neighbors of old label
        {
            VertexIdx old_nbr = nbors[pos];
            VertexIdx new_nbr = mapping[old_nbr]; //corresponding new neighbor
            ret.nbors[current] = new_nbr; // insert new neighbor in nbors of output
            current++;
        }
        ret.offsets[new_label+1] = current; // all neighbors of new_label have been added, so we set offset for new_label+1
    }

    return ret;
}


CGraph Escape::makeCSR(Graph g, bool inPlace)
{
  //create a temporary graph for sorting unless the user requests in-place
  //operation.
  auto tmpG = inPlace ? g : g.copy();

  auto begin = JSIterator<VertexIdx, VertexIdx>{tmpG.srcs, tmpG.dsts};
  auto end   = begin + tmpG.nEdges;
  std::sort(begin, end);

  //We steal the dsts array from tmpG.
  CGraph ret = {g.nVertices, g.nEdges, new EdgeIdx[g.nVertices + 1], tmpG.dsts};

  //Now we have everything sorted by src, compress:
  VertexIdx cv = 0;
  for (EdgeIdx i = 0; i < tmpG.nEdges; ++i)
  {
    auto src = tmpG.srcs[i];
    while (cv <= src)
      ret.offsets[cv++] = i;
  }
  while (cv <= g.nVertices)
    ret.offsets[cv++] = g.nEdges;

  delete[] tmpG.srcs; //we retain tmpG.dsts in the output

  ret.sortById();

  return ret;
}

VertexIdx Escape::find_index(VertexIdx *arr, VertexIdx nVertices, VertexIdx value){
    for(VertexIdx i = 0; i<nVertices; ++i){
        if(arr[i]==value){
            return i;
        }
    }
    return -1;
}

CGraph Escape::makeCSC(Graph g, bool inPlace)
{
  return makeCSR({g.nVertices, g.nEdges, g.dsts, g.srcs}, inPlace);
}
void CGraph::saveGraphToFile(std::string graph_name) const{
    std::ofstream outFile("bin/"+graph_name+"_CGraph.bin", std::ios::binary);

    // Write the number of vertices and edges
    outFile.write(reinterpret_cast<const char*>(&nVertices), sizeof(VertexIdx));
    outFile.write(reinterpret_cast<const char*>(&nEdges), sizeof(EdgeIdx));

    // Write the offsets array
    outFile.write(reinterpret_cast<const char*>(offsets), sizeof(EdgeIdx) * (nVertices + 1));

    // Write the neighbors array
    outFile.write(reinterpret_cast<const char*>(nbors), sizeof(VertexIdx) * offsets[nVertices]);

    outFile.close();
}
void CGraph::loadGraphFromFile(std::string graph_name) {
    std::ifstream inFile("bin/"+graph_name+"_CGraph.bin", std::ios::binary);

    // Read the number of vertices and edges
    inFile.read(reinterpret_cast<char*>(&nVertices), sizeof(VertexIdx));
    inFile.read(reinterpret_cast<char*>(&nEdges), sizeof(EdgeIdx));

    // Allocate memory for offsets array and read the data
    offsets = new EdgeIdx[nVertices + 1];
    inFile.read(reinterpret_cast<char*>(offsets), sizeof(EdgeIdx) * (nVertices + 1));

    // Allocate memory for neighbors array and read the data
    nbors = new VertexIdx[offsets[nVertices]];
    inFile.read(reinterpret_cast<char*>(nbors), sizeof(VertexIdx) * offsets[nVertices]);

    inFile.close();
}

  void CGraph::sanitizeForL0(std::string filename, bool* L0, bool* L1, EdgeIdx new_nEdges){
    EdgeIdx counter = 0;
    printf("start iterating\n");
    for (VertexIdx i = 0; i < nVertices; i++){
      EdgeIdx new_offset = counter;
      for(EdgeIdx e = offsets[i]; e < offsets[i+1]; e++){
        VertexIdx nbor = nbors[e];
        //if(counter>14100000) printf("through vertex %ld", counter);
        if(L0[i] && L1[nbor]) {
            // printf("L0 -> L1: %ld  %ld \n", i, nbor);
            continue;
        }
        else if (L1[i] && L1[nbor]) continue;
        else if (L1[i] && !L0[nbor] && !L1[nbor]) continue;
        else{
            nbors[counter] = nbor;
            counter++;
        }
      }
      offsets[i] = new_offset;
    }
    printf("done iterating %ld\n", counter);
    offsets[nVertices] = counter;
    nEdges = counter;
    printf("WHY AM I BEING AN ASSHOLE %ld\n", nEdges);
  }
  void CGraph::createSanitizedFile(std::string filename){
    std::ofstream file(filename);
    file<<std::endl;
    EdgeIdx sEdges = 0;
    for(VertexIdx v = 0; v<nVertices; v++){
      std::unordered_set<VertexIdx> nbor_set;
      for(EdgeIdx e = offsets[v]; e<offsets[v+1]; e++){
        VertexIdx nbor = nbors[e];
        if(nbor_set.find(nbor) == nbor_set.end() && nbor > v){
            nbor_set.insert(nbor);
            file  << v << " " << nbor << "\n";
            sEdges ++;
        }
      }
    }
    std::cout << nVertices << " " << sEdges << std::endl;

    file.seekp(0, std::ios::beg);
    file << nVertices << " " << sEdges;
    file.close();
  }
// Returns the distance between p1 and p2
int CGraph::BFS(VertexIdx p1, VertexIdx p2) const
{
    // Mark all the vertices as not visited
    bool* visited = new bool[nVertices];
    std::fill_n(visited, nVertices, false);
 
    // Create a queue for BFS
    std::queue<VertexIdx> q;
    bool found = false;
 
    // Mark the current node as visited and enqueue it
    visited[p1] = true;
    q.push(p1);
    std::unordered_map<VertexIdx, VertexIdx> parent; // map to retrace BFS (node -> parent)
    
    while(!q.empty() && !found)
    {
        VertexIdx s = q.front();
        q.pop();
        // Queue unexplored neighbors
        for (EdgeIdx j = offsets[s]; j < offsets[s+1]; j++){
            VertexIdx nbor = nbors[j];
            if(!visited[nbor]){
                visited[nbor] = true;
                q.push(nbor);
                parent[nbor] = s;
            }
            if (nbor==p2){
                found = true;
                break;
            }
        }
    }
    if (!found) return -1;

    // Rebuild shortest path to find distance
    int distance = 0;
    VertexIdx current = p2;
    while (current!=p1){
        current = parent[current];
        distance++;
    }
    return distance;
}
SamplingInformation CGraph::BidirectionalBFS(VertexIdx p1, VertexIdx p2) const
{
    bool* visited = new bool[nVertices];
    std::fill_n(visited, nVertices, false);

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
    
    int mu = 2*nVertices;


    while(!q1.empty() && !q2.empty()){
        VertexIdx s1 = q1.front();
        q1.pop();
        VertexIdx s2 = q2.front();
        q2.pop();
        if(distance1[s1]+distance2[s2]>mu){
            return SamplingInformation(visited, mu, nVertices);
        }
        for (EdgeIdx j = offsets[s1]; j < offsets[s1+1]; j++){
            VertexIdx nbor = nbors[j];
            if(!visited[nbor]){
                visited[nbor] = true;
                q1.push(nbor);
                distance1[nbor] = distance1[s1]+1;
            } else if(distance1.find(nbor) == distance1.end()){
                distance1[nbor] = distance1[s1]+1;
                q1.push(nbor);
                if(distance1[nbor] + distance2[nbor]< mu){
                    mu = distance1[nbor] + distance2[nbor];
                }
            }
        }
        for (EdgeIdx j = offsets[s2]; j < offsets[s2+1]; j++){
            VertexIdx nbor = nbors[j];
            if(!visited[nbor]){
                visited[nbor] = true;
                q2.push(nbor);
                distance2[nbor] = distance2[s2] + 1;
            } else if(distance2.find(nbor) == distance2.end()){
                distance2[nbor] = distance2[s2]+1;
                q2.push(nbor);
                if(distance1[nbor] + distance2[nbor]< mu){
                    mu = distance1[nbor] + distance2[nbor];
                }
            } 
        }
    }
    if(mu != 2*nVertices) return SamplingInformation(visited, mu, nVertices);

    return SamplingInformation(visited, nVertices);
}