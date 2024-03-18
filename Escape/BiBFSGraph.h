#ifndef ESCAPE_BIBFSGRAPH_H_
#define ESCAPE_BIBFSGRAPH_H_

#include <cstdlib>
#include <cstdio>
#include <set>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include "Escape/Graph.h"
#include "Escape/Config.h"

using namespace std;

namespace Escape
{

    // class SamplingInformation
    // {
    // public:
    //     bool *visited;
    //     int estimate_distance = -1;
    //     VertexIdx nVertices;
    //     double p1_time_to_L0 = -1;
    //     double p2_time_to_L0 = -1;
    //     EdgeIdx p1_dist_to_L0 = -1;
    //     EdgeIdx p2_dist_to_L0 = -1;
    //     EdgeIdx p1_L0_hops = -1;
    //     EdgeIdx p2_L0_hops = -1;

    //     SamplingInformation() {}
    //     // init for when a path couldn't be found
    //     SamplingInformation(bool *verticesTouched, VertexIdx nv)
    //     {
    //         visited = verticesTouched;
    //         nVertices = nv;
    //     }
    //     // init for when a path is found
    //     SamplingInformation(bool *verticesTouched, int distance, VertexIdx nv)
    //     {
    //         estimate_distance = distance;
    //         visited = verticesTouched;
    //         nVertices = nv;
    //     }
    //     void insert(bool *new_visited)
    //     {
    //         std::transform(visited, visited + nVertices, new_visited, visited, [](bool &a, bool b)
    //                        { return a |= b; });
    //         delete[] new_visited;
    //     }
    //     int countQueries()
    //     {
    //         return std::count(visited, visited + (int64_t)nVertices, true);
    //     }
    // };

    struct BiBFSGraph : CGraph
    {
        // BiBFS stuff
        bool *visited;
        VertexIdx *queue1; // first two slots are reserved for start index and end index
        VertexIdx *queue2; // first two slots are reserved for start index and end index

        const VertexIdx Q_START_IDX = 0; // index for the queue start pointer
        const VertexIdx Q_END_IDX = 1;   // index for the queue end pointer
        const VertexIdx Q_START = 2;     // actual start of the queue

        VertexIdx *level1Record; // first slot reserved for current level
        VertexIdx *level2Record; // first slot reserved for current level

        const VertexIdx LEVEL_IDX = 0;
        const VertexIdx LEVEL_START = 1;

        int8_t *distance;
        bool collision;
        const int INFI = 126;

        BiBFSGraph() {}

        BiBFSGraph(const CGraph &graph) : CGraph(graph)
        {

            // printf("%ld %ld \n", nVertices, nEdges);
            queue1 = new VertexIdx[nVertices + 2];
            queue2 = new VertexIdx[nVertices + 2];
            distance = new int8_t[nVertices];

            memset(distance, INFI, nVertices);

            level1Record = new VertexIdx[INFI + 1];
            level2Record = new VertexIdx[INFI + 1];

            visited = new bool[nVertices];
            memset(visited, false, nVertices);

            collision = false;
        };
        bool *getVisited() const { return visited; }

        void reset()
        {
            memset(distance, INFI, nVertices);
            memset(visited, false, nVertices);
            collision = false;
        }

        void queueVertex(VertexIdx v, VertexIdx *queue, int d)
        {
            visited[v] = true;
            distance[v] = d;
            queue[queue[Q_END_IDX]++] = v;
        }

        bool Q1empty() const { return queue1[Q_START_IDX] == queue1[Q_END_IDX]; }
        bool Q2empty() const { return queue2[Q_START_IDX] == queue2[Q_END_IDX]; }

        bool virtual exploreVertexQ1(VertexIdx nbor, VertexIdx parent) { return true; }
        bool virtual exploreVertexQ2(VertexIdx nbor, VertexIdx parent) { return true; }

        bool virtual continueSearch() { return !collision && !Q1empty() && !Q2empty(); }

        int bidirectionalBFS();

        int BidirectionalBFS(VertexIdx p1, VertexIdx p2);
    };

}
#endif
