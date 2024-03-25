#include "Escape/BiBFSGraph.h"
#include <algorithm>
#include <string>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

using namespace Escape;

int BiBFSGraph::bidirectionalBFS()
{
    int8_t pathDistance = INFI;

    VertexIdx v, nbor;

    bool exploreQ1Level;

    // printf("queue indices %ld %ld %ld %ld\n", queue1[Q_START_IDX], queue1[Q_END_IDX], queue2[Q_START_IDX], queue2[Q_END_IDX]);

    while (continueSearch())
    {
        // printf("queue indices %ld %ld %ld %ld\n", queue1[Q_START_IDX], queue1[Q_END_IDX], queue2[Q_START_IDX], queue2[Q_END_IDX]);
        exploreQ1Level = (level1Record[LEVEL_IDX] <= level2Record[LEVEL_IDX] && !Q1empty()) || Q2empty();
        if (exploreQ1Level)
        {
            // printf("1: starting level %ld \n", level1Record[LEVEL_IDX]);
            while (queue1[Q_START_IDX] < level1Record[level1Record[LEVEL_IDX]])
            {
                v = queue1[queue1[Q_START_IDX]++];
                // printf("1: checking vertex %ld start now at %ld \n", v, queue1[Q_START_IDX]);
                for (EdgeIdx j = offsets[v]; j < offsets[v + 1]; j++)
                {
                    nbor = nbors[j];
                    if (exploreVertexQ1(nbor, v))
                    {
                        // printf("1: checking neighbor %ld \n", nbor);
                        if (!visited[nbor])
                        {
                            // printf("1: queueing %ld\n", nbor);
                            queueVertex(nbor, queue1, distance[v] + 1);
                            // printf("1: done queueing \n");
                        }
                        else if (distance[nbor] < 0) // v has been visited from side 2
                        {
                            // printf("1: resolving merge %ld\n", nbor);
                            int temp_distance = (-1 * distance[nbor]) + distance[v] - 1;
                            // printf("%d %d %d %d\n", temp_distance, pathDistance, distance[nbor], distance[v]);
                            if (temp_distance < pathDistance)
                                pathDistance = temp_distance;
                            collision = true;
                            // printf("1: done resolving merge \n");
                        }
                    }
                }
            }
            // printf("1: updating record\n");
            level1Record[++level1Record[LEVEL_IDX]] = queue1[Q_END_IDX];
            // printf("1: finished updating record %ld %ld %ld\n", level1Record[LEVEL_IDX], queue1[Q_END_IDX], level1Record[level1Record[LEVEL_IDX]]);
        }
        else
        {
            // printf("2: starting level %ld \n", level2Record[LEVEL_IDX]);
            while (queue2[Q_START_IDX] < level2Record[level2Record[LEVEL_IDX]])
            {
                v = queue2[queue2[Q_START_IDX]++];
                // printf("2: checking vertex %ld start now at %ld \n", v, queue2[Q_START_IDX]);
                for (EdgeIdx j = offsets[v]; j < offsets[v + 1]; j++)
                {
                    nbor = nbors[j];
                    if (exploreVertexQ2(nbor, v))
                    {
                        // printf("2:checking neighbor %ld \n", nbor);
                        if (!visited[nbor])
                        {
                            // printf("2: queueing %ld\n", nbor);
                            queueVertex(nbor, queue2, distance[v] - 1);
                            // printf("2: done queueing \n");
                        }
                        else if (distance[nbor] > 0) // v has been visited from side 1
                        {
                            // printf("2: resolving merge %ld\n", nbor);
                            int temp_distance = distance[nbor] + (-1 * distance[v]) - 1;
                            // printf("%d %d %d %d\n", temp_distance, pathDistance, distance[nbor], distance[v]);
                            if (temp_distance < pathDistance)
                                pathDistance = temp_distance;
                            collision = true;
                            // printf("2: done resolving merge \n");
                        }
                    }
                }
            }
            // printf("2: updating record\n");
            level2Record[++level2Record[LEVEL_IDX]] = queue2[Q_END_IDX];
            // printf("2: finished updating record %ld %ld %ld\n", level2Record[LEVEL_IDX], queue2[Q_END_IDX], level2Record[level2Record[LEVEL_IDX]]);
        }
    }
    if (pathDistance == INFI)
        return -1;
    return pathDistance;
};

int BiBFSGraph::BidirectionalBFS(VertexIdx p1, VertexIdx p2)
{
    // printf("resetting \n");
    reset();
    // printf("setting distance \n");

    // printf("setting q1 \n");
    queue1[Q_START_IDX] = Q_START;
    queue1[Q_END_IDX] = Q_START;
    queueVertex(p1, queue1, 1);
    // printf("finished setting q1 %ld %ld %ld\n", queue1[Q_START_IDX], queue1[Q_END_IDX], queue1[Q_START]);

    // printf("setting q2 \n");
    queue2[Q_START_IDX] = Q_START;
    queue2[Q_END_IDX] = Q_START;
    queueVertex(p2, queue2, -1);
    // printf("finished setting q2 %ld %ld %ld\n", queue2[Q_START_IDX], queue2[Q_END_IDX], queue2[Q_START]);

    level1Record[LEVEL_START] = queue1[Q_END_IDX];
    level1Record[LEVEL_IDX] = LEVEL_START;

    // printf("finished level record 1 %ld %ld %ld\n", level1Record[LEVEL_IDX], level1Record[LEVEL_START]);

    level2Record[LEVEL_START] = queue2[Q_END_IDX];
    level2Record[LEVEL_IDX] = LEVEL_START;

    // printf("finished level record 2 %ld %ld %ld\n", level2Record[LEVEL_IDX], level2Record[LEVEL_START]);

    return bidirectionalBFS();
};