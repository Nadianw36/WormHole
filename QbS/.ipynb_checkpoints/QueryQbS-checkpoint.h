#ifndef QBS_QUERYQBS_H_
#define QBS_QUERYQBS_H_

#include <iostream>
#include <stdint.h>
#include <cstring>
#include <stdlib.h>
#include <vector>
#include <map>
#include <set>
#include <math.h>
#include <queue>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <sys/time.h>

#define INF 99

using namespace std;

namespace QbS
{

    struct pairLM
    {
        int lm1;
        int lm2;

        int to1;
        int to2;
        pairLM(int a, int b, int c, int d)
        {
            lm1 = a;
            lm2 = b;
            to1 = c;
            to2 = d;
        };
    };

    class spQuery
    {

    private:
        int countVertex;
        int diameter;
        int k;

        int *queueA, frontA, endA; // simulate a queue "queueA"
        int *queueB, frontB, endB; // simulate a queue "queueB"
        int *levelARecord, levelA; // simulate a queue "levelARecord"
        int *levelBRecord, levelB; // simulate a queue "levelBRecord"
        int *meet, endMeet;        // simulate a queue "meet"

        int8_t *dist; // a distance list

        vector<vector<int>> adjList;
        uint8_t **labels;
        uint8_t **highway;
        vector<vector<string>> highwayPath;

        int maxDistance;
        vector<pairLM> combination;

        vector<vector<int>> pairs;
        long totalTime;
        map<int, int> id;

        void loadMapFromFile(const std::string &filename)
        {

            std::ifstream inFile(filename, std::ios::binary);
            if (!inFile.is_open())
            {
                std::cerr << "Error opening file for reading" << std::endl;
                return; // Return an empty map on error
            }

            // Read the size of the map
            size_t mapSize;
            inFile.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));

            // Read each key-value pair and populate the map
            for (size_t i = 0; i < mapSize; ++i)
            {
                int key, value;
                inFile.read(reinterpret_cast<char *>(&key), sizeof(key));
                inFile.read(reinterpret_cast<char *>(&value), sizeof(value));
                id[key] = value;
            }

            inFile.close();
        }

        void initialization()
        {

            countVertex = adjList.size();

            queueA = new int[countVertex];
            queueB = new int[countVertex];

            levelARecord = new int[diameter];
            levelBRecord = new int[diameter];

            meet = new int[countVertex / 10];
            dist = new int8_t[countVertex];
            memset(dist, INF, countVertex);
        }

        void loadSparsifiedGraph(string graph_name, string graph_directory)
        {
            string graphFile = graph_directory + graph_name + "_spsf_" + to_string(k) + ".txt";
            int u, d;
            string line;
            stringstream ss;
            ifstream graph(graphFile);
            while (getline(graph, line))
            {
                ss << line;
                ss >> u;
                ss >> d;
                vector<int> temp(d, 0);
                for (int i = 0; i < d; i++)
                {
                    ss >> temp[i];
                }
                adjList.push_back(temp);
                ss.clear();
                ss.str("");
            }
            graph.close();
            loadMapFromFile(graph_directory + graph_name + "_id-map.bin");
        }

        void loadLabels(string graph_name, string graph_directory)
        {
            string graphFile = graph_directory + graph_name + "_spsf_" + to_string(k) + ".txt";
            string indexFile = graph_directory + graph_name + "_indexInfo_" + to_string(k) + ".txt";

            int u, v, d;
            labels = new uint8_t *[countVertex];
            for (int i = 0; i < countVertex; i++)
            {
                labels[i] = new uint8_t[k];
            }
            highway = new uint8_t *[k];
            vector<string> tempp(k, "");
            for (int i = 0; i < k; i++)
            {
                highway[i] = new uint8_t[k];
                highwayPath.push_back(tempp);
            }
            stringstream ss;
            string line;
            ifstream labelFile(indexFile);
            while (getline(labelFile, line))
            {
                if (line == "###")
                    break;
                ss << line;
                ss >> u;
                ss >> v;
                ss >> d;
                highway[u][v] = d;
                highway[v][u] = d;
                highwayPath[u][v] = " "; // ss.str();
                ss.clear();
                ss.str("");
            }
            while (getline(labelFile, line))
            {
                ss << line;
                ss >> u;
                for (int i = 0; i < k; i++)
                {
                    ss >> v;
                    // if(v==0) cout << ss.str() << endl;
                    labels[u][i] = v;
                }
                ss.clear();
                ss.str("");
            }
            cout << "|V|: " << countVertex << endl;
            labelFile.close();
        }

        unsigned int rand_32()
        {
            return (rand() & 0x3) << 30 | rand() << 15 | rand();
        }

    public:
        spQuery(int k, int diameter, string graph_name, string graph_directory)
        {
            spQuery::k = k;
            spQuery::diameter = diameter;
            loadSparsifiedGraph(graph_name, graph_directory);
            initialization();
            loadLabels(graph_name, graph_directory);
        }

        void generateTest(int testTime, string inputFileLoc)
        {
            ifstream inputFile;
            inputFile.open(inputFileLoc);
            vector<int> p(2);
            for (int i = 0; i < testTime; i++)
            {
                inputFile >> p[0] >> p[1];

                pairs.push_back(p);
            }
        }

        void seriesTest(int testTime)
        {
            timeval st, ed;
            gettimeofday(&st, NULL);
            for (int i = 0; i < testTime; i++)
            {
                queryEfficient(i);
            }
            gettimeofday(&ed, NULL);
            long sumTime = (ed.tv_sec - st.tv_sec) * 1000000 + (ed.tv_usec - st.tv_usec);
        }
        void queryEfficient(int test)
        {

            int a = pairs[test][0];
            int b = pairs[test][1];

            if (id.find(a) == id.end() || id.find(b) == id.end())
            {
                // printf("-1\n");
                return;
            }

            a = id[a];
            b = id[b];
            // // ADDED
            // printf("pair %d %d\n", a, b);
            // printf("neighbors of a \n");
            // for (int i = 0; i < adjList[a].size(); i++){
            //     // printf("%d ", adjList[a][i]);
            // }
            // // END ADDED
            // get the upper bound and the graph sketch based on the label information
            if (a == b)
            {
                // printf("0\n");
                return;
            }
            if (a < k && b < k)
            {
                // printf("%ld\n", highway[a][b]);
                return;
            }
            else if (a < k)
            {
                // printf("%ld\n", labels[b][a]);
                return;
            }
            else if (b < k)
            {
                // printf("%ld\n", labels[a][b]);
                return;
            }

            vector<pairLM>().swap(combination);
            maxDistance = 99;
            int maxASide = 0;
            int maxBSide = 0;
            int d;

            for (int i = 0; i < k; i++)
            {
                if (labels[a][i] == INF)
                    continue;

                // printf("label %d dist %d ", i, labels[a][i]);
                for (int j = 0; j < k; j++)
                {
                    if (labels[b][j] == INF)
                        continue;
                    d = labels[a][i] + labels[b][j] + highway[i][j];
                    if (d < maxDistance)
                    {
                        vector<pairLM>().swap(combination);
                        pairLM p(i, j, labels[a][i], labels[b][j]);
                        combination.push_back(p);
                        maxDistance = d;
                        maxASide = labels[a][i] - 1;
                        maxBSide = labels[b][j] - 1;
                    }
                    else if (d == maxDistance)
                    {
                        pairLM p(i, j, labels[a][i], labels[b][j]);
                        combination.push_back(p);
                        maxASide = max(maxASide, labels[a][i] - 1);
                        maxBSide = max(maxBSide, labels[b][j] - 1);
                    }
                }
            }
            // printf("\n");

            if (maxDistance >= diameter)
            {
                // printf("-1\n");
                return;
            }

            dist[a] = 1;  // from 1 to 2, 3, 4...
            dist[b] = -1; // from -1 to -2, -3, -4...

            frontA = 0;    // initialize the queueA by setting the front pointer to 0;
            frontB = 0;    // initialize the queueB by setting the front pointer to 0;
            queueA[0] = a; // push a into queueA;
            endA = 1;
            queueB[0] = b; // push b into queueB;
            endB = 1;
            levelARecord[0] = 1;
            levelBRecord[0] = 1;
            levelA = 0;
            levelB = 0;
            endMeet = 0;

            int u, v;

            // part A: bi-BFS search on the sparsifed graph
            int choice;
            int distance = maxDistance;
            while (endMeet == 0 && (levelA + levelB) < maxDistance) // the level-bounded bi-BFS (based on the upper bound)
            {
                choice = 1;
                if (levelB < maxDistance - maxASide)
                {
                    if (levelA < maxDistance - maxBSide && endA > endB)
                        choice = 0;
                    else if (levelA >= maxDistance - maxBSide)
                        choice = 0;
                }
                if (choice == 1)
                {
                    while (frontA < levelARecord[levelA])
                    {
                        u = queueA[frontA++];
                        for (int i = 0; i < adjList[u].size(); i++)
                        {
                            v = adjList[u][i];
                            if (dist[v] == INF) // v is unvisited before
                            {
                                dist[v] = dist[u] + 1;
                                queueA[endA++] = v; // push v into next level of queueA
                            }
                            else if (dist[v] < 0) // v has been visited from side b
                            {
                                meet[endMeet++] = v; // v is on the level where the bi-BFS meets
                                // printf("from side A %d %d\n", dist[v], dist[u]);
                                int temp_distance = (-1 * dist[v]) + dist[u] + 1;
                                if (temp_distance < distance)
                                    distance = temp_distance;
                                dist[v] = dist[u] + 1; // set dist[v]>0 for consistency
                                queueA[endA++] = v;
                            }
                        }
                    }
                    levelARecord[++levelA] = endA;
                }
                else // search from side b.
                {
                    while (frontB < levelBRecord[levelB])
                    {
                        u = queueB[frontB++];
                        for (int i = 0; i < adjList[u].size(); i++)
                        {
                            v = adjList[u][i];
                            if (dist[v] == INF) // v is unvisited before
                            {
                                dist[v] = dist[u] - 1;
                                queueB[endB++] = v;
                            }
                            else if (dist[v] > 0) // v has been visited from side a
                            {
                                meet[endMeet++] = v; // v is on the level where the bi-BFS meets
                                // printf("from side B %d %d\n", dist[v], dist[u]);
                                int temp_distance = dist[v] + (-1 * dist[u]) + 1;
                                if (temp_distance < distance)
                                    distance = temp_distance;
                                dist[v] = dist[u] - 1;
                                queueB[endB++] = v;
                            }
                        }
                    }
                    levelBRecord[++levelB] = endB;
                }
            }

            // printf("%d\n", distance);

            // part b: go back and find edges
            frontA = endA;
            frontB = endB;
            vector<vector<int>> sideA(maxDistance);
            vector<vector<int>> sideB(maxDistance);

            if (levelA + levelB == maxDistance)
            {
                for (int i = 0; i < combination.size(); i++)
                {
                    sideA[(combination[i].to1) - 1].push_back(combination[i].lm1);
                    sideB[(combination[i].to2) - 1].push_back(combination[i].lm2);
                }
            }

            if (levelA + levelB == maxDistance && maxASide > levelA) // if the steps we take are not enough to recover the shortest paths to landmarks
            {
                for (int i = 0; i < combination.size(); i++)
                {
                    if (combination[i].to1 - 1 > levelA)
                    {
                        queue<int> addi;
                        for (int j = levelARecord[levelA - 1]; j < levelARecord[levelA]; j++)
                        {
                            if (labels[queueA[j]][combination[i].lm1] == combination[i].to1 - levelA)
                            {
                                meet[endMeet++] = queueA[j];
                                addi.push(queueA[j]);
                            }
                        }
                        while (!addi.empty())
                        {
                            int qf = addi.front();
                            addi.pop();
                            if (labels[qf][combination[i].lm1] == 1)
                            {
                                // cout << qf << " " << combination.lm1 << " ";
                                continue;
                            }
                            for (int z = 0; z < adjList[qf].size(); z++)
                            {
                                int zp = adjList[qf][z];
                                if (labels[zp][combination[i].lm1] == labels[qf][combination[i].lm1] - 1)
                                {
                                    // cout << qf << " " << zp << " ";
                                    addi.push(zp);
                                }
                            }
                        }
                    }
                }
            }
            if (levelA + levelB == maxDistance && maxBSide > levelB)
            {
                for (int i = 0; i < combination.size(); i++)
                {
                    if (combination[i].to2 - 1 > levelB)
                    {
                        queue<int> addi;
                        for (int j = levelBRecord[levelB - 1]; j < levelBRecord[levelB]; j++)
                        {
                            if (labels[queueB[j]][combination[i].lm2] == combination[i].to2 - levelB)
                            {
                                meet[endMeet++] = queueB[j];
                                addi.push(queueB[j]);
                            }
                        }
                        while (!addi.empty())
                        {
                            int qf = addi.front();
                            addi.pop();
                            if (labels[qf][combination[i].lm2] == 1)
                            {
                                // cout << qf << " " << combination.lm1 << " ";
                                continue;
                            }
                            for (int z = 0; z < adjList[qf].size(); z++)
                            {
                                int zp = adjList[qf][z];
                                if (labels[zp][combination[i].lm2] == labels[qf][combination[i].lm2] - 1)
                                {
                                    // cout << qf << " " << zp << " ";
                                    addi.push(zp);
                                }
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < endMeet; i++) // push all predecessors of v w.r.t a and b to queueA and queueB separately
            {
                u = meet[i]; // dist[u] = (levelA+1) or (-1-levelB)
                bool foundPredA = false;
                bool foundPredB = false;
                for (int j = 0; j < adjList[u].size(); j++)
                {
                    v = adjList[u][j];
                    if (dist[v] == levelA && !foundPredA) // predecessor of v w.r.t a
                    {
                        // cout << u << " " << v << " ";
                        queueA[endA++] = v;
                        foundPredA = true;
                    }
                    else if (dist[v] == 0 - levelB && !foundPredB) // predecessors of v w.r.t b
                    {
                        // cout << u << " " << v << " ";
                        queueB[endB++] = v;
                        foundPredB = true;
                    }
                }
                dist[u] = INF;
            }
            // printf("%d %d\n", endA, endB);

            // printf("%d %d \n", endA, endB);

            if (!sideA[levelA].empty())
            {
                for (int i = levelARecord[levelA - 1]; i < levelARecord[levelA]; i++)
                {
                    u = queueA[i];
                    for (int j = 0; j < sideA[levelA].size(); j++)
                    {
                        if (labels[u][sideA[levelA][j]] == 1)
                        {
                            // cout << u << " " << sideA[levelA][j] << " " ;
                            for (int z = j + 1; z < sideA[levelA].size(); z++)
                            {
                                if (labels[u][sideA[levelA][z]] == 1)
                                {
                                    // cout << u << sideA[levelA][z] << " " ;
                                }
                            }
                            if (dist[u] == INF)
                                break;
                            for (int z = 0; z < adjList[u].size(); z++)
                            {
                                if (dist[adjList[u][z]] == levelA)
                                {
                                    queueA[endA++] = adjList[u][z];
                                }
                            }
                            break;
                        }
                    }
                    dist[u] = INF;
                }
            }
            else
            {
                for (int i = levelARecord[levelA - 1]; i < levelARecord[levelA]; i++)
                    dist[queueA[i]] = INF;
            }

            if (!sideB[levelB].empty())
            {
                for (int i = levelBRecord[levelB - 1]; i < levelBRecord[levelB]; i++)
                {
                    u = queueB[i];
                    for (int j = 0; j < sideB[levelB].size(); j++)
                    {
                        if (labels[u][sideB[levelB][j]] == 1)
                        {
                            // cout << u << " " << sideB[levelB][j] << " " ;
                            for (int z = j + 1; z < sideB[levelB].size(); z++)
                            {
                                if (labels[u][sideB[levelB][z]] == 1)
                                {
                                    // cout << u << sideB[levelB][z] << " " ;
                                }
                            }
                            if (dist[u] == INF)
                                break;
                            for (int z = 0; z < adjList[u].size(); z++)
                            {
                                if (dist[adjList[u][z]] == 0 - levelB)
                                {
                                    queueB[endB++] = adjList[u][z];
                                }
                            }
                            break;
                        }
                    }
                    dist[u] = INF;
                }
            }
            else
            {
                for (int i = levelBRecord[levelB - 1]; i < levelBRecord[levelB]; i++)
                    dist[queueB[i]] = INF;
            }
            levelA--;
            levelB--;

            while (levelA > 0)
            {
                int p = endA;
                while (frontA != p) // queueA not empty
                {
                    u = queueA[frontA++];
                    if (dist[u] == INF)
                        continue; // a vertex may appear in the queue more than once.
                    for (int i = 0; i < adjList[u].size(); i++)
                    {
                        v = adjList[u][i];
                        if (dist[v] + 1 == dist[u]) // predecessors w.r.t a
                        {
                            queueA[endA++] = v;
                            // cout << u << " " << v << " ";
                        }
                    }
                    dist[u] = INF; // avoid visiting a vertex twice by setting its dist=99
                }
                if (!sideA[levelA].empty())
                {
                    for (int i = levelARecord[levelA - 1]; i < levelARecord[levelA]; i++)
                    {
                        u = queueA[i];
                        for (int j = 0; j < sideA[levelA].size(); j++)
                        {
                            if (labels[u][sideA[levelA][j]] == 1)
                            {
                                // cout << u << " " << sideA[levelA][j] << " " ;
                                for (int z = j + 1; z < sideA[levelA].size(); z++)
                                {
                                    if (labels[u][sideA[levelA][z]] == 1)
                                    {
                                        // cout <<
                                    }
                                }
                                if (dist[u] == INF)
                                    break;
                                for (int z = 0; z < adjList[u].size(); z++)
                                {
                                    if (dist[adjList[u][z]] == levelA)
                                    {
                                        queueA[endA++] = adjList[u][z];
                                    }
                                }
                                break;
                            }
                        }
                        dist[u] = INF;
                    }
                }
                else
                {
                    for (int i = levelARecord[levelA - 1]; i < levelARecord[levelA]; i++)
                        dist[queueA[i]] = INF;
                }
                levelA--;
            }

            while (levelB > 0)
            {
                int p = endB;
                while (frontB != p)
                {
                    u = queueB[frontB++];
                    if (dist[u] == INF)
                        continue;
                    for (int i = 0; i < adjList[u].size(); i++)
                    {
                        v = adjList[u][i];
                        if (dist[v] == dist[u] + 1)
                        {
                            queueB[endB++] = v;
                            // cout << u << " " << v << " ";
                        }
                    }
                    dist[u] = INF;
                }
                if (!sideB[levelB].empty())
                {
                    for (int i = levelBRecord[levelB - 1]; i < levelBRecord[levelB]; i++)
                    {
                        u = queueB[i];
                        for (int j = 0; j < sideB[levelB].size(); j++)
                        {
                            if (labels[u][sideB[levelB][j]] == 1)
                            {
                                // cout << u << " " << sideB[levelB][j] << " " ;
                                for (int z = j + 1; z < sideB[levelB].size(); z++)
                                {
                                    if (labels[u][sideB[levelB][z]] == 1)
                                    {
                                        // cout <<
                                    }
                                }
                                if (dist[u] == INF)
                                    break;
                                for (int z = 0; z < adjList[u].size(); z++)
                                {
                                    if (dist[adjList[u][z]] == 0 - levelB)
                                    {
                                        queueB[endB++] = adjList[u][z];
                                    }
                                }
                                break;
                            }
                        }
                        dist[u] = INF;
                    }
                }
                else
                {
                    for (int i = levelBRecord[levelB - 1]; i < levelBRecord[levelB]; i++)
                        dist[queueB[i]] = INF;
                }
                levelB--;
            }

            dist[a] = INF;
            dist[b] = INF;

            return;
        }
    };
}
#endif