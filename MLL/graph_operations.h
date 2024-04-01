// Created by 24148 on 6/20/2021.

#ifndef PRUNEDTD_GRAPH_OPERATIONS_H
#define PRUNEDTD_GRAPH_OPERATIONS_H

#include "header.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
namespace graph_ops
{
  void convert_to_undirected(std::vector<std::vector<unsigned int>> &g)
  {
    for (unsigned int src = 0; src < g.size(); src++)
    {
      for (unsigned int i = 0; i < g[src].size(); i++)
      {
        auto &dst = g[src][i];
        if (src != dst)
        {
          g[dst].emplace_back(src);
        }
      }
    }

    for (unsigned int src = 0; src < g.size(); src++)
    {
      // remove self loops
      for (unsigned int i = 0; i < g[src].size(); i++)
      {
        if (src == g[src][i])
        {
          g[src][i] = g[src].back();
          g[src].pop_back();
          i--;
        }
      }
      // remove duplicated edges
      auto &v = g[src];

      sort(v.begin(), v.end());
      auto it = unique(v.begin(), v.end());
      v.erase(it, v.end());
    }
  }

  void convert_to_single_directed_graph(std::vector<std::vector<unsigned int>> &g)
  {
    std::vector<std::vector<unsigned int>> n_g;
    n_g.resize(g.size());

    for (int i = 0; i < g.size(); ++i)
    {
      for (int j = 0; j < g[i].size(); ++j)
      {
        if (g[i][j] > i)
        {
          n_g[i].push_back(g[i][j]);
        }
      }
    }
    swap(g, n_g);
  }

  std::vector<std::vector<unsigned int>>
  convert2unweighted_graph(std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &&g)
  {
    vector<vector<unsigned int>> graph(g.size());
    for (int i = 0; i < g.size(); ++i)
    {
      for (int j = 0; j < g[i].size(); ++j)
      {
        graph[i].push_back(g[i][j].first);
      }
    }
    g.clear();
    return graph;
  }

  /*
   * returned id_map, <0 equ1, >=n equ2
   */
  std::vector<int> equivalance_reduction(std::vector<std::vector<unsigned int>> &g)
  {

    unsigned long m = 0;
    unsigned int n = g.size();
    for (const auto &v : g)
    {
      m += v.size();
    }
    vector<int> id_map(n);
    for (int i = 0; i < id_map.size(); ++i)
    {
      id_map[i] = i;
    }

    vector<long> f1(n, 0);
    vector<long> newid(m + n + 1, 0);
    long inc_nid = 0, previous_max_nid;
    for (uint32_t u = 0; u < g.size(); u++)
    {
      previous_max_nid = inc_nid;
      for (auto v : g[u])
      {
        if (newid[f1[v]] <= previous_max_nid)
        { // old group id, and separate out a new group
          newid[f1[v]] = ++inc_nid;
        }
        f1[v] = newid[f1[v]];
      }
    }
    fill(newid.begin(), newid.end(), n);
    for (uint32_t u = 0; u < f1.size(); ++u)
    {
      if (newid[f1[u]] == n)
      {
        newid[f1[u]] = u;
      }
      else
      {
        id_map[u] = -newid[f1[u]] - 1;
      }
      f1[u] = newid[f1[u]];
    }

    // cal f2
    vector<long> f2(n, 0);
    fill(newid.begin(), newid.end(), 0);
    inc_nid = 0;
    for (uint32_t u = 0; u < g.size(); u++)
    {
      previous_max_nid = inc_nid;

      newid[f2[u]] = ++inc_nid;
      f2[u] = newid[f2[u]];

      for (auto v : g[u])
      {
        if (newid[f2[v]] <= previous_max_nid)
        { // old group id, and separate out a new group
          newid[f2[v]] = ++inc_nid;
        }
        f2[v] = newid[f2[v]]; // separate v from old to new (assign a new group id)
      }
    }
    fill(newid.begin(), newid.end(), n);
    for (uint32_t u = 0; u < f2.size(); ++u)
    {
      if (newid[f2[u]] == n)
      {
        newid[f2[u]] = u;
      }
      else
      {
        id_map[u] = id_map.size() + newid[f2[u]];
      }
      f2[u] = newid[f2[u]];
    }

    // reduct graph
    inc_nid = 0;
    for (uint32_t u = 0; u < g.size(); ++u)
    {
      if (f1[u] != u || f2[u] != u)
      {
        g[u].clear();
        continue;
      }
      uint32_t p = 0;
      for (uint32_t j = 0; j < g[u].size(); j++)
      {
        int v = g[u][j];
        if (f1[v] == v && f2[v] == v)
          g[u][p++] = v;
      }
      g[u].resize(p);
      newid[u] = inc_nid++;
    }
    for (int i = 0; i < id_map.size(); ++i)
    {
      if (id_map[i] >= 0 && id_map[i] < id_map.size())
      {
        id_map[i] = newid[id_map[i]];
      }
    }
    // assign new id;
    vector<vector<unsigned int>> new_graph(inc_nid);
    for (uint32_t u = 0; u < g.size(); ++u)
    {
      if (f1[u] != u || f2[u] != u)
        continue;
      for (auto v : g[u])
      {
        new_graph[newid[u]].push_back(newid[v]);
      }
    }
    g = std::move(new_graph);

    return id_map;
  }

  /*
   * id->r,
   * if corresponding r is 0, means vid has max degree.
   */
  vector<unsigned int> order_vertexes(const vector<vector<unsigned int>> &g)
  {
    vector<pair<unsigned int, unsigned int>> r2degvid;
    r2degvid.reserve(g.size());
    for (int i = 0; i < g.size(); ++i)
    {
      r2degvid.emplace_back(g[i].size(), i);
    }
    sort(r2degvid.begin(), r2degvid.end(), greater<>());
    vector<unsigned int> id2r(g.size());
    for (int i = 0; i < r2degvid.size(); ++i)
    {
      id2r[r2degvid[i].second] = i;
    }
    return id2r;
  }

  void relabel_graphs(vector<vector<unsigned int>> &g, vector<unsigned int> &r)
  {
    vector<vector<unsigned int>> new_graph(g.size());
    for (int i = 0; i < g.size(); ++i)
    {
      for (auto neib : g[i])
      {
        new_graph[r[i]].push_back(r[neib]);
      }
    }
    g = std::move(new_graph);
    // remove zero degree vertices
    while (g.back().empty())
      g.pop_back();
  }

  void store_graph_bin(const string &path, const vector<vector<unsigned int>> &graph)
  {
    ofstream ofs(path);
    if (!ofs)
    {
      mlog("open file %s failed.", path.c_str());
      exit(-1);
    }

    unsigned int size = graph.size();
    ofs.write(reinterpret_cast<const char *>(&size), sizeof(size));
    for (const auto &i : graph)
    {
      size = i.size();
      ofs.write(reinterpret_cast<const char *>(&size), sizeof(size));
      ofs.write(reinterpret_cast<const char *>(i.data()), sizeof(unsigned int) * size);
    }
  }

  template <typename T>
  void store_vector_bin(const string &path, const vector<T> &arr)
  {
    FILE *fo = fopen(path.c_str(), "wb");
    if (fo == nullptr)
    {
      printf("write file %s failed.\n", path.c_str());
    }
    int size = arr.size();
    fwrite(&size, sizeof(size), 1, fo);
    fwrite(arr.data(), sizeof(T), arr.size(), fo);
    fclose(fo);
  }

  template <typename T>
  bool load_vector_bin(const string &path, vector<T> &arr)
  {
    FILE *fi = fopen(path.c_str(), "rb");
    if (fi == nullptr)
    {
      printf("read file %s failed.\n", path.c_str());
      return false;
    }
    int size;
    fread(&size, sizeof(size), 1, fi);
    arr.resize(size);
    fread(arr.data(), sizeof(T), arr.size(), fi);
    fclose(fi);
    return true;
  }

  bool load_graph_bin(const string &path, vector<vector<unsigned int>> &graph)
  {
    ifstream ifs(path);
    if (!ifs)
    {
      mlog("open file %s failed.", path.c_str());
      return false;
    }
    unsigned int size;
    ifs.read(reinterpret_cast<char *>(&size), sizeof(size));
    graph.resize(size);
    for (auto &i : graph)
    {
      ifs.read(reinterpret_cast<char *>(&size), sizeof(size));
      i.resize(size);
      ifs.read(reinterpret_cast<char *>(i.data()), sizeof(unsigned int) * size);
    }
    return true;
  }
} // namespace graph_ops

#endif // PRUNEDTD_GRAPH_OPERATIONS_H
