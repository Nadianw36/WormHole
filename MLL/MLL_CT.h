
#ifndef MLL_CT_H_
#define MLL_CT_H_

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <functional>
#include <iostream>
#include <malloc.h>
#include <map>
#include <omp.h>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "header.h"
#include "file_io.h"
#include "graph_operations.h"
#include "obj2str.h"

typedef char dint;

using namespace std::chrono;

// template<>
std::string str(const dint &obj)
{
  std::ostringstream oss;
  oss << int(obj);
  return oss.str();
}

struct
{
  string sub_cmd; //
  string dataset{"lctd_example"};
  string index_folder;   // folder for storing generated indexes.
  int max_w{0};          // max bag size in tree reduction
  int num_threads{1};    //
  int query_cnt{1000};   //
  int report_level{2};   //
  bool clear{false};     //
  int n_roots{4};        //
  int num_qgroups{5};    //
  int qgroup_size{1000}; //
  int sub_pct{100};      // sub_pct% of original graph, used for scalability.
  string query_path;
} params;

struct vd_bitseg
{
  static unsigned int dbits;
  static unsigned int MASK;

  static void init_static(unsigned int vbits)
  {
    dbits = 32 - vbits;
    MASK = (1u << dbits) - 1u;
  }

  static void init_static_byncore(unsigned int ncore)
  {
    ncore--;
    int vbits = 0;
    for (int i = 0; i < 32; ++i)
    {
      if ((1u << i) & ncore)
      {
        vbits = i + 1;
      }
    }
    init_static(vbits);
  }

  unsigned int value{0};

  vd_bitseg() = default;

  vd_bitseg(unsigned int vid, unsigned int d)
      : value(vid << dbits | min(d, 31u))
  {
    //        mlogn("c({},{})", vid, d)
  }

  unsigned int v() const { return value >> (dbits); }

  unsigned int d() const { return value & MASK; }
};

unsigned int vd_bitseg::dbits;

unsigned int vd_bitseg::MASK;

std::string str(const vd_bitseg &obj)
{
  ostringstream oss;
  oss << "(" << obj.v() << "," << obj.d() << ")";
  return oss.str();
  //  return fmt::format("({},{})", obj.v(), obj.d());
}

using namespace std;

const unsigned int MAXT = 65000;

const unsigned int MAXD = 120; // max possible distance in experiment (small world graphs, less than 32)

// select N_ROOTS root vertices together with its neighbors to build BP labels.
// #define N_ROOTS 4
const int MAX_N_ROOTS = 4;

const int MAX_BP_THREADS = 8;

typedef unsigned short tint;

template <typename T>
int locate_first1_in_bitmap(T bitmap)
{
  for (int i = 0; i < sizeof(bitmap) * 8; ++i)
  {
    if (1lu << i & bitmap)
    {
      return i;
    }
  }
  return -1;
}

/*
 * 0: report nothing
 * 1: basic report
 * 2: progressive report
 * 3: output_all
 * 4: develop&debug
 */

enum report_level
{
  nothing = 0,
  basic = 1,
  progress = 2,
  alldata = 3,
  debug = 4
};

/*
struct vdv_bitseg64 {
    static unsigned int dbits;
    static unsigned int vbits;
    static unsigned long dMASK;
    static unsigned long vMASK;

    static void init_static_byn(unsigned int n) {
        n--;
        vbits = 0;
        for (int i = 0; i < 64; ++i) {
            if ((1lu << i) & n) {
                vbits = i + 1;
            }
        }
        dbits = 64 - 2 * vbits;
        dMASK = (1lu << dbits) - 1lu;
        vMASK = (1lu << vbits) - 1lu;
    }

    unsigned long value{0};

    vdv_bitseg64() = default;

    vdv_bitseg64(unsigned long vid, unsigned long d) : value(vid << dbits | d) {}

    vdv_bitseg64(unsigned long vid, unsigned long d, unsigned long neib) : value(
            neib << (dbits + vbits) | vid << dbits | d) {}

    unsigned int v() const {
        return (value >> (dbits)) & vMASK;
    }

    unsigned int d() const {
        return value & dMASK;
    }

    unsigned int neib() const {
        return value >> (dbits + vbits);
    }
};

template<>
std::string str(const vdv_bitseg64 &obj) {
    return fmt::format("({},{},{})", obj.v(), obj.d(), obj.neib());
}

unsigned int vdv_bitseg64::dbits;
unsigned int vdv_bitseg64::vbits;
unsigned long vdv_bitseg64::dMASK;
unsigned long vdv_bitseg64::vMASK;
*/
struct TreeNode
{
  /* id: corresponding vertex id(also indexed by id)
   * f: parrent, if -1, it is a root or core.
   * h: nodes# from root to this node, root is 1,
   * rid: root id,
   * rsize: root bag size
   * w: node size.
   */
  int f, rid, rsize;
  vector<int> nbr; // of size <=w
  vector<vd_bitseg> cost;
  vector<int> anc;  // of size h
  vector<dint> dis; // of size h
  vector<int> ch;   // children
};

struct BPLabel
{

  BPLabel()
  {
    for (int i = 0; i < MAX_N_ROOTS; ++i)
    {
      bpspt_d[i] = MAXD;
      bpspt_s[i] = {0, 0};
    }
  }

  array<uint8_t, MAX_N_ROOTS> bpspt_d{};
  array<pair<uint64_t, uint64_t>, MAX_N_ROOTS> bpspt_s;
};

vector<vector<unsigned int>> load_graph(const string &path)
{
  vector<vector<unsigned int>> g;
  ifstream ifs(path);
  if (!ifs)
  {
    mlog("open file %s failed.", path.c_str());
    exit(EXIT_FAILURE);
  }
  // skip comments
  string line;
  ifstream::pos_type pos;
  do
  {
    pos = ifs.tellg();
    getline(ifs, line);
  } while (line[0] == '#');
  ifs.seekg(pos);
  // read graph
  int n, m, s, t;
  ifs >> n >> m;
  g.resize(n);
  while (ifs >> s >> t)
  {
    g[s].push_back(t);
  }
  return g;
}
/**
 * document outside
 */
class mll_ct
{
public:
  // dataset and graph
  vector<vector<unsigned int>> g;
  vector<int> id_map;
  vector<int> n2o_id;
  unsigned int n{0};
  unsigned long m{0};

  // bp labels
  array<vector<unsigned int>, MAX_N_ROOTS> rses;
  vector<unsigned char> usd_bp; // if a vertex is BP label
  vector<BPLabel> label_bp{};

  // tree label
  vector<TreeNode> tree{};
  vector<int> ord{}, rank{};
  bool tree_label_global{false};

  // core label
  int n_core{0};
  //    vector<int> cid2vid{};
  //    vector<int> vid2cid{};
  vector<vector<vd_bitseg>> corL{};

  // other
  vector<vector<pair<int, vd_bitseg>>> ctGraph;
  vector<unsigned int> r2v, v2r;

  // pTD label {v,d}-{vertex on the path}
  vector<vector<pair<vd_bitseg, unsigned int>>> ptdL;

  void load_graph()
  {
    mlog("load graph %s", params.dataset.c_str());
    string path = filesystem::path(params.index_folder).append(params.dataset + "-dis.bin").string();
    string path2 = filesystem::path(params.index_folder).append(params.dataset + "-idmap.bin").string();

    if (!graph_ops::load_graph_bin(path, g))
    {
      g = ::load_graph("../graphs/" + params.dataset + ".edges");
      graph_ops::convert_to_undirected(g);
      //      print_graph();
      id_map = graph_ops::equivalance_reduction(g);
      auto r = graph_ops::order_vertexes(g);
      graph_ops::relabel_graphs(g, r);
      //      print_graph();
      graph_ops::store_graph_bin(path, g);
      for (int &v : id_map)
      {
        if (v >= 0 && v < id_map.size())
        {
          v = r[v];
        }
      }
      graph_ops::store_vector_bin(path2, id_map);
    }
    else
    {
      graph_ops::load_vector_bin(path2, id_map);
    }
    n = g.size();
    m = 0;
    for (auto &v : g)
      m += v.size();
  }

public:
  explicit mll_ct(bool lgraph = true)
  { // load graph?
    mlog("CT: dataset is %s.", params.dataset.c_str());
    if (lgraph)
    {
      load_graph();
      mlog("CT: %s contains %d vertices and %lu edges after procession.", params.dataset.c_str(), n, m);
    }
  }

  vector<int> scc;

  int scc_cnt()
  {
    int sccid = 0;
    scc.resize(n, -1);

    queue<int> q;
    for (int i = 0; i < n; ++i)
    {
      if (scc[i] == -1)
      {
        scc[i] = sccid;
        q.push(i);
        while (!q.empty())
        {
          auto u = q.front();
          q.pop();
          for (const auto &v : g[u])
          {
            if (scc[v] == -1)
            {
              scc[v] = sccid;
              q.push(v);
            }
          }
        }
        sccid++;
      }
    }
    mlog("there are %d sccs in the graph", sccid);
    return sccid;
  }

  void construct_ctl_ptdl_seperately()
  {
    std::ifstream resultsFile;
    resultsFile.open("../results/" + params.dataset + "/MLL/" + params.dataset + "_MLL_initialization.txt");

    high_resolution_clock::time_point start_clock;
    high_resolution_clock::time_point end_clock;
    duration<double, std::milli> runtime = 0;
    duration<double, std::milli> tempRuntime = 0;
    if (!load_label_bp())
    {
      start_clock = high_resolution_clock::now();
      construct_bp_label();
      end_clock = high_resolution_clock::now();
      save_label_bp();
      runtime += end_clock - start_clock;
      tempRuntime = end_clock - start_clock;
      resultsFile << "construct_bp_label microseconds " << (double)tempRuntime.count() << "\n";
    }
    // decompose_tree
    if (!load_label_tree())
    {
      start_clock = high_resolution_clock::now();
      reduce();
      create_tree();
      end_clock = high_resolution_clock::now();
      save_label_tree();
      save_tmp_graph();
      runtime += end_clock - start_clock;
      tempRuntime = end_clock - start_clock;
      resultsFile << "decompose_tree microseconds " << (double)tempRuntime.count() << "\n";
    }
    else
    {
      load_tmp_graph();
    }
    // may remove the memory used by graph
    g.clear(), g.shrink_to_fit();
    // decompose_core
    if (!load_label_core())
    {
      //      compute_ct_ptd_label_in_core();
      start_clock = high_resolution_clock::now();
      compute_ctL_in_core();
      compute_ptdL_in_core();
      compute_ptdL_for_bpcore();
      end_clock = high_resolution_clock::now();
      // may remove the memory used by ctGraph
      ctGraph.clear(), ctGraph.shrink_to_fit();
      save_label_core();
      runtime += end_clock - start_clock;
      tempRuntime = end_clock - start_clock;
      resultsFile << "decompose_core microseconds " << (double)tempRuntime.count() << "\n";
    }
    else
    {
      load_ptdL();
    }
    // compute global tree label
    if (!tree_label_global)
    {
      start_clock = high_resolution_clock::now();
      vector<vector<char>> dmtxes;
      compute_ctL_in_tree(dmtxes);
      compute_ptdL_in_tree(dmtxes);
      sort_ptdL();
      end_clock = high_resolution_clock::now();
      save_label_tree();
      save_ptdL();
      runtime += end_clock - start_clock;
      tempRuntime = end_clock - start_clock;
      resultsFile << "compute_global_tree_label microseconds " << (double)tempRuntime.count() << "\n";
    }
    resultsFile << "total " << runtime.count() << "\n";
    resultsFile.close();
  }

  void construct_index()
  {
    // decompose_bp
    if (!load_label_bp())
    {
      construct_bp_label();
      save_label_bp();
    }
    // decompose_tree
    if (!load_label_tree())
    {
      reduce();
      create_tree();
      save_label_tree();
      save_tmp_graph();
    }
    else
    {
      load_tmp_graph();
    }
    // may remove the memory used by graph
    g.clear(), g.shrink_to_fit();
    // decompose_core
    if (!load_label_core())
    {
      compute_ct_ptd_label_in_core();
      compute_ptdL_for_bpcore();
      // may remove the memory used by ctGraph
      ctGraph.clear(), ctGraph.shrink_to_fit();
      save_label_core();
    }
    else
    {
      load_ptdL();
    }
    // compute global tree label
    if (!tree_label_global)
    {
      compute_ctL_ptL_in_tree();
      sort_ptdL();
      save_label_tree();
      save_ptdL();
    }
    //    scc_cnt();
  }

  void construct_bp_label()
  {
    mlog("CT: constructing BP labels of %d roots", params.n_roots);
    auto stime = omp_get_wtime();
    usd_bp.resize(n, false);
    //        vector<int> rses[N_ROOTS];// r and maximum 64 neighbors.

    // select N_ROOTS root vertices and max 64 neighbors for each root
    int r = 0;
    for (int i_bpspt = 0; i_bpspt < params.n_roots; ++i_bpspt)
    {
      while (r < n && usd_bp[r])
        ++r;
      if (r == n)
        break;
      usd_bp[r] = true;
      rses[i_bpspt].push_back(r);
      int ns = 0;
      for (auto v : g[r])
      {
        if (!usd_bp[v])
        {
          usd_bp[v] = true;
          rses[i_bpspt].push_back(v);
          if (++ns == 64)
            break;
        }
      }
    }

    if (params.report_level >= report_level::alldata)
    {
      mlog("BP roots and Srs:");
      for (int i = 0; i < params.n_roots; ++i)
      {
        printf("r%i=%s\n", i, str(rses[i]).c_str());
      }
    }

    // Construct BP labels.
    label_bp.resize(n);
    omp_set_num_threads(min(min(params.num_threads, params.n_roots), MAX_BP_THREADS));
#pragma omp parallel
    {
      //            if (omp_get_thread_num() == 0) mlog("{} threads used in creating BP labels", omp_get_num_threads())
      // #pragma omp for
      //            for (int i = 0; i < n; ++i) {
      //                for (int j = 0; j < N_ROOTS; ++j) {
      //                    label_bp[i].bpspt_d[j] = MAXD;
      //                }
      //            }

      vector<int> Q(n);

#pragma omp for schedule(dynamic)
      for (int i_bpspt = 0; i_bpspt < params.n_roots; ++i_bpspt)
      {
        if (rses[i_bpspt].empty())
          continue;

        vector<uint8_t> tmp_d(n, MAXD);
        vector<pair<uint64_t, uint64_t>> tmp_s(n, {0, 0});

        // init queue
        r = rses[i_bpspt][0];
        int que_t0 = 0, que_t1 = 0, que_h = 0;
        Q[que_h++] = r;
        tmp_d[r] = 0;
        que_t1 = que_h;

        for (size_t i = 1; i < rses[i_bpspt].size(); ++i)
        {
          int v = rses[i_bpspt][i];
          Q[que_h++] = v;
          tmp_d[v] = 1;
          tmp_s[v].first = 1ULL << (i - 1);
        }

        for (int d = 0; que_t0 < que_h; ++d)
        {
          for (int i = que_t0; i < que_t1; ++i)
          {
            int v = Q[i];
            int td = d + 1;

            for (auto tv : g[v])
            {
              if (d == tmp_d[tv])
              {
                tmp_s[tv].second |= tmp_s[v].first;
                tmp_s[v].second |= tmp_s[tv].first;
              }
            }
            for (const auto &tv : g[v])
            {
              if (d < tmp_d[tv])
              {
                if (tmp_d[tv] == MAXD)
                {
                  Q[que_h++] = tv;
                  tmp_d[tv] = td;
                }
                tmp_s[tv].first |= tmp_s[v].first;
                tmp_s[tv].second |= tmp_s[v].second;
              }
            }
          }

          que_t0 = que_t1;
          que_t1 = que_h;
        }

        for (int v = 0; v < n; ++v)
        {
          label_bp[v].bpspt_d[i_bpspt] = tmp_d[v];
          label_bp[v].bpspt_s[i_bpspt].first = tmp_s[v].first;
          label_bp[v].bpspt_s[i_bpspt].second = tmp_s[v].second & ~tmp_s[v].first;
        }
      }
    }
    //        exit(-1);
    mlog("CT: BP Label Constructed, bp_size=%fMB, t =%fsecs", sizeof(BPLabel) * n / (1024.0 * 1024.0),
         omp_get_wtime() - stime);
  }

  /* usd_bp, label_bp, rses
   */
  void save_label_bp()
  {
    mlog("CT: Saving BP Label...");
    FILE *fout = write_file(filesystem::path(params.index_folder)
                                .append(params.dataset + ".label-bp-" + to_string(params.n_roots) + ".bin")
                                .string());

    write(fout, usd_bp);
    write(fout, label_bp);
    for (int i = 0; i < params.n_roots; ++i)
    {
      write(fout, rses[i]);
    }

    fclose(fout);
    mlog("CT: BP Label Saved!");
  }

  /* usd_bp, label_bp, rses
   */
  bool load_label_bp()
  {
    mlog("Loading BP Label...");
    FILE *fin = read_file(filesystem::path(params.index_folder)
                              .append(params.dataset + ".label-bp-" + to_string(params.n_roots) + ".bin")
                              .string(),
                          false);

    if (fin == nullptr)
      return false;
    read(fin, usd_bp);
    n = usd_bp.size();
    read(fin, label_bp);
    for (int i = 0; i < params.n_roots; ++i)
    {
      read(fin, rses[i]);
    }
    check_eof(fin);
    return true;
  }

  bool prune_by_bp(int u, int v, int d)
  {
    BPLabel &idx_u = label_bp[u], &idx_v = label_bp[v];
    for (int i = 0; i < params.n_roots; ++i)
    {
      int td = idx_u.bpspt_d[i] + idx_v.bpspt_d[i];
      if (td - 2 <= d)
        td += (idx_u.bpspt_s[i].first & idx_v.bpspt_s[i].first) ? -2
              : ((idx_u.bpspt_s[i].first & idx_v.bpspt_s[i].second) | (idx_u.bpspt_s[i].second & idx_v.bpspt_s[i].first))
                  ? -1
                  : 0;
      if (td <= d)
        return true;
    }
    return false;
  }

  void reduce()
  {
    mlog("CT: reduce graph with max_w %d", params.max_w);
    auto stime = omp_get_wtime();
    // #define old
    // output structures<<<<<<<<<<<<<<<<<<<<<<
    ctGraph.resize(n);
    rank.resize(n, -1);
    // ord
    n_core = 0;
    unsigned long m_core = 0;
    //>>>>>>>>>>>>>>>
    vd_bitseg::init_static_byncore(n);

    vector<unordered_map<int, vd_bitseg>> E2(n);
    omp_set_num_threads(params.num_threads);
    int r = 0;

    // fill deg2vs, cur_deg and E2
    vector<vector<unsigned int>> deg2vs(params.max_w);
    vector<int> cur_deg(n);
    for (unsigned int i = 0; i < n; ++i)
    {
      cur_deg[i] = g[i].size();
      if (g[i].size() < params.max_w)
      {
        deg2vs[g[i].size()].push_back(i);
      }

      for (int j = 0; j < g[i].size(); ++j)
      {
        E2[i][g[i][j]] = {i, 1};
      }
    }

    int cur_min = 0;
    while (true)
    {
      /*
       * find a vertex to reduce.
       */
      unsigned int x;
      while (true)
      {
        while (cur_min < deg2vs.size() && deg2vs[cur_min].empty())
          cur_min++;
        if (cur_min >= params.max_w)
          goto out;
        x = deg2vs[cur_min].back();
        deg2vs[cur_min].pop_back();
        if (usd_bp[x])
          continue;
        if (cur_deg[x] == cur_min && rank[x] == -1)
          break; // found
      }

      ord.push_back(x);
      rank[x] = n - ++r;

      /*
       * reduce this vertex
       * 1. create the node E2[x]
       * 2. create clique
       */
      // step1
      for (auto &e : E2[x])
      {
        if (rank[e.first] == -1)
        {
          ctGraph[x].emplace_back(e);
        }
      }
      E2[x].clear();
      // step2
      for (const auto &s : ctGraph[x])
      {
        int deg_inc = -1;
        for (const auto &item : ctGraph[x])
        {
          if (s.first != item.first)
          {
            auto it = E2[s.first].find(item.first);
            if (it == E2[s.first].end())
            {
              E2[s.first][item.first] = {x, s.second.d() + item.second.d()};
              deg_inc++;
            }
            else
            {
              if (E2[s.first][item.first].d() >= s.second.d() + item.second.d())
              {
                E2[s.first][item.first] = {x, s.second.d() + item.second.d()};
              }
            }
          }
        }
        if (deg_inc != 0)
        {
          cur_deg[s.first] += deg_inc;
          if (cur_deg[s.first] < params.max_w)
          {
            deg2vs[cur_deg[s.first]].push_back(s.first);
          }
        }
      }
    }

  out:

    mlog("reduce finished. now collecting core graph...");

#pragma omp parallel
    {

#pragma omp for schedule(dynamic)
      for (int u = 0; u < n; ++u)
      {
        if (rank[u] != -1)
          continue;
        vector<pair<int, vd_bitseg>> tmp_edges;
        for (auto &e : E2[u])
        {
          if (rank[e.first] != -1)
            continue;
          tmp_edges.emplace_back(e);
        }
        //                sort(tmp_edges.begin(), tmp_edges.end());
        swap(tmp_edges, ctGraph[u]);
      }
    }

    for (int u = 0; u < n; ++u)
    {
      if (rank[u] != -1)
        continue;
      ++n_core;
      m_core += ctGraph[u].size();
    }
    mlog("collection finished, t=%f secs. n_core=%d, m_core=%lu, node_rate=%f, edge_rate=%f",
         omp_get_wtime() - stime, n_core, m_core, n_core * 1.0 / n, m_core * 1.0 / m);
  }

  /*
   * literally, create the tree structure.
   */
  void create_tree()
  {
    mlog("Creating Tree...");
    auto stime = omp_get_wtime();
    tree.resize(n);
    int cnt_root = 0, maxdep = 0, max_sub_tree = 1;
    vector<int> tcnt(n, 0);
    double tw = 0;
    for (int i = (int)ord.size() - 1; i >= 0; --i)
    {
      int x = ord[i];
      TreeNode &tn = tree[x];
      tn.nbr.resize(ctGraph[x].size());
      tn.cost.resize(ctGraph[x].size());
      for (int j = 0; j < (int)ctGraph[x].size(); ++j)
      {
        tn.nbr[j] = ctGraph[x][j].first;
        tn.cost[j] = ctGraph[x][j].second;
      }
      ctGraph[x].clear(), ctGraph[x].shrink_to_fit();
      tn.f = -1;
      for (auto &u : tn.nbr)
      {
        if (rank[u] != -1 && (tn.f == -1 || rank[u] > rank[tn.f]))
          tn.f = u;
      }
      if (tn.f == -1)
      {
        ++cnt_root;
        ++tcnt[x];
        tn.rid = x;
        tn.rsize = tn.nbr.size();
        tn.anc.push_back(x);
      }
      else
      {
        tree[tn.f].ch.push_back(x);
        tn.rid = tree[tn.f].rid;
        ++tcnt[tn.rid];
        max_sub_tree = max(max_sub_tree, tcnt[tn.rid]);
        tn.rsize = tree[tn.f].rsize;
        tn.anc = tree[tn.f].anc;
        tn.anc.push_back(x);
      }
      tw += tn.rsize;
      maxdep = max(maxdep, (int)tn.anc.size());
    }
    mlog("Core tree constructed, maxdep=%d, cnt_root=%d, max_stree=%d, avg_rsize=%f, t=%f secs", maxdep,
         cnt_root, max_sub_tree, tw / (n - n_core), omp_get_wtime() - stime);
  }

  vector<char> cache_interface_distances(const vector<int> &interface, vector<char> &tmp_dis)
  {
    int mtxn = interface.size();
    vector<char> d(mtxn * mtxn, -1);
    for (int i = 0; i < mtxn; ++i)
    {
      d[i * mtxn + i] = 0;
      for (const auto &l : corL[interface[i]])
        tmp_dis[l.v()] = l.d();

      for (int j = i + 1; j < mtxn; ++j)
      {
        int mind = 127;
        for (const auto &l : corL[interface[j]])
        {
          if (tmp_dis[l.v()] >= 0)
          {
            mind = min((unsigned int)mind, tmp_dis[l.v()] + l.d());
          }
        }
        int d2 = query_by_bp(interface[i], interface[j]);
        if (d2 < MAXD || mind < 127)
          mind = min(d2, mind);
        else
          mind = -1;
        d[i * mtxn + j] = mind;
        d[j * mtxn + i] = mind;
        //                assert(query_by_core(inter)
      }
      for (const auto &l : corL[interface[i]])
        tmp_dis[l.v()] = -1;
    }
    return d;
  }
  void compute_ptdL_in_tree(const vector<vector<char>> &dmtxes)
  {
    mlog("Computing ptd label in tree ...");
    auto stime = omp_get_wtime();
#pragma omp parallel
    {
      vector<char> dis_tmp(n_core, -1);
      vector<int> idx_s(n); // vid-> position in s
      vector<int> s;        // position in tree-> vid
#pragma omp for schedule(dynamic)
      for (int v = 0; v < n; ++v)
      {
        if (rank[v] >= 0 && tree[v].f == -1)
        {

          for (int neib : tree[v].nbr)
          {
            idx_s[neib] = s.size();
            s.push_back(neib);
          }
          compute_ptdL_in_tree(v, tree[v].rsize, s, idx_s, dmtxes[v]);
          s.clear();
        }
      }
    }

    tree_label_global = true;
    auto etime = omp_get_wtime();
    printf("ptd label in tree Computed, t=%0.3lf secs\n", etime - stime);
  }

  void compute_ptdL_in_tree(int x, int rsize, vector<int> &s, vector<int> &idx_s, const vector<char> &dmtx)
  {
    idx_s[x] = s.size();
    s.push_back(x);
    TreeNode &tn = tree[x];

    {
      // set cost
      vector<pair<vd_bitseg, unsigned int>> d2v(tn.nbr.size());
      for (int i = 0; i < tn.nbr.size(); ++i)
      {
        d2v[i] = {vd_bitseg(tn.cost[i].v(), tn.dis[idx_s[tn.nbr[i]]]), tn.nbr[i]};
      }
      sort(d2v.begin(), d2v.end(),
           [](const pair<vd_bitseg, unsigned int> &a, const pair<vd_bitseg, unsigned int> &b)
           {
             return a.first.d() < b.first.d();
           });
      vector<bool> pruned(d2v.size(), false);
      for (int i = 0; i < d2v.size(); ++i)
      {
        for (int j = 0; j < i; ++j)
        {
          if (pruned[j])
            continue;
          if (d2v[j].first.d() == d2v[i].first.d())
            break;
          int nowdis = -1;
          int a = max(idx_s[d2v[i].second], idx_s[d2v[j].second]);
          int b = min(idx_s[d2v[i].second], idx_s[d2v[j].second]);
          if (a >= rsize)
            nowdis = tree[s[a]].dis[b];
          else
            nowdis = dmtx[a * rsize + b];
          if (nowdis >= 0 && d2v[j].first.d() + nowdis == d2v[i].first.d())
          {
            pruned[i] = true;
            break;
          }
        }
        if (!pruned[i])
        {
          ptdL[x].emplace_back(vd_bitseg(d2v[i].second, d2v[i].first.d()), d2v[i].first.v());
        }
      }
    }

    tn.dis.back() = 0;
    for (int &u : tree[x].ch)
    {
      compute_ptdL_in_tree(u, rsize, s, idx_s, dmtx);
    }
    s.pop_back();
  }

  void compute_ctL_in_tree(vector<vector<char>> &dmtxes)
  {
    mlog("Computing Tree Label...");
    auto stime = omp_get_wtime();
    dmtxes.resize(n);
#pragma omp parallel
    {
      vector<char> dis_tmp(n_core, -1);
      vector<int> idx_s(n); // vid-> position in s
      vector<int> s;        // position in tree-> vid
#pragma omp for schedule(dynamic)
      for (int v = 0; v < n; ++v)
      {
        if (rank[v] >= 0 && tree[v].f == -1)
        {

          dmtxes[v] = cache_interface_distances(tree[v].nbr, dis_tmp);

          for (int neib : tree[v].nbr)
          {
            idx_s[neib] = s.size();
            s.push_back(neib);
          }
          compute_ctL_in_tree(v, tree[v].rsize, s, idx_s, dmtxes[v]);
          s.clear();
        }
      }
    }
    double t_size = 0;
    int maxdis = 0;

    for (int v = 0; v < n; ++v)
    {
      if (rank[v] >= 0)
      {
        t_size += tree[v].dis.size() * 1.0 * (sizeof(int) + sizeof(dint));
        for (auto &d : tree[v].dis)
          maxdis = max(maxdis, (int)d);
      }
    }
    tree_label_global = true;
    auto etime = omp_get_wtime();
    printf("ct label in tree Computed, t=%0.3lf secs, maxdis=%d, tree label size=%0.3lf MB\n",
           etime - stime, maxdis, t_size / (1024.0 * 1024.0));
  }

  void compute_ctL_in_tree(int x, int rsize, vector<int> &s, vector<int> &idx_s, const vector<char> &dmtx)
  {
    idx_s[x] = s.size();
    s.push_back(x);
    TreeNode &tn = tree[x];
    tn.dis.resize(s.size());

    for (int i = 0; i < tn.dis.size() - 1; ++i)
    {
      tn.dis[i] = -1;
      for (int j = 0; j < tn.nbr.size(); ++j)
      {
        int w = tn.cost[j].d(), k = idx_s[tn.nbr[j]], nowdis = -1;
        int a = max(i, k), b = min(i, k);
        if (a >= rsize)
        {
          nowdis = tree[s[a]].dis[b];
        }
        else
        {
          nowdis = dmtx[a * rsize + b];
        }

        if (nowdis >= 0 && (tn.dis[i] == -1 || nowdis + w < tn.dis[i]))
          tn.dis[i] = min(nowdis + w, int(MAXD));
      }
    }

    tn.dis.back() = 0;
    for (int &u : tree[x].ch)
    {
      compute_ctL_in_tree(u, rsize, s, idx_s, dmtx);
    }
    s.pop_back();
  }

  void compute_ctL_ptL_in_tree()
  {
    mlog("Computing Tree Label...");
    auto stime = omp_get_wtime();
#pragma omp parallel
    {
      vector<char> dis_tmp(n_core, -1);
      vector<int> idx_s(n); // vid-> position in s
      vector<int> s;        // position in tree-> vid
#pragma omp for schedule(dynamic)
      for (int v = 0; v < n; ++v)
      {
        if (rank[v] >= 0 && tree[v].f == -1)
        {
          vector<char> dmtx = cache_interface_distances(tree[v].nbr, dis_tmp);
          //                vector<char> dmtx;

          for (int neib : tree[v].nbr)
          {
            idx_s[neib] = s.size();
            s.push_back(neib);
          }
          compute_ctL_ptL_in_tree(v, tree[v].rsize, s, idx_s, dmtx);
          s.clear();
        }
      }
    }
    double t_size = 0;
    int maxdis = 0;

    for (int v = 0; v < n; ++v)
    {
      if (rank[v] >= 0)
      {
        t_size += tree[v].dis.size() * 1.0 * (sizeof(int) + sizeof(dint));
        //        ctGraph[v].clear();
        for (auto &d : tree[v].dis)
          maxdis = max(maxdis, (int)d);
      }
      //      ctGraph[v].shrink_to_fit();
    }
    tree_label_global = true;
    printf("Tree Label Computed, t=%0.3lf secs, maxdis=%d, tree label size=%0.3lf MB\n", omp_get_wtime() - stime,
           maxdis, t_size / (1024.0 * 1024.0));
  }

  void compute_ctL_ptL_in_tree(int x, int rsize, vector<int> &s, vector<int> &idx_s, const vector<char> &dmtx)
  {
    idx_s[x] = s.size();
    s.push_back(x);
    TreeNode &tn = tree[x];
    tn.dis.resize(s.size());

    for (int i = 0; i < tn.dis.size() - 1; ++i)
    {
      tn.dis[i] = -1;
      for (int j = 0; j < tn.nbr.size(); ++j)
      {
        int w = tn.cost[j].d(), k = idx_s[tn.nbr[j]], nowdis = -1;
        int a = max(i, k), b = min(i, k);
        if (a >= rsize)
          nowdis = tree[s[a]].dis[b];
        else
          nowdis = dmtx[a * rsize + b];

        if (nowdis >= 0 && (tn.dis[i] == -1 || nowdis + w < tn.dis[i]))
          tn.dis[i] = min(nowdis + w, int(MAXD));
      }
    }

    {
      // set cost
      vector<pair<vd_bitseg, unsigned int>> d2v(tn.nbr.size());
      for (int i = 0; i < tn.nbr.size(); ++i)
      {
        d2v[i] = {vd_bitseg(tn.cost[i].v(), tn.dis[idx_s[tn.nbr[i]]]), tn.nbr[i]};
      }
      sort(d2v.begin(),
           d2v.end(),
           [](const pair<vd_bitseg, unsigned int> &a, const pair<vd_bitseg, unsigned int> &b)
           {
             return a.first.d() < b.first.d();
           });
      vector<bool> pruned(d2v.size(), false);
      for (int i = 0; i < d2v.size(); ++i)
      {
        for (int j = 0; j < i; ++j)
        {
          if (pruned[j])
            continue;
          if (d2v[j].first.d() == d2v[i].first.d())
            break;
          int nowdis = -1;
          int a = max(idx_s[d2v[i].second], idx_s[d2v[j].second]);
          int b = min(idx_s[d2v[i].second], idx_s[d2v[j].second]);
          if (a >= rsize)
            nowdis = tree[s[a]].dis[b];
          else
            nowdis = dmtx[a * rsize + b];
          if (nowdis >= 0 && d2v[j].first.d() + nowdis == d2v[i].first.d())
          {
            pruned[i] = true;
            break;
          }
        }
        if (!pruned[i])
        {
          ptdL[x].emplace_back(vd_bitseg(d2v[i].second, d2v[i].first.d()), d2v[i].first.v());
        }
      }
    }

    tn.dis.back() = 0;
    for (int &u : tree[x].ch)
    {
      compute_ctL_ptL_in_tree(u, rsize, s, idx_s, dmtx);
    }
    s.pop_back();
  }

  void sort_ptdL()
  {
#pragma omp parallel for schedule(dynamic, 16) num_threads(params.num_threads)
    for (int u = 0; u < n; ++u)
    {
      sort(ptdL[u].begin(), ptdL[u].end(), [](auto const &a, auto const &b)
           { return a.first.v() < b.first.v(); });
    }
  }

  /*rank, tree,
   */
  void save_label_tree()
  {
    mlog("Saving Tree Label...");
    FILE *fout = write_file(filesystem::path(params.index_folder)
                                .append(params.dataset + ".label-tree-" + to_string(params.max_w) + ".bin")
                                .string());

    write(fout, tree_label_global);
    write(fout, rank);
    for (int i = 0; i < n; ++i)
    {
      if (rank[i] >= 0)
      {
        TreeNode &tn = tree[i];
        //        write(fout, tn.f);
        write(fout, tn.rid);
        write(fout, tn.rsize);
        write(fout, tn.nbr);
        write(fout, tn.cost);
        write(fout, tn.anc);
        write(fout, tn.dis);
        //        write(fout, tn.ch);
      }
    }
    fclose(fout);

    mlog("Tree Label Saved!");
  }

  /*rank, tree,
   */
  bool load_label_tree()
  {
    mlog("Loading Tree Label...");
    FILE *fin =
        read_file(filesystem::path(params.index_folder)
                      .append(params.dataset + ".label-tree-" + to_string(params.max_w) + ".bin")
                      .string(),
                  false);
    if (fin == nullptr)
      return false;

    read(fin, tree_label_global);
    read(fin, rank);
    tree.resize(rank.size());
    for (int i = 0; i < n; ++i)
    {
      if (rank[i] >= 0)
      {
        TreeNode &tn = tree[i];
        //        read(fin, tn.f);
        read(fin, tn.rid);
        read(fin, tn.rsize);

        read(fin, tn.nbr);
        read(fin, tn.cost);
        read(fin, tn.anc);
        read(fin, tn.dis);
        //        read(fin, tn.ch);
      }
    }
    check_eof(fin);
    return true;
  }

  /* n_core, ctGraph
   */
  void save_tmp_graph()
  {
    printf("Saving Tmp Graph...\n");
    FILE *fout = write_file(filesystem::path(params.index_folder)
                                .append(params.dataset + ".tmp-" + to_string(params.max_w) + ".bin")
                                .string());

    write(fout, n);
    write(fout, n_core);
    write(fout, rank);
    write(fout, usd_bp);
    write(fout, label_bp);

    for (int i = 0; i < n; ++i)
    {
      if (rank[i] >= 0)
        continue;
      write(fout, ctGraph[i]);
    }

    fclose(fout);
    printf("Tmp Graph Saved!\n");
  }

  /* n_core, ctGraph
   */
  void load_tmp_graph()
  {
    mlog("Loading Tmp Graph...");

    FILE
        *fin = read_file(filesystem::path(params.index_folder)
                             .append(params.dataset + ".tmp-" + to_string(params.max_w) + ".bin")
                             .string());

    read(fin, n);
    read(fin, n_core);
    read(fin, rank);
    read(fin, usd_bp);
    read(fin, label_bp);

    ctGraph.resize(n);
    for (int i = 0; i < n; ++i)
    {
      if (rank[i] >= 0)
        continue;
      read(fin, ctGraph[i]);
    }

    check_eof(fin);
    mlog("Tmp Graph Loaded!");
  }

  void compute_ct_ptd_label_in_core()
  {
    mlog("compute core label with %d threads.", params.num_threads);
    auto stime1 = omp_get_wtime();

    //        omp_set_num_threads(params.num_threads);
    if (n_core == 0)
    {
      mlog("No core nodes!");
      return;
    }

    vd_bitseg::init_static_byncore(n);

    vector<vector<int>> pos(n);
    vector<vector<int>> min_level(n);

    r2v.resize(n);
    v2r.resize(n);
    corL.resize(n);
    ptdL.resize(n);

    int p = 0;
    for (int i = 0; i < rses.size(); ++i)
    {
      for (int j = 0; j < rses[i].size(); ++j)
      {
        r2v[p] = rses[i][j];
        v2r[rses[i][j]] = p;
        p++;
      }
    }

    for (int u = 0; u < n; ++u)
    {
      if (usd_bp[u])
        continue;
      if (rank[u] == -1)
      {
        r2v[p] = u;
        v2r[u] = p;

        corL[u].emplace_back(p, 0);
        min_level[u].emplace_back(n);
        pos[u].push_back(1);

        p++;
      }
      else
      {
        r2v[rank[u]] = u;
        v2r[u] = rank[u];
      }
    }

    int dis = 0;

    vector<vector<pair<vd_bitseg, unsigned int>>> label_new(n);
    long long cnt = 1;
    double stime2;
#pragma omp parallel num_threads(params.num_threads)
    { // parallel
      vector<int> min_level_local(n_core, -2);
      vector<int> cand;
      vector<char> nowdis(n_core, -1);

      while (cnt > 0)
      { // for dis
#pragma omp master
        {
          if (dis == 1)
            mlog("%d threads is used for creating core label", omp_get_num_threads());
          dis++;
          stime2 = omp_get_wtime();
        }
#pragma omp barrier

#pragma omp master
        cnt = 0;

        long long local_cnt = 0;

#pragma omp for schedule(static, 1)
        for (int u = 0; u < n; ++u)
        { // u
          if (rank[u] >= 0 || usd_bp[u])
            continue;
          for (auto &l : corL[u])
            nowdis[l.v()] = l.d();
          for (int i = 0; i < ctGraph[u].size(); ++i)
          { // neib
            int v = ctGraph[u][i].first, d = ctGraph[u][i].second.d();
            if (d > dis || usd_bp[v])
              continue;
            for (int j = d == dis ? 0 : pos[v][dis - d - 1]; j < pos[v][dis - d]; ++j)
            { // corv
              int w = corL[v][j].v();
              if (r2v[w] >= u)
                break;
              if (min_level_local[w] == -2)
              {
                if (!prune_by_bp(u, r2v[w], dis) && can_update(r2v[w], dis, nowdis))
                {
                  if (dis == d)
                  {
                    if (d == 1)
                    { // edge is from original graph
                      min_level_local[w] = v2r[u];
                    }
                    else
                    { // edge is created by vertex reduction
                      min_level_local[w] = rank[ctGraph[u][i].second.v()];
                    }
                  }
                  else
                  {
                    min_level_local[w] = min_level[v][j];
                  }
                }
                else
                {
                  min_level_local[w] = -1; // pruned.
                }
                cand.push_back(w);
              }
              else if (min_level_local[w] >= 0)
              {
                min_level_local[w] = min(min_level_local[w], min_level[v][j]);
              }
            }
          } // neib

          int n_cand = 0;
          for (int i = 0; i < cand.size(); ++i)
          {
            if (min_level_local[cand[i]] >= 0)
              cand[n_cand++] = cand[i];
            else
              min_level_local[cand[i]] = -2;
          }

          cand.resize(n_cand);
          sort(cand.begin(), cand.end());

          label_new[u].reserve(label_new[u].size() + cand.size());
          for (auto v : cand)
          {
            label_new[u].emplace_back(vd_bitseg(v, dis), min_level_local[v]);
            min_level_local[v] = -2;
          }
          local_cnt += cand.size();
          cand.clear();
          for (int i = 0; i < corL[u].size(); ++i)
            nowdis[corL[u][i].v()] = -1;
        } // u

#pragma omp critical
        {
          cnt += local_cnt;
        }

#pragma omp for schedule(static, 1)
        for (int u = 0; u < n; ++u)
        {
          if (rank[u] >= 0 || usd_bp[u])
            continue;
          sort(label_new[u].begin(), label_new[u].end(),
               [](const auto &a, const auto &b)
               { return a.first.v() < b.first.v(); });
          corL[u].reserve(corL[u].size() + label_new[u].size());
          min_level[u].reserve(min_level[u].size() + label_new[u].size());
          for (int i = 0; i < label_new[u].size(); ++i)
          {
            corL[u].emplace_back(label_new[u][i].first);
            min_level[u].emplace_back(min(label_new[u][i].second, v2r[u]));

            if (label_new[u][i].second >= v2r[u])
            {
              ptdL[u].emplace_back(vd_bitseg(r2v[label_new[u][i].first.v()], label_new[u][i].first.d()),
                                   r2v[label_new[u][i].second]);
            }
          }
          ptdL[u].shrink_to_fit();
          // vector<unsigned>(label[u]).swap(label[u]);
          vector<pair<vd_bitseg, unsigned int>>().swap(label_new[u]);
          pos[u].push_back((int)corL[u].size());
        }
        //
#pragma omp master
        {
          mlog("done! dis=%d,cnt=%lld,t=%f secs", dis, cnt, omp_get_wtime() - stime2);
        }
      } // while dis
    }   // omp parallel

    double tt = 0, max_label = 0;
    for (int i = 0; i < n; ++i)
    {
      if (rank[i] == -1)
      {
        tt += (double)corL[i].size() * 4;
        max_label = max(max_label, (double)corL[i].size());
      }
    }

    printf("Core label size=%0.3lfMB, Max core label size=%0.0lf, Avg core label size=%0.3lf, Time = %0.3lf sec\n",
           tt / (1024 * 1024.0), max_label, tt * 0.25 / n, omp_get_wtime() - stime1);
  }

  void compute_ctL_in_core()
  {
    mlog("compute ct core label with %d threads.", params.num_threads);
    auto stime1 = omp_get_wtime();

    if (n_core == 0)
    {
      mlog("No core nodes!");
      return;
    }

    vd_bitseg::init_static_byncore(n);

    vector<vector<int>> pos(n);

    corL.resize(n);

    // root vertices 2 rank, small value<->higher rank
    r2v.resize(n);
    v2r.resize(n);
    int p = 0;
    for (int i = 0; i < rses.size(); ++i)
    {
      for (int j = 0; j < rses[i].size(); ++j)
      {
        r2v[p] = rses[i][j];
        v2r[rses[i][j]] = p;
        p++;
      }
    }

    for (int u = 0; u < n; ++u)
    {
      if (usd_bp[u])
        continue; // bp vertex
      if (rank[u] == -1)
      { // core vertex
        r2v[p] = u;
        v2r[u] = p;
        corL[u].emplace_back(p, 0);
        pos[u].push_back(1);
        p++;
      }
      else
      { // tree vertex
        r2v[rank[u]] = u;
        v2r[u] = rank[u];
      }
    }

    int dis = 0;
    // u->[{v,d}]
    vector<vector<vd_bitseg>> label_new(n);
    long long cnt = 1;
    double stime2;
#pragma omp parallel num_threads(params.num_threads)
    { // parallel
      vector<char> processed(n_core, 0);
      vector<int> cand;
      vector<char> nowdis(n_core, -1); // map index (v,d) to an array nowdis[v]=d

      while (cnt > 0)
      { // for dis
#pragma omp single
        {
          if (dis == 0)
            mlog("%d threads is used for creating core label", omp_get_num_threads());
          dis++;
          stime2 = omp_get_wtime();
        }
#pragma omp barrier

#pragma omp master //--------------Barrier------------------------
        cnt = 0;

        long long local_cnt = 0;

#pragma omp for schedule(static, 1)
        for (int u = 0; u < n; ++u)
        { // u
          if (rank[u] >= 0 || usd_bp[u])
            continue;
          for (auto &l : corL[u])
            nowdis[l.v()] = (char)l.d();
          for (int i = 0; i < ctGraph[u].size(); ++i)
          { // neib
            int v = ctGraph[u][i].first, d = (int)ctGraph[u][i].second.d();
            if (d > dis || usd_bp[v])
              continue;
            for (int j = d == dis ? 0 : pos[v][dis - d - 1]; j < pos[v][dis - d]; ++j)
            { // corv
              int w = (int)corL[v][j].v();
              if (r2v[w] >= u)
                break; // w has lower rank than u;
              if (!processed[w])
              {
                if (!prune_by_bp(u, (int)r2v[w], dis) && can_update((int)r2v[w], dis, nowdis))
                {
                  processed[w] = 1; //
                }
                else
                {
                  processed[w] = 2; // pruned.
                }
                cand.push_back(w);
              }
            }
          } // neib

          int n_cand = 0;
          for (int i = 0; i < cand.size(); ++i)
          {
            if (processed[cand[i]] == 1)
            {
              cand[n_cand++] = cand[i];
            }
            processed[cand[i]] = 0;
          }
          cand.resize(n_cand);
          sort(cand.begin(), cand.end());

          label_new[u].reserve(cand.size());
          for (auto v : cand)
          {
            label_new[u].emplace_back(vd_bitseg(v, dis));
          }
          local_cnt += cand.size();
          cand.clear();
          for (auto &l : corL[u])
            nowdis[l.v()] = -1;
        } // u

#pragma omp critical
        {
          cnt += local_cnt;
        }

#pragma omp for schedule(static, 1)
        for (int u = 0; u < n; ++u)
        {
          if (rank[u] >= 0 || usd_bp[u])
            continue;
          // seems already sorted in sort(cand)
          //           sort(label_new[u].begin(), label_new[u].end(),
          //                [](const auto &a, const auto &b) { return a.first.v() < b.first.v(); });
          corL[u].reserve(corL[u].size() + label_new[u].size());
          corL[u].insert(corL[u].end(), label_new[u].begin(), label_new[u].end());
          vector<vd_bitseg>().swap(label_new[u]);
          pos[u].push_back((int)corL[u].size());
          pos[u].shrink_to_fit();
        }
        //
#pragma omp master
        {
          mlog("done! dis=%d,cnt=%lld,t=%f secs", dis, cnt, omp_get_wtime() - stime2);
        }
      } // while dis
    }   // omp parallel
    auto etime1 = omp_get_wtime();

    double tt = 0, max_label = 0;
    for (int i = 0; i < n; ++i)
    {
      if (rank[i] == -1)
      {
        tt += corL[i].size() * 4;
        max_label = max(max_label, corL[i].size() * 1.0);
      }
    }

    printf("Core label size=%0.3lfMB, Max core label size=%0.0lf, Avg core label size=%0.3lf, Time = %0.3lf sec\n",
           tt / (1024 * 1024.0), max_label, tt * 0.25 / n, etime1 - stime1);
  }

  void compute_ptdL_in_core()
  {
    mlog("compute ptd label in core with %d threads.", params.num_threads);
    auto stime1 = omp_get_wtime();

    if (n_core == 0)
    {
      mlog("No core nodes!");
      return;
    }

    vd_bitseg::init_static_byncore(n);

    vector<vector<int>> pos(n);
    vector<vector<int>> min_level(n);

    //    corL.resize(n);
    ptdL.resize(n);

    // root vertices 2 rank, small value<->higher rank
    //     r2v.resize(n);
    //     v2r.resize(n);
    //     int p = 0;
    //     for (int i = 0; i < rses.size(); ++i) {
    //       for (int j = 0; j < rses[i].size(); ++j) {
    //         r2v[p] = rses[i][j];
    //         v2r[rses[i][j]] = p;
    //         p++;
    //       }
    //     }

    for (int u = 0; u < n; ++u)
    {
      if (usd_bp[u])
        continue; // bp vertex
      if (rank[u] == -1)
      { // core vertex
        //        r2v[p] = u;
        //        v2r[u] = p;
        //        corL[u].emplace_back(p, 0);
        min_level[u].emplace_back(n);
        pos[u].push_back(1);
        //        p++;
      }
      else
      { // tree vertex
        //        r2v[rank[u]] = u;
        //        v2r[u] = rank[u];
      }
    }

    int dis = 0;
    // u->[{v,d},midw]
    vector<vector<pair<vd_bitseg, unsigned int>>> label_new(n);
    long long cnt = 1;
    double stime2;
#pragma omp parallel num_threads(params.num_threads)
    { // parallel
      vector<int> min_level_local(n_core, -2);
      vector<int> cand;
      //      vector<char> nowdis(n_core, -1);//map index (v,d) to an array nowdis[v]=d

      while (cnt > 0)
      { // for dis
#pragma omp single
        {
          if (dis == 0)
            mlog("%d threads is used for creating core label", omp_get_num_threads());
          dis++;
          stime2 = omp_get_wtime();
        }
#pragma omp barrier

#pragma omp master //--------------Barrier------------------------
        cnt = 0;

        long long local_cnt = 0;

#pragma omp for schedule(static, 1)
        for (int u = 0; u < n; ++u)
        { // u
          if (rank[u] >= 0 || usd_bp[u])
            continue;
          //          for (auto &l : corL[u]) nowdis[l.v()] = l.d();
          for (size_t i = pos[u].back(); i < corL[u].size() && corL[u][i].d() == dis; ++i)
          {
            min_level_local[corL[u][i].v()] = n;
            cand.push_back(corL[u][i].v());
          }
          if (!cand.empty())
          {
            for (int i = 0; i < ctGraph[u].size(); ++i)
            { // neib
              int v = ctGraph[u][i].first, d = ctGraph[u][i].second.d();
              if (d > dis || usd_bp[v])
                continue;
              for (int j = d == dis ? 0 : pos[v][dis - d - 1]; j < pos[v][dis - d]; ++j)
              { // corv
                int w = corL[v][j].v();
                if (min_level_local[w] >= 0)
                {
                  if (dis == d)
                  {
                    if (d == 1)
                    { // edge is from original graph
                      min_level_local[w] = v2r[u];
                    }
                    else
                    { // edge is created by vertex reduction
                      min_level_local[w] = rank[ctGraph[u][i].second.v()];
                    }
                  }
                  else
                  {
                    min_level_local[w] = min(min_level_local[w], min_level[v][j]);
                  }
                }
              }
            } // neib
          }

          label_new[u].reserve(cand.size());
          for (auto v : cand)
          {
            label_new[u].emplace_back(vd_bitseg(v, dis), min_level_local[v]);
            min_level_local[v] = -2;
          }
          local_cnt += cand.size();
          cand.clear();
          //          for (int i = 0; i < corL[u].size(); ++i) nowdis[corL[u][i].v()] = -1;
        } // u

#pragma omp critical
        {
          cnt += local_cnt;
        }

#pragma omp for schedule(static, 1)
        for (int u = 0; u < n; ++u)
        {
          if (rank[u] >= 0 || usd_bp[u])
            continue;
          // seems already sorted
          //           sort(label_new[u].begin(), label_new[u].end(),
          //                [](const auto &a, const auto &b) { return a.first.v() < b.first.v(); });
          //           corL[u].reserve(corL[u].size() + label_new[u].size());
          min_level[u].reserve(min_level[u].size() + label_new[u].size());
          for (int i = 0; i < label_new[u].size(); ++i)
          {
            //            corL[u].emplace_back(label_new[u][i].first);
            min_level[u].emplace_back(min(label_new[u][i].second, v2r[u]));

            if (label_new[u][i].second >= v2r[u])
            {
              ptdL[u].emplace_back(vd_bitseg(r2v[label_new[u][i].first.v()], label_new[u][i].first.d()),
                                   r2v[label_new[u][i].second]);
            }
          }
          ptdL[u].shrink_to_fit();
          // vector<unsigned>(label[u]).swap(label[u]);
          vector<pair<vd_bitseg, unsigned int>>().swap(label_new[u]);
          pos[u].push_back((int)min_level[u].size());
          pos[u].shrink_to_fit();
        }
        //
#pragma omp master
        {
          mlog("done! dis=%d,cnt=%lld,t=%f secs", dis, cnt, omp_get_wtime() - stime2);
        }
      } // while dis
    }   // omp parallel
    auto etime1 = omp_get_wtime();

    printf("ptd label in Core construction done. Time = %0.3lf sec\n", etime1 - stime1);
  }

  void save_label_core()
  {
    printf("Saving Core Label...\n");
    FILE *fout =
        write_file(filesystem::path(params.index_folder)
                       .append(params.dataset + ".label-core-" + to_string(params.max_w) + ".bin")
                       .string());

    write(fout, n);
    write(fout, n_core);
    write(fout, v2r);
    write(fout, r2v);
    for (int i = 0; i < n; ++i)
    {
      if (rank[i] >= 0 || usd_bp[i])
        continue;
      write(fout, corL[i]);
    }

    fclose(fout);
    printf("Core Label Saved!\n");
  }

  bool load_label_core()
  {
    mlog("Loading Core Label...");
    FILE *fin =
        read_file(filesystem::path(params.index_folder)
                      .append(params.dataset + ".label-core-" + to_string(params.max_w) + ".bin")
                      .string(),
                  false);
    if (fin == nullptr)
      return false;

    read(fin, n);
    read(fin, n_core);
    read(fin, v2r);
    read(fin, r2v);
    corL.resize(n);

    //    double tme = omp_get_wtime();
    for (int i = 0; i < n; ++i)
    {
      if (rank[i] >= 0 || usd_bp[i])
        continue;

      read(fin, corL[i]);
      //      if ((omp_get_wtime() - tme) > 1) {
      //        mlog("reading {} finishd size is {} tme{}", i, corL[i].size(), tme)
      //        tme = omp_get_wtime();
      //      }
    }
    vd_bitseg::init_static_byncore(n);

    check_eof(fin);
    return true;
  }

  void save_ptdL() const
  {
    printf("Saving ptd Label...\n");
    FILE *fout =
        write_file(filesystem::path(params.index_folder)
                       .append(params.dataset + ".label-ptd-" + to_string(params.max_w) + ".bin")
                       .string());

    write(fout, ptdL);
    fclose(fout);
  }

  void load_ptdL()
  {
    mlog("Loading ptdL Label...");
    FILE *fin =
        read_file(filesystem::path(params.index_folder)
                      .append(params.dataset + ".label-ptd-" + to_string(params.max_w) + ".bin")
                      .string(),
                  false);
    if (fin == nullptr)
    {
      return;
    }

    read(fin, ptdL);

    check_eof(fin);
  }

  // return {dis,ridx,sidx}
  bool q_bpr2v_by_bp(pair<unsigned int, char> bpr, unsigned int v, char d)
  {
    BPLabel &idx_u = label_bp[bpr.first], &idx_v = label_bp[v];
    //        if (bpr.first == 0 && v == 3) {
    //            printf("=====================================>>><<<<<<=====================\n");
    //        }
    for (int i = 0; i <= bpr.second; ++i)
    {
      char td = idx_u.bpspt_d[i] + (int)idx_v.bpspt_d[i];
      if (td - 2 <= d)
      {
        if (idx_u.bpspt_s[i].first & idx_v.bpspt_s[i].first)
        {
          if (i < bpr.second)
            return true;
        }
        else if (idx_u.bpspt_s[i].first & idx_v.bpspt_s[i].second && idx_u.bpspt_s[i].second & idx_v.bpspt_s[i].first)
        {
          if (td - 1 <= d)
          {
            if (i < bpr.second)
            {
              return true;
            }
            else
            {
              unsigned long tmpbit = idx_u.bpspt_s[i].second & idx_v.bpspt_s[i].first;
              tmpbit = tmpbit - (tmpbit & (tmpbit - 1));
              if (tmpbit < idx_u.bpspt_s[i].first)
                return true;
            }
          }
        }
        else
        {
          if (td <= d && idx_u.bpspt_s[i].first)
            return true;
        }
      }
    }
    return false;
  }

  void compute_ptdL_for_bpcore()
  {
    mlog("compute ptdL for bp_core");
    auto stime = omp_get_wtime();
    vector<pair<unsigned int, char>> bpvs;
    for (int i = 0; i < rses.size(); ++i)
    {
      for (int j = 0; j < rses[i].size(); ++j)
      {
        bpvs.emplace_back(rses[i][j], i);
      }
    }

    vector<vector<pair<vd_bitseg, unsigned int>>> result(bpvs.size());

#pragma omp parallel for num_threads(params.num_threads)
    for (int i = 0; i < bpvs.size(); ++i)
    {
      auto u = bpvs[i];
      //            unsigned long bitrep = 1;
      //            mlog("index for vertex {}", u.first)
      vector<pair<unsigned int, unsigned int>> dnmin_level(n, {MAXD, n});

      vector<vector<unsigned int>> d2vQ;
      d2vQ.resize(1);
      d2vQ[0].push_back(u.first);
      dnmin_level[u.first] = {0, n};

      int d = 0;
      while (d < d2vQ.size())
      {
        for (const auto &v : d2vQ[d])
        { // Q.top->v
          //   mlog("{} toped from queue with d={}", v, d)
          if (d > dnmin_level[v].first)
            continue; // if old value, processed, skip the vertex
          if (q_bpr2v_by_bp(u, v, d))
          { // shortest path pass other bp vertexes
            continue;
          }
          for (int j = 0; j < ctGraph[v].size(); ++j)
          { // v.neib->w
            auto w = ctGraph[v][j];
            if (v2r[w.first] <= v2r[u.first])
              continue;
            //                        if (usd_bp[w.first])continue;
            if (dnmin_level[w.first].first > (d + w.second.d()))
            {
              // update
              dnmin_level[w.first].first = d + w.second.d();
              if (v == u.first)
              { // direct edge.
                if (w.second.d() > 1)
                  dnmin_level[w.first].second = v2r[w.second.v()];
              }
              else
              {
                dnmin_level[w.first].second = min(dnmin_level[v].second, v2r[v]);
              }
              // push to q
              d2vQ.resize(max((unsigned int)d2vQ.size(), d + w.second.d() + 1));
              d2vQ[d + w.second.d()].push_back(w.first);
            }
            else if (dnmin_level[w.first].first == (d + w.second.d()))
            {
              dnmin_level[w.first].second = min(dnmin_level[v].second, dnmin_level[w.first].second);
              dnmin_level[w.first].second = min(v2r[v], dnmin_level[w.first].second);
            }
          }
          // check and add to ptdL
          if (d > 0 && dnmin_level[v].second >= v2r[v])
          {
            if (dnmin_level[v].second < n)
              result[i].emplace_back(vd_bitseg(v, d), r2v[dnmin_level[v].second]);
            else
              result[i].emplace_back(vd_bitseg(v, d), v);
            //                        ptdL[v].emplace_back(vd_bitseg(u.first, d), r2v[dnmin_level[v].second]);
          }
        }
        vector<unsigned int>().swap(d2vQ[d]);
        d++;
      }
    } // end parallel

    for (int i = 0; i < result.size(); ++i)
    {
      auto u = bpvs[i];
      for (int j = 0; j < result[i].size(); ++j)
      {
        ptdL[result[i][j].first.v()]
            .emplace_back(vd_bitseg(u.first, result[i][j].first.d()), result[i][j].second);
      }
    }

    auto etime = omp_get_wtime();
    mlog("compute ptdL for bp_core done in %f seconds.", etime - stime);
  }
  //  static

  void load_label()
  {
    load_label_tree();
    load_label_bp();
    load_label_core();
    load_graph();
    load_ptdL();
  }

  int query_by_bp(int u, int v)
  {
    BPLabel &idx_u = label_bp[u], &idx_v = label_bp[v];
    int d = MAXD;
    for (int i = 0; i < params.n_roots; ++i)
    {
      int td = idx_u.bpspt_d[i] + (int)idx_v.bpspt_d[i];
      if (td - 2 <= d)
        td += (idx_u.bpspt_s[i].first & idx_v.bpspt_s[i].first) ? -2
              : ((idx_u.bpspt_s[i].first & idx_v.bpspt_s[i].second) | (idx_u.bpspt_s[i].second & idx_v.bpspt_s[i].first))
                  ? -1
                  : 0;
      if (td < d)
        d = td;
    }
    return d == MAXD ? MAXD : d;
  }

  int query_with_original_id(int u, int v)
  {
    //        fflush(stdout);
    if (u == v)
      return 0;
    int fu, fv;
    int eqtype = 0;
    if (id_map[u] < 0)
    {
      fu = -id_map[u] - 1;
      eqtype = 1;
    }
    else if (id_map[u] < id_map.size())
    {
      fu = u;
    }
    else
    {
      fu = id_map[u] - id_map.size();
      eqtype = 2;
    }

    if (id_map[v] < 0)
    {
      fv = -id_map[v] - 1;
      eqtype = 1;
    }
    else if (id_map[v] < id_map.size())
    {
      fv = v;
    }
    else
    {
      fv = id_map[v] - id_map.size();
      eqtype = 2;
    }

    if (fu == fv)
    {
      if (eqtype == 1)
        return 2;
      else if (eqtype == 2)
        return 1;
    }
    return query(id_map[fu], id_map[fv]);
  }

  int query(int u, int v)
  {
    if (u >= n || v >= n)
      return MAXD;

    int d = query_by_bp(u, v);
    if (usd_bp[u] || usd_bp[v])
      return d;
    if (rank[u] >= 0 && rank[v] >= 0 && tree[u].rid == tree[v].rid)
    {
      return min(d, query_by_tree(u, v));
    }
    else
    {
      return query_by_core(u, v, d);
    }
  }

  int query_by_core(int u, int v, int d)
  {
    //        if (usd_bp[u] || usd_bp[v]) return d;
    static tint nowt = 0;
    static vector<tint> last_t(n_core, 0);
    static vector<int> dis(n_core);
    ++nowt;
    if (nowt == MAXT)
    {
      fill(last_t.begin(), last_t.end(), 0);
      nowt = 1;
    }
    if (rank[v] == -1 && rank[u] >= 0)
      swap(v, u);
    else if (rank[v] >= 0 && rank[u] >= 0)
    {
      if (tree[v].rsize < tree[u].rsize)
        swap(v, u);
    }
    else if (rank[v] == -1 && rank[u] == -1)
    {
      if (corL[v].size() < corL[u].size())
        swap(v, u);
    }
    int mind = d;
    if (rank[u] == -1)
    {
      for (auto l : corL[u])
      {
        last_t[l.v()] = nowt;
        dis[l.v()] = l.d();
      }
    }
    else
    {
      TreeNode &tu = tree[u];
      TreeNode &r = tree[tu.rid];
      for (size_t i = 0; i < tu.rsize; ++i)
      {
        int x = r.nbr[i], w = tu.dis[i];
        // if(w >= mind) continue;
        for (auto l : corL[x])
        {
          int y = l.v(), nowd = l.d() + w;
          if (nowd >= mind)
            break;
          if (last_t[y] != nowt)
          {
            last_t[y] = nowt;
            dis[y] = nowd;
          }
          else if (nowd < dis[y])
            dis[y] = nowd;
        }
      }
    }
    if (rank[v] == -1)
    {
      for (auto l : corL[v])
      {
        if (last_t[l.v()] == nowt)
        {
          mind = min(mind, (int)(l.d()) + dis[l.v()]);
        }
      }
    }
    else
    {
      TreeNode &tv = tree[v];
      TreeNode &r = tree[tv.rid];
      for (size_t i = 0; i < tv.rsize; ++i)
      {
        int x = r.nbr[i], w = tv.dis[i];
        // if(w >= mind) continue;
        for (auto l : corL[x])
        {
          int nowd = l.d() + w;
          if (nowd >= mind)
            break;
          if (last_t[l.v()] == nowt)
          {
            mind = min(mind, nowd + dis[l.v()]);
          }
        }
      }
    }
    return mind;
  }

  int query_by_tree(int u, int v)
  {
    TreeNode &tu = tree[u], &tv = tree[v];
    int len = min(tu.dis.size(), tv.dis.size());
    int d = MAXD;
    for (int i = 0; i < len; ++i)
    {
      if (i >= tu.rsize && tu.anc[i - tu.rsize] != tv.anc[i - tv.rsize])
        break;
      d = min(d, tu.dis[i] + tv.dis[i]);
    }
    return d;
  }

  deque<tuple<int, int, int>> retrieve_ptd_path_to_p_on_core(unsigned int u, pair<int, int> dc)
  {
    deque<tuple<int, int, int>> pu;
    int inc_dis = 0;
    unsigned int cu = u;
    while (inc_dis < dc.first)
    {
      for (int i = 0; i < ptdL[cu].size(); ++i)
      {
        auto u0 = ptdL[cu][i].first.v(), d0 = ptdL[cu][i].first.d(), i0 = ptdL[cu][i].second;
        auto c_it =
            lower_bound(corL[u0].begin(), corL[u0].end(), make_pair(dc.first - inc_dis - int(d0), dc.second),
                        [](const vd_bitseg &a, const pair<int, int> &b)
                        {
                          if (a.d() < b.first)
                            return true;
                          else if (a.d() == b.first && a.v() < b.second)
                            return true;
                          return false;
                        });
        //                for (int j = 0; j < corL[u0].size(); ++j) //done:accelerate with ordered corL
        if (c_it != corL[u0].end() && c_it->v() == dc.second && dc.first == inc_dis + d0 + c_it->d())
        {
          pu.emplace_back(u0, d0, i0);
          cu = u0;
          inc_dis += d0;
          goto out;
        }
      }
    out:;
    }
    return pu;
  }

  //    deque<tuple<int, int, int>>
  //    retrieve_ptd_path_to_p_on_core(unsigned int u, int cr, int cd) {//done:pass distance to here
  ////        for (int i = 0; i < corL[u].size(); ++i) {//done:accelerate with sort and distance information
  //        auto c_it = lower_bound(corL[u].begin(), corL[u].end(), make_pair(cd, cr),
  //                                [](const vd_bitseg &a, const pair<int, int> &b) {
  //                                    if (a.d() < b.first)return true;
  //                                    else if (a.d() == b.first && a.v() < b.second)return true;
  //                                    return false;
  //                                });
  //        if (cr == corL[u][i].v())return retrieve_ptd_path_to_p_on_core(u, {corL[u][i].d(), cr});
  ////        }
  //    }

  /**
   *
   * @param u
   * @param v
   * @return {distance, common vertex rank in core}
   */
  tuple<int, int, int> retrieve_ptd_dist_on_core(unsigned int u, unsigned int v)
  {
    static unsigned short nowt = 0;
    static vector<unsigned short> last_t(n_core, 0);
    static vector<unsigned char> dis(n_core);
    static unsigned int previous_u = UINT32_MAX;

    if (previous_u != u)
    {
      if (nowt == MAXT)
      {
        fill(last_t.begin(), last_t.end(), 0);
        nowt = 0;
      }
      nowt++;
      previous_u = u;
      unsigned int pur;
      for (const auto &e : corL[u])
      {
        dis[e.v()] = e.d();
        last_t[e.v()] = nowt;
      }
    }

    unsigned char d = 255, ud = 255, vd = 255;
    unsigned int pur;

    for (const auto &e : corL[v])
    {
      if (last_t[e.v()] == nowt)
      {
        if (d > dis[e.v()] + e.d())
        {
          d = dis[e.v()] + e.d();
          ud = dis[e.v()];
          vd = e.d();
          pur = e.v();
        }
      }
    }
    return {ud, vd, pur};
  }

  deque<tuple<int, int, int>> retrieve_ptd_path_on_core(unsigned int u, unsigned int v, tuple<int, int, int> &dc)
  {
    auto u_p = retrieve_ptd_path_to_p_on_core(u, {get<1>(dc), get<2>(dc)});
    auto v_p = retrieve_ptd_path_to_p_on_core(v, {get<0>(dc), get<2>(dc)});
    for (int i = v_p.size() - 1; i > 0; --i)
    {
      u_p.emplace_back(get<0>(v_p[i - 1]), get<1>(v_p[i]), get<2>(v_p[i]));
    }
    if (v_p.size() > 0)
    {
      u_p.emplace_back(v, get<1>(v_p[0]), get<2>(v_p[0]));
    }
    return u_p;
  }

  deque<tuple<int, int, int>> retrieve_ptd_path_on_tree(unsigned int u, unsigned int u_itfc)
  {
    deque<tuple<int, int, int>> pu;
    auto cu = u;
    unsigned int inc_dis = 0;
    while (inc_dis < tree[u].dis[u_itfc])
    {
      if (rank[cu] != -1)
      {
        bool found = false;
        for (int i = 0; i < ptdL[cu].size(); ++i)
        {
          auto u0 = ptdL[cu][i].first.v(), d0 = ptdL[cu][i].first.d(), i0 = ptdL[cu][i].second;
          if (rank[u0] != -1 && tree[u0].dis.size() > u_itfc)
          {
            if (tree[u0].dis[u_itfc] + inc_dis + d0 == tree[u].dis[u_itfc])
            {
              cu = u0;
              inc_dis += d0;
              pu.emplace_back(u0, d0, i0);
              found = true;
              break;
            }
          }
        }
        if (!found)
        {
          for (int i = 0; i < ptdL[cu].size(); ++i)
          {
            auto u0 = ptdL[cu][i].first.v(), d0 = ptdL[cu][i].first.d(), i0 = ptdL[cu][i].second;
            if (rank[u0] == -1)
            {
              auto dc = retrieve_ptd_dist_on_core(tree[tree[u].rid].nbr[u_itfc], u0);
              if (get<0>(dc) + get<1>(dc) + inc_dis + d0 == tree[u].dis[u_itfc])
              {
                cu = u0;
                inc_dis += d0;
                pu.emplace_back(u0, d0, i0);
                inc_dis += get<0>(dc) + get<1>(dc);
                auto ptmp = retrieve_ptd_path_on_core(u0, tree[tree[u].rid].nbr[u_itfc], dc);
                pu.insert(pu.end(), ptmp.begin(), ptmp.end());
                break;
              }
            }
          }
        }
      }
    }
    return pu;
  }

  deque<tuple<int, int, int>> retrieve_ptd_path_tree2core(unsigned int u, unsigned int u_itfc, int cr, int u2cd)
  {
    deque<tuple<int, int, int>> put2c;
    if (rank[u] != -1)
    {
      put2c = retrieve_ptd_path_on_tree(u, u_itfc);
      auto pucore =
          retrieve_ptd_path_to_p_on_core(tree[tree[u].rid].nbr[u_itfc], {u2cd - tree[u].dis[u_itfc], cr});
      put2c.insert(put2c.end(), pucore.begin(), pucore.end());
    }
    else
    {
      put2c = retrieve_ptd_path_to_p_on_core(u, {u2cd, cr});
    }
    return put2c;
  }

  deque<tuple<int, int, int>> retrieve_ptd_path_on_bp(unsigned int u, int bpri, unsigned long bp_vbin, int bpud)
  {
    deque<tuple<int, int, int>> pu;
    int inc_dis = 0;
    int cu = u;

    while (inc_dis < bpud)
    {
      for (int i = 0; i < ptdL[cu].size(); ++i)
      {
        auto u0 = ptdL[cu][i].first.v(), d0 = ptdL[cu][i].first.d(), i0 = ptdL[cu][i].second;
        if ((((label_bp[u0].bpspt_d[bpri] + inc_dis + d0) == bpud + 1) && label_bp[u0].bpspt_s[bpri].first & bp_vbin) || (((label_bp[u0].bpspt_d[bpri] + inc_dis + d0) == bpud) && (label_bp[u0].bpspt_s[bpri].second & bp_vbin || !bp_vbin)))
        {
          cu = u0;
          inc_dis += d0;
          pu.emplace_back(u0, d0, i0);
          break;
        }
      }
    }
    //        if (inc_dis == bpud - 2) {
    //            pu.emplace_back(rses[bpri][1 + locate_first1_in_bitmap(label_bp[cu].bpspt_s[bpri].first)], 1,
    //                            rses[bpri][1 + locate_first1_in_bitmap(label_bp[cu].bpspt_s[bpri].first)]);
    //        }
    //        if (inc_dis < bpud) {
    //            pu.emplace_back(rses[bpri][1 + locate_first1_in_bitmap(bp_vbin)], 1,
    //                            rses[bpri][1 + locate_first1_in_bitmap(bp_vbin)]);
    //        }
    return pu;
  }

  void unfold_ptd_edge(unsigned int left, unsigned int right, unsigned int middle, vector<unsigned int> &path)
  {
    //        int idxl, idxr;
    auto idxl_it = lower_bound(ptdL[middle].begin(), ptdL[middle].end(), left,
                               [](const pair<vd_bitseg, unsigned int> &a, int b)
                               { return a.first.v() < b; });
    auto idxr_it = lower_bound(ptdL[middle].begin(), ptdL[middle].end(), right,
                               [](const pair<vd_bitseg, unsigned int> &a, int b)
                               { return a.first.v() < b; });
    //        for (int i = 0; i < ptdL[middle].size(); ++i) {//done:search here can be accerlerated with sort
    //            if (ptdL[middle][i].first.v() == left) {
    //                idxl = i;
    //            } else if (ptdL[middle][i].first.v() == right) {
    //                idxr = i;
    //            }
    //        }
    if (idxl_it->first.d() == 1)
    {
      path.push_back(left);
    }
    else if (idxl_it->first.d() == 2)
    {
      path.push_back(left);
      path.push_back(idxl_it->second);
    }
    else if (idxl_it->first.d() > 2)
    {
      unfold_ptd_edge(left, middle, idxl_it->second, path);
    }

    if (idxr_it->first.d() == 1)
    {
      path.push_back(middle);
    }
    else if (idxr_it->first.d() == 2)
    {
      path.push_back(middle);
      path.push_back(idxr_it->second);
    }
    else if (idxr_it->first.d() > 2)
    {
      unfold_ptd_edge(middle, right, idxr_it->second, path);
    }
  }

  vector<unsigned int> unfold_ptd_path(unsigned int u, unsigned int v, const deque<tuple<int, int, int>> &pu,
                                       const deque<tuple<int, int, int>> &pv)
  {
    vector<unsigned int> path;
    unsigned int left = u;
    for (int i = 0; i < pu.size(); ++i)
    {
      if (get<1>(pu[i]) == 1)
      {
        path.push_back(left);
      }
      else if (get<1>(pu[i]) == 2)
      {
        path.push_back(left);
        path.push_back(get<2>(pu[i]));
      }
      else if (get<1>(pu[i]) > 2)
      {
        unfold_ptd_edge(left, get<0>(pu[i]), get<2>(pu[i]), path);
      }
      left = get<0>(pu[i]);
    }

    for (int i = pv.size() - 1; i >= 0; --i)
    {
      if (get<1>(pv[i]) == 1)
      {
        path.push_back(get<0>(pv[i]));
      }
      else if (get<1>(pv[i]) == 2)
      {
        path.push_back(get<0>(pv[i]));
        path.push_back(get<2>(pv[i]));
      }
      else if (get<1>(pv[i]) > 2)
      {
        if (i > 0)
        {
          unfold_ptd_edge(get<0>(pv[i]), get<0>(pv[i - 1]), get<2>(pv[i]), path);
        }
        else
        {
          unfold_ptd_edge(get<0>(pv[i]), v, get<2>(pv[i]), path);
        }
      }
    }
    path.push_back(v);
    return path;
  }

  vector<unsigned int> query_path(int u, int v)
  {
    if (u >= n || v >= n)
      return {};

    int bpri, bpud, bpvd;
    unsigned long bp_vbin;
    int bpd = MAXD;
    { //=================================================bp==========================================
      BPLabel &idx_u = label_bp[u], &idx_v = label_bp[v];
      for (int i = 0; i < params.n_roots; ++i)
      {
        int td = idx_u.bpspt_d[i] + idx_v.bpspt_d[i];
        if (td - 2 < bpd && idx_u.bpspt_s[i].first & idx_v.bpspt_s[i].first)
        {
          bp_vbin = idx_u.bpspt_s[i].first & idx_v.bpspt_s[i].first;
          bp_vbin = bp_vbin - ((bp_vbin - 1) & bp_vbin);
          bpri = i;
          bpd = td - 2;
          bpud = idx_u.bpspt_d[i] - 1;
          bpvd = idx_v.bpspt_d[i] - 1;
        }
        if (td - 1 < bpd)
        {
          auto bp_vbin_tmp =
              (idx_u.bpspt_s[i].first & idx_v.bpspt_s[i].second) | (idx_u.bpspt_s[i].second & idx_v.bpspt_s[i].first);
          if (bp_vbin_tmp)
          {
            bp_vbin = bp_vbin_tmp - ((bp_vbin_tmp - 1) & bp_vbin_tmp);
            bpri = i;
            bpd = td - 1;
            if (bp_vbin & idx_u.bpspt_s[i].first & idx_v.bpspt_s[i].second)
            {
              bpud = idx_u.bpspt_d[i] - 1;
              bpvd = idx_v.bpspt_d[i];
            }
            else
            {
              bpud = idx_u.bpspt_d[i];
              bpvd = idx_v.bpspt_d[i] - 1;
            }
          }
        }
        if (td < bpd)
        {
          bpri = i;
          bpd = td;
          bpud = idx_u.bpspt_d[i];
          bpvd = idx_v.bpspt_d[i];
          bp_vbin = 0l;
        }
      }
    }

    int cored = MAXD;
    unsigned int core_vr, u_itfc, v_itfc, u2cd, v2cd;
    if (!usd_bp[u] && !usd_bp[v])
    {
      //==============================================core=============================================
      static unsigned short nowt = 0;
      static vector<unsigned short> last_t(n_core, 0);
      static vector<pair<int, int>> dis(n_core);
      ++nowt;

      if (nowt == MAXT)
      {
        fill(last_t.begin(), last_t.end(), 0);
        nowt = 1;
      }
      if (rank[u] == -1)
      {
        for (auto l : corL[u])
        {
          last_t[l.v()] = nowt;
          dis[l.v()] = {l.d(), u};
        }
      }
      else
      {
        TreeNode &tu = tree[u];
        TreeNode &r = tree[tu.rid];
        for (size_t i = 0; i < tu.rsize; ++i)
        {
          int x = r.nbr[i], w = tu.dis[i];
          // if(w >= cored) continue;
          for (auto l : corL[x])
          {
            int y = l.v(), nowd = l.d() + w;
            if (nowd >= bpd)
              break;
            if (last_t[y] != nowt)
            {
              last_t[y] = nowt;
              dis[y] = {nowd, i};
            }
            else if (nowd < dis[y].first)
              dis[y] = {nowd, i};
          }
        }
      }
      if (rank[v] == -1)
      {
        for (auto l : corL[v])
        {
          if (last_t[l.v()] == nowt)
          {
            if (cored > l.d() + dis[l.v()].first)
            {
              cored = l.d() + dis[l.v()].first;
              core_vr = l.v();
              v2cd = l.d();
              u_itfc = dis[l.v()].second;
            }
          }
        }
      }
      else
      {
        TreeNode &tv = tree[v];
        TreeNode &r = tree[tv.rid];
        for (size_t i = 0; i < tv.rsize; ++i)
        {
          int x = r.nbr[i], w = tv.dis[i];
          // if(w >= cored) continue;
          for (auto l : corL[x])
          {
            int nowd = l.d() + w;
            if (nowd >= bpd)
              break;
            if (last_t[l.v()] == nowt)
            {
              if (cored > nowd + dis[l.v()].first)
              {
                cored = nowd + dis[l.v()].first;
                core_vr = l.v();
                v2cd = nowd;
                u_itfc = dis[l.v()].second;
                v_itfc = i;
              }
            }
          }
        }
      }
      u2cd = cored - v2cd;
    }

    int treed = MAXD, treeidx;
    if (!usd_bp[u] && !usd_bp[v])
    {
      if (rank[u] >= 0 && rank[v] >= 0 && tree[u].rid == tree[v].rid)
      {
        //====================================tree===================================
        TreeNode &tu = tree[u], &tv = tree[v];
        int len = min(tu.dis.size(), tv.dis.size());
        for (int i = 0; i < len; ++i)
        {
          if (i >= tu.rsize && tu.anc[i - tu.rsize] != tv.anc[i - tv.rsize])
            break;
          if (treed > tu.dis[i] + tv.dis[i])
          {
            treed = tu.dis[i] + tv.dis[i];
            treeidx = i;
          }
        }
      }
    }

    // calculate path on ptd
    deque<tuple<int, int, int>> pu, pv; // v,d,intermediate v
    if (bpd <= cored && bpd <= treed)
    { // answer only by bp
      if (bpd == MAXD)
        return {}; // no path
      pu = retrieve_ptd_path_on_bp(u, bpri, bp_vbin, bpud);
      pv = retrieve_ptd_path_on_bp(v, bpri, bp_vbin, bpvd);
    }
    else if (cored <= treed)
    { // answer only by core
      if (cored == MAXD)
        return {}; // no path
      pu = retrieve_ptd_path_tree2core(u, u_itfc, core_vr, u2cd);
      pv = retrieve_ptd_path_tree2core(v, v_itfc, core_vr, v2cd);
    }
    else
    { // answer only by tree
      if (treed == MAXD)
        return {}; // no path
      pu = retrieve_ptd_path_on_tree(u, treeidx);
      pv = retrieve_ptd_path_on_tree(v, treeidx);
    }
    // unfold ptd path
    return unfold_ptd_path(u, v, pu, pv);
  }

  void query_init()
  {
    n2o_id.resize(n);
    int nid;
    for (int oid = 0; oid < id_map.size(); ++oid)
    {
      if (id_map[oid] < 0)
        continue;
      if (id_map[oid] >= id_map.size())
        continue;
      // oid is not removed
      n2o_id[id_map[oid]] = oid;
    }
  }

  vector<unsigned int> query_path_with_original_id(unsigned int u, unsigned int v)
  {
    if (u == v)
      return {u};
    int fu, fv;
    int eqtype = 0;
    if (id_map[u] < 0)
    {
      fu = -id_map[u] - 1;
      eqtype = 1;
    }
    else if (id_map[u] < id_map.size())
    {
      fu = u;
    }
    else
    {
      fu = id_map[u] - id_map.size();
      eqtype = 2;
    }

    if (id_map[v] < 0)
    {
      fv = -id_map[v] - 1;
      eqtype = 1;
    }
    else if (id_map[v] < id_map.size())
    {
      fv = v;
    }
    else
    {
      fv = id_map[v] - (int)id_map.size();
      eqtype = 2;
    }

    if (fu == fv)
    {
      if (eqtype == 1)
      {
        if (g.size() > id_map[fu])
          return {u, g[id_map[fu]][0], v};
        else
          return {};
      }
      else if (eqtype == 2)
        return {u, v};
    }

    n2o_id[id_map[fu]] = u;
    n2o_id[id_map[fv]] = v;

    auto p = query_path(id_map[fu], id_map[fv]);
    for (auto &nid : p)
    {
      nid = n2o_id[nid];
    }
    return p;
  }

  unsigned diameter_estimation()
  {
    unsigned int max_dis = 0;
    unsigned int furthest_v = 0;
    for (int i = 0; i < 3; ++i)
    {

      unsigned int s = furthest_v;
      unsigned int max_dis_local;
      BFS(s, max_dis_local, furthest_v);
      max_dis = max(max_dis_local, max_dis);
      //      mlog("max_dis{},furthest_v{}", max_dis, furthest_v)
    }
    return max_dis;
  }
  /**
   *
   * @param s: in
   * @param max_dis: out
   * @param furthest_v: out
   * @return dis_map: out
   */
  vector<unsigned char> BFS(unsigned int s, unsigned int &max_dis, unsigned int &furthest_v)
  {
    vector<unsigned char> dis_map(n, MAXD);
    queue<unsigned int> Q;
    Q.push(s);
    dis_map[s] = 0u;
    while (!Q.empty())
    { // Q is not empty
      unsigned int x = Q.front();
      Q.pop();
      unsigned char d = dis_map[x];
      for (auto y : g[x])
      {
        if (dis_map[y] != MAXD)
          continue;
        dis_map[y] = d + 1;
        max_dis = d + 1;
        furthest_v = y;
        Q.push(y);
      }
    }
    return dis_map;
  }

  bool can_update(int v, int dis, vector<char> &nowdis)
  {
    for (auto l : corL[v])
    {
      int d = l.d(), u = l.v();
      if (nowdis[u] >= 0 && nowdis[u] + d <= dis)
        return false;
    }
    return true;
  }

  void print_graph()
  {
    printf("=================print graph=================\n");
    for (size_t i = 0; i < g.size(); ++i)
    {
      cout << i << ": ";
      for (const auto &neb : g[i])
      {
        cout << neb << ", ";
      }
      cout << endl;
    }
    printf("=================print graph done=================\n");
  }
  bool edge_exist(unsigned int u, unsigned int v)
  {
    if (u >= n || v >= n)
      return false;
    if (g[u].size() > g[v].size())
      swap(u, v);
    for (const auto &item : g[u])
    {
      if (item == v)
        return true;
    }
    return false;
  }
};

#endif
