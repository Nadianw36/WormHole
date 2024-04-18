#include "MLL_CT.h"

#include <fstream>
using namespace std::chrono;

void parse_args(int argc, char *argv[])
{
  char *arg;

  printf("program args:");
  for (int i = 0; i < argc; ++i)
    printf(" %s", argv[i]);
  printf("\n");

  if (argc > 1 && argv[1][0] != '-')
    params.sub_cmd = argv[1];
  else
    params.sub_cmd = "";

  // dataset
  if ((arg = getCmdOption(argc, argv, "-d")) != nullptr)
    params.dataset = arg;
  else
    params.dataset = "lctd_example";

  // config the degree threshold for tree decomposition in core tree.
  if ((arg = getCmdOption(argc, argv, "-w")) != nullptr)
    params.max_w = stoi(arg);
  else
    params.max_w = 3;

  // specify the number of threads
  if ((arg = getCmdOption(argc, argv, "-t")) != nullptr)
    params.num_threads = stoi(arg);
  else
    params.num_threads = 1;

  if ((arg = getCmdOption(argc, argv, "-qc")) != nullptr)
    params.query_cnt = stoi(arg);
  else
    params.query_cnt = 1000;

  if ((arg = getCmdOption(argc, argv, "-rl")) != nullptr)
    params.report_level = stoi(arg);
  else
    params.report_level = 2;

  if (cmdOptionExists(argc, argv, "-c"))
    params.clear = true;
  else
    params.clear = false;

  if ((arg = getCmdOption(argc, argv, "-nr")) != nullptr)
    params.n_roots = stoi(arg);
  else
    params.n_roots = 0;

  if ((arg = getCmdOption(argc, argv, "-ng")) != nullptr)
    params.num_qgroups = stoi(arg);
  else
    params.num_qgroups = 5;

  if ((arg = getCmdOption(argc, argv, "-gs")) != nullptr)
    params.qgroup_size = stoi(arg);
  else
    params.qgroup_size = 1000;

  if ((arg = getCmdOption(argc, argv, "-sub")) != nullptr)
    params.sub_pct = stoi(arg);
  else
    params.sub_pct = 100;
  if (params.sub_pct != 100)
  {
    params.dataset = params.dataset + "." + to_string(params.sub_pct);
  }

  CMTOPTSET_MACROF("-idx_prefix", params.index_folder, "indexes");
  CMTOPTSET_MACROF("-q", params.query_path, "queries/toyG.queries.txt");
}

int main(int argc, char *argv[])
{
  parse_args(argc, argv);
  if (params.sub_cmd == "index")
  {
    mll_ct mc;
    mc.construct_ctl_ptdl_seperately();
  }
  else if (params.sub_cmd == "query")
  {
    std::string resultsPrefix = "../results/" + params.dataset + "/MLL/" + params.dataset + "_MLL_";
    if (params.dataset.find("seed") != std::string::npos)
    {
      std::string graph_name = params.dataset;

      std::string delimiter = "_";
      std::size_t d_pos = params.dataset.find(delimiter);
      if (d_pos != std::string::npos)
        graph_name = params.dataset.substr(0, d_pos);

      std::string graphL0Name = params.dataset.substr(0, params.dataset.find("_core_COO"));
      std::string input_type = "_random-L0_";
      if (params.query_path.find("random-L1") != std::string::npos)
      {
        input_type = "_random-L0-random-L1_";
      }
      else if (params.query_path.find("highest-deg-L0") != std::string::npos)
      {
        input_type = "_highest-deg-L0_";
      }
      else if (params.query_path.find("all-pairs-L0") != std::string::npos)
      {
        input_type = "_all-pairs-L0_";
      }
      resultsPrefix = "../results/" + graph_name + "/L0/MLL/" + graphL0Name + input_type + "MLL_";
    }

    std::ofstream distancesFile;
    distancesFile.open(resultsPrefix + "distances.txt");
    std::ofstream pathsFile;
    pathsFile.open(resultsPrefix + "paths.txt");
    std::ofstream runtimesFile;
    runtimesFile.open(resultsPrefix + "runtimes.txt");

    mll_ct mc;
    mc.load_label();
    mc.query_init();

    ifstream ifs(params.query_path);
    if (!ifs)
    {
      mlog("open file %s failed.", params.query_path.c_str());
      exit(EXIT_FAILURE);
    }
    string line;
    ifstream::pos_type pos;
    do
    {
      pos = ifs.tellg();
      getline(ifs, line);
    } while (line[0] == '#');
    ifs.seekg(pos);

    high_resolution_clock::time_point start_clock;
    high_resolution_clock::time_point end_clock;

    int s, t;
    long long int totalRuntime = 0;
    int totalQueries = 0;
    while (ifs >> s >> t)
    {
      if (s < -1 || t < -1)
      {
        runtimesFile << -2 << "\n";
        distancesFile << -2 << "\n";
        pathsFile << -2 << "\n";
        continue;
      }
      else if (s < 0 || t < 0)
      {
        runtimesFile << -1 << "\n";
        distancesFile << -1 << "\n";
        pathsFile << -1 << "\n";
        continue;
      }

      start_clock = high_resolution_clock::now();
      auto p = mc.query_path_with_original_id(s, t);
      end_clock = high_resolution_clock::now();

      auto runtime = duration_cast<std::chrono::nanoseconds>(end_clock - start_clock);
      runtimesFile << (long long int)runtime.count() << "\n";
      distancesFile << p.size() << "\n";
      totalRuntime += (long long int)runtime.count();
      totalQueries++;

      for (unsigned int i : p)
      {
        pathsFile << i << ",";
      }
      pathsFile << "\n";
    }
    printf("MLL on %s completed %d inqueries in %lld nanoseconds\n", params.dataset, totalQueries, totalRuntime);
    distancesFile.close();
    pathsFile.close();
    runtimesFile.close();
  }

  return 0;
}