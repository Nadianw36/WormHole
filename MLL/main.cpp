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

    std::ofstream distancesFile;
    distancesFile.open(resultsPrefix + "distances.txt");
    std::ofstream pathsFile;
    distancesFile.open(resultsPrefix + "paths.txt");
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
    while (ifs >> s >> t)
    {
      start_clock = high_resolution_clock::now();
      auto p = mc.query_path_with_original_id(s, t);
      end_clock = high_resolution_clock::now();

      auto runtime = duration_cast<std::chrono::microseconds>(end_clock - start_clock);
      runtimesFile << (double)runtime.count() << "\n";
      distancesFile << p.size() << "\n";

      for (unsigned int i : p)
      {
        pathsFile << i << ",";
      }
      pathsFile << "\n";
    }

    distancesFile.close();
    pathsFile.close();
    runtimesFile.close();
  }

  return 0;
}