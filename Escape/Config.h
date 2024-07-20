#ifndef ESCAPE_CONFIG_H_
#define ESCAPE_CONFIG_H_

#include <cstdlib>
#include <cstdio>
#include <set>
#include <map>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace Escape

{
    extern string GRAPH_FOLDER;
    extern string INPUT_FOLDER;
    extern string RESULTS_FOLDER;
    extern string BIN_FOLDER;

    extern string BIBFS_FOLDER;
    extern string MLL_FOLDER;
    extern string PLL_FOLDER;
    extern string MLL_INPUT_FOLDER;
    extern string QBS_FOLDER;
    extern string L0_FOLDER;

    extern string DEGREES_FOLDER;

    extern string BIBFS_RANDOML0_FOLDER;
    extern string BIBFS_RANDOML0FROMRANDOML1_FOLDER;
    extern string BIBFS_HIGHESTDEGL0_FOLDER;

    void checkSetupFor(string graph);

    std::string getSubfolderName(std::string command);

    inline char *getCmdOption(int argc, char **argv, const std::string &option)
    {
        char **itr = std::find(argv, argv + argc, option);
        if (itr != (argv + argc) && ++itr != (argv + argc))
        {
            return *itr;
        }
        return nullptr;
    }
    
    inline bool cmdOptionExists(int argc, char **argv, const std::string &option)
    {
        return std::find(argv, argv + argc, option) != (argv + argc);
    }
}

#endif
