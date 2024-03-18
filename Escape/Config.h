#ifndef ESCAPE_CONFIG_H_
#define ESCAPE_CONFIG_H_

#include <cstdlib>
#include <cstdio>
#include <set>
#include <map>
#include <string>
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
    extern string QBS_FOLDER;
    extern string L0_FOLDER;
    extern string DEGREES_FOLDER;

    void checkSetupFor(string graph);

    std::string getSubfolderName(std::string command);
}

#endif
