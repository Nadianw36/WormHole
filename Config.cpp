#include "Escape/Config.h"

// global folders
string Escape::BIN_FOLDER = "bin/";
string Escape::RESULTS_FOLDER = "results/";
string Escape::INPUT_FOLDER = "inputs/";
string Escape::GRAPH_FOLDER = "graphs/";
string Escape::BIBFS_FOLDER = "BiBFS/";
string Escape::MLL_FOLDER = "MLL/";
string Escape::PLL_FOLDER = "PLL/";
string Escape::MLL_INPUT_FOLDER = "MLL-inputs/";
string Escape::QBS_FOLDER = "QbS/";
string Escape::L0_FOLDER = "L0/";
string Escape::DEGREES_FOLDER = "core-stats/";

string Escape::BIBFS_RANDOML0_FOLDER = "randomL0/";
string Escape::BIBFS_RANDOML0FROMRANDOML1_FOLDER = "randomL0FromRandomL1/";
string Escape::BIBFS_HIGHESTDEGL0_FOLDER = "highestDegL0/";

/** intended folder structure
 * /bin
 *      /graph_name
 *          /L0
 * /input
 *  /MLL-input
 *      /graphname
 * /graph
 * /results
 *      /graph_name
 *          /BiBFS
 *          /QbS
 *          /L0
 *              /BiBFS
 *                  /manysL0s
 *                  /randomL0
 *                  /randomL0FromRandomL1
 *              /QbS
 *
 */

void Escape::checkL0SetupFor(string graph_name)
{
    create_directory(INPUT_FOLDER);
    create_directory(INPUT_FOLDER + MLL_INPUT_FOLDER);
    create_directory(BIN_FOLDER);
    create_directory(GRAPH_FOLDER);
    create_directory(RESULTS_FOLDER);

    string graph_name_folder = graph_name + "/";

    create_directories(BIN_FOLDER + graph_name_folder + "L0");
    create_directories(RESULTS_FOLDER + graph_name_folder);
    create_directory(INPUT_FOLDER + MLL_INPUT_FOLDER + graph_name_folder);

    string resultsSubDirectories[5] = {
        BIBFS_FOLDER,
        QBS_FOLDER,
        L0_FOLDER,
        MLL_FOLDER,
        PLL_FOLDER};

    string resultsL0SubDirectories[4] = {
        BIBFS_FOLDER,
        QBS_FOLDER,
        DEGREES_FOLDER,
        MLL_FOLDER};

    string resultsL0BiBFSSubDirectories[3] = {
        BIBFS_RANDOML0_FOLDER,
        BIBFS_RANDOML0FROMRANDOML1_FOLDER,
        BIBFS_HIGHESTDEGL0_FOLDER};

    for (int i = 0; i < 5; i++)
        create_directory(RESULTS_FOLDER + graph_name_folder + resultsSubDirectories[i]);

    for (int i = 0; i < 4; i++)
        create_directory(RESULTS_FOLDER + graph_name_folder + L0_FOLDER + resultsL0SubDirectories[i]);

    for (int i = 0; i < 3; i++)
        create_directory(RESULTS_FOLDER + graph_name_folder + L0_FOLDER + BIBFS_FOLDER + resultsL0BiBFSSubDirectories[i]);
}

std::string Escape::getSubfolderName(std::string command)
{
    if (!command.compare("L0-BiBFS"))
        return L0_FOLDER + BIBFS_FOLDER;
    if (!command.compare("BiBFS"))
        return BIBFS_FOLDER;
    if (!command.compare("PLL"))
        return PLL_FOLDER;
    if (!command.compare("L0-BiBFS-RandomL0"))
        return L0_FOLDER + BIBFS_FOLDER + BIBFS_RANDOML0_FOLDER;
    if (!command.compare("L0-BiBFS-RandomL0FromRandomL1"))
        return L0_FOLDER + BIBFS_FOLDER + BIBFS_RANDOML0FROMRANDOML1_FOLDER;
    if (!command.compare("L0-BiBFS-HighestDegL0"))
        return L0_FOLDER + BIBFS_FOLDER + BIBFS_HIGHESTDEGL0_FOLDER;
    return "";
}