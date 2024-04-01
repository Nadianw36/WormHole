#include "Escape/Config.h"

// global folders
string Escape::BIN_FOLDER = "bin/";
string Escape::RESULTS_FOLDER = "results/";
string Escape::INPUT_FOLDER = "inputs/";
string Escape::GRAPH_FOLDER = "graphs/";
string Escape::BIBFS_FOLDER = "BiBFS/";
string Escape::QBS_FOLDER = "QbS/";
string Escape::L0_FOLDER = "L0/";
string Escape::DEGREES_FOLDER = "core-stats/";

/** intended folder structure
 * /bin
 *      /graph_name
 *          /L0
 * /input
 * /graph
 * /results
 *      /graph_name
 *          /BiBFS
 *          /QbS
 *          /L0
 *              /BiBFS
 *              /QbS
 *
 */

void Escape::checkSetupFor(string graph_name)
{
    create_directory(INPUT_FOLDER);
    create_directory(GRAPH_FOLDER);

    string graph_name_folder = graph_name + "/";
    create_directories(BIN_FOLDER + graph_name_folder + "L0");

    create_directories(RESULTS_FOLDER + graph_name_folder);

    string resultsSubDirectories[3] = {
        BIBFS_FOLDER,
        QBS_FOLDER,
        L0_FOLDER};
    string resultsL0SubDirectories[3] = {
        BIBFS_FOLDER,
        QBS_FOLDER,
        DEGREES_FOLDER};

    for (int i = 0; i < 3; i++)
        create_directory(RESULTS_FOLDER + graph_name_folder + resultsSubDirectories[i]);

    for (int i = 0; i < 3; i++)
        create_directory(RESULTS_FOLDER + graph_name_folder + L0_FOLDER + resultsL0SubDirectories[i]);
}

std::string Escape::getSubfolderName(std::string command)
{
    if (!command.compare("L0-BiBFS"))
        return L0_FOLDER + BIBFS_FOLDER;
    if (!command.compare("BiBFS"))
        return BIBFS_FOLDER;
    return "";
}