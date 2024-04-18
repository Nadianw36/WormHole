import os

scripts_path = "scripts/"

large_graphs = [
    "wikipedia",
    "soc-twitter"
]

mid_graphs = ["pokec",
              "livejournal",
              "orkut",
              "skitter",
              "large-dblp",
              "slashdot",
              "epinions",
              "higgs-twitter",
              "dblp"
              ]

L0_sizes_large = [0.5, 1, 1.5, 2, 4]

L0_sizes_not_large = [2*i for i in range(1, 6)]


def size_to_str(size):
    size_str = str(size)
    if '.' in size_str:
        size_str = size_str.replace(".", "-")
    else:
        size_str += "-0"
    return size_str


def create_MLL_script(graphname):
    MLL_script = f"""#!/bin/bash
#SBATCH -J {graphname}_mll          # Job Name 
#SBATCH -p sched_mit_hill
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH -t 11:59:00   
#SBATCH -o ../output/{graphname}_mll_%j.txt

cd ../MLL
if [ ! -d ../results/{graphname}/MLL ]; then
mkdir ../results/{graphname}/MLL
fi

module rm gcc/4.8.4
module load gcc/11.2.0
./mll index -d {graphname} > ../results/{graphname}/MLL/{graphname}_MLL_setup.txt
./mll query -d {graphname} -q '../inputs/{graphname}_input.txt'
"""
    mll_filepath = scripts_path + f"MLL_{graphname}.sh"
    with open(mll_filepath, 'w') as file:
        file.write(MLL_script)
    os.chmod(mll_filepath, 0o777)


def create_MLL_on_core_script(graphname, L0_sizes):
    for size in L0_sizes:
        size_str = size_to_str(size)
        MLL_script = f"""#!/bin/bash
#SBATCH -J {graphname}_{size}_mll-on-core      # Job Name 
#SBATCH -p sched_mit_hill
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH -t 11:59:00   
#SBATCH -o ../output/{graphname}_mll_%j.txt

cd ../MLL
module rm gcc/4.8.4
module load gcc/11.2.0
./mll index -d {graphname}_seed_{size_str}_core-COO_sanitized > ../results/{graphname}/L0/MLL/{graphname}_{size_str}_MLL_setup.txt

./mll query -d {graphname}_seed_{size_str}_core-COO_sanitized -q ../inputs/MLL-inputs/{graphname}/{graphname}_seed_{size_str}_random-L0.txt
./mll query -d {graphname}_seed_{size_str}_core-COO_sanitized -q ../inputs/MLL-inputs/{graphname}/{graphname}_seed_{size_str}_random-L0-random-L1.txt
"""
        mll_filepath = scripts_path + f"MLL-on-core_{graphname}_{size_str}.sh"
        with open(mll_filepath, 'w') as file:
            file.write(MLL_script)
        os.chmod(mll_filepath, 0o777)


def create_setup_script(graphname):
    setup_script = f"""#!/bin/bash
#SBATCH -J {graphname}_setup        # Job Name
#SBATCH -p sched_mit_hill
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH -t 11:59:00
#SBATCH -o ../output/{graphname}_setup_%j.txt
cd ..
module rm gcc/4.8.4
module load gcc/11.2.0
module load python/3.9.4

exe/generateL0seed {graphname} seed
exe/generateInput {graphname}
"""
    setup_filepath = scripts_path + f"setup_{graphname}.sh"
    with open(setup_filepath, 'w') as file:
        file.write(setup_script)
    os.chmod(setup_filepath, 0o777)


def create_BiBFS_script(graphname):
    BiBFS_script = f"""#!/bin/bash
#SBATCH -J {graphname}_bibfs           # Job Name 
#SBATCH -p sched_mit_hill
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH -t 11:59:00  
#SBATCH -o ../output/{graphname}_bibfs_%j.txt

cd ..
module rm gcc/4.8.4
module load gcc/11.2.0
exe/sample BiBFS {graphname} 0 10000
"""
    bibfs_filepath = scripts_path + f"BiBFS_{graphname}.sh"
    with open(bibfs_filepath, 'w') as file:
        file.write(BiBFS_script)
    os.chmod(bibfs_filepath, 0o777)


def create_L0_scripts(graphname, L0_sizes):
    for size in L0_sizes:
        size_str = size_to_str(size)
        L0_run_script = f"""#!/bin/bash
#SBATCH -J {graphname}_L0_{size_str}           # Job Name 
#SBATCH -p sched_mit_hill
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH -t 11:59:00   
#SBATCH -o ../output/{graphname}_L0_{size_str}_%j.txt

cd ..
module rm gcc/4.8.4
module load gcc/11.2.0
module load python/3.9.4

exe/sample L0-BiBFS {graphname}_seed_{size_str} 0 10000
exe/queryOnCore {graphname}_seed_{size_str} 0 10000

python python/remap_inputs.py {graphname} {size_str}
"""
        L0_filepath = scripts_path + f"L0_{size_str}_{graphname}.sh"
        with open(L0_filepath, 'w') as file:
            file.write(L0_run_script)
        os.chmod(L0_filepath, 0o777)


def create_scripts(graphname, params_dict):
    if params_dict["mll"]:
        create_MLL_script(graphname)

    if params_dict["mll-on-core"]:
        create_MLL_on_core_script(graphname, params["l0_sizes"])

    if params_dict["setup"]:
        create_setup_script(graphname)

    if params_dict["bibfs"]:
        create_BiBFS_script(graphname)

    if params_dict["l0"]:
        create_L0_scripts(graphname, params["l0_sizes"])


# for graphname in graphs:
#     create_scripts(graphname)

def create_submit_script(graph_dict, name):
    sbatch_setup_scripts = '#!/bin/bash'
    for graphname, params in graph_dict.items():
        if params["setup"]:
            sbatch_setup_scripts += f"\nsbatch setup_{graphname}.sh"
    if sbatch_setup_scripts != '#!/bin/bash':
        sbatch_setup_filepath = scripts_path + f"{name}_sbatch_setup.sh"
        with open(sbatch_setup_filepath, 'w') as file:
            file.write(sbatch_setup_scripts)
        os.chmod(sbatch_setup_filepath, 0o777)

    sbatch_run_scripts = '#!/bin/bash'
    for graphname, params in graph_dict.items():
        if params["bibfs"]:
            sbatch_run_scripts += f"\nsbatch BiBFS_{graphname}.sh"
        if params["mll"]:
            sbatch_run_scripts += f"\nsbatch MLL_{graphname}.sh"
        if params["l0"]:
            for size in params["l0_sizes"]:
                size_str = size_to_str(size)
                sbatch_run_scripts += f"\nsbatch L0_{size_str}_{graphname}.sh"
        if params["l0"]:
            for size in params["l0_sizes"]:
                size_str = size_to_str(size)
                sbatch_run_scripts += f"\nsbatch L0_{size_str}_{graphname}.sh"

    sbatch_run_filepath = scripts_path + f"{name}_sbatch_run.sh"
    with open(sbatch_run_filepath, 'w') as file:
        file.write(sbatch_run_scripts)
    os.chmod(sbatch_run_filepath, 0o777)

    sbatch_run_scripts_mll_on_core = '#!/bin/bash'
    for graphname, params in graph_dict.items():
        if params["mll-on-core"]:
            for size in params["l0_sizes"]:
                size_str = size_to_str(size)
                sbatch_run_scripts_mll_on_core += f"\nsbatch MLL-on-core_{graphname}_{size_str}.sh"

    sbatch_run_scripts_mll_on_core_filepath = scripts_path + \
        f"{name}_sbatch_run_mll-on-core.sh"
    with open(sbatch_run_scripts_mll_on_core_filepath, 'w') as file:
        file.write(sbatch_run_scripts_mll_on_core)
    os.chmod(sbatch_run_scripts_mll_on_core_filepath, 0o777)


large_graph_dict = {}
for graph in large_graphs:
    large_graph_dict[graph] = {
        "mll": False,
        "setup": False,
        "bibfs": False,
        "l0": True,
        "mll-on-core": True,
        "l0_sizes": L0_sizes_large
    }

mid_graph_dict = {}
for graph in mid_graphs:
    mid_graph_dict[graph] = {
        "mll": False,
        "setup": False,
        "bibfs": False,
        "l0": True,
        "mll-on-core": True,
        "l0_sizes": L0_sizes_not_large
    }
    if graph in ["dblp", "slashdot", "orkut", "pokec"]:
        mid_graph_dict[graph]["mll-on-core"] = True

for graphname, params in large_graph_dict.items():
    create_scripts(graphname, params)

for graphname, params in mid_graph_dict.items():
    create_scripts(graphname, params)
create_submit_script(mid_graph_dict, "mid-graphs")

for graphname, params in large_graph_dict.items():
    create_scripts(graphname, params)
create_submit_script(large_graph_dict, "large-graphs")
