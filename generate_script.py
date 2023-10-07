import os
import numpy as np
import sys
graph = "wikipedia"
ROUNDS = 10000
ROUND_INCREMENT = 200
BiBFS_ROUND_INC = 50
L0_sizes =["0.5", "1.0", "1.5", "2.0", "4.0"]
def make_scripts():
    graphs = ["epinions", "slashdot", "skitter", "higgs-twitter"]
    L0s = ["2.0", "4.0", "6.0", "8.0", "10.0"]
    for g in graphs:
        for L0 in L0s:
            file_path = f'scripts/{g}_{L0}.sh'
            script_content = f'''#!/bin/bash
#SBATCH -J {g}_{L0}               # Job Name 
#SBATCH -p sched_mit_hill
#SBATCH -t 3:00:00   
module load gcc/11.2.0
cd ..
exe/sample L0 {g}_seed_{L0} 0 10000'''
            with open(file_path, 'w') as file:
                file.write(script_content)
            os.chmod(file_path, 0o777)
        file_path = f'scripts/{g}_BiBFS.sh'
        script_content = f'''#!/bin/bash
#SBATCH -J {g}_{L0}               # Job Name 
#SBATCH -p sched_mit_hill
#SBATCH -t 3:00:00   
module load gcc/11.2.0
cd ..
exe/sample BiBFS {g} 0 10000'''
        with open(file_path, 'w') as file:
            file.write(script_content)
        os.chmod(file_path, 0o777)
    sbatch_str = "#!/bin/bash\n"
    sbatch_file = "scripts/sbatch.sh"
    for g in graphs:
        for L0 in L0s:
            sbatch_str+= f"sbatch {g}_{L0}.sh\n"
        sbatch_str+= f"sbatch {g}_BiBFS.sh\n"
    with open(sbatch_file, 'w') as file:
        file.write(sbatch_str)
    os.chmod(sbatch_file, 0o777)
def checkL0_files():
    files = ["queries", "failed", "distances", "runtimes", "distance_to_core",  "runtime_to_core",  "queries_list",  "hops_in_core"]
    L0_sizes =["0.5", "1.0", "1.5", "2.0"]
    for L0_size in L0_sizes:
        print(L0_size)
        for file in files:
            for i in range(0, ROUNDS, ROUND_INCREMENT):
                f = f"results/{i}_{i+ROUND_INCREMENT}_soc-twitter_seed_{L0_size}_L0_{file}.txt"
                if (not os.path.isfile(f)):
                    print(f'{i}')
def checkBiBFS_files():
    files = ["queries", "failed", "distances", "runtimes"]
    for file in files:
        for i in range(0, ROUNDS, 100):
            f = f"results/{i}_{i+100}_soc-twitter_BiBFS_{file}.txt"
            if (not os.path.isfile(f)):
                    print(f'{i}')
def generate_L0_scripts():
    for L0_size in L0_sizes:
        for i in range(0, ROUNDS, ROUND_INCREMENT):
            file_path = f'scripts/runL0_{L0_size}_{i}.sh'
            script_content = f'''#!/bin/bash
#SBATCH -J {graph}_{L0_size}_{i}             # Job Name 
#SBATCH -p sched_mit_hill
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH -t 5:00:00   
#SBATCH -o ../output/%j_{i}_{L0_size}.txt
module load gcc/11.2.0
cd ..
exe/sample L0 {graph}_seed_{L0_size} {i} {i+ROUND_INCREMENT}'''

            with open(file_path, 'w') as file:
                file.write(script_content)
            os.chmod(file_path, 0o777)
def generate_cat_script():
    files = ["queries", "failed", "distances", "runtimes", "distance_to_core",  "runtime_to_core",  "queries_list",  "hops_in_core"]
    bash_script = "#!/bin/bash\n"
    # L0_sizes =["0.5", "1.0", "1.5", "2.0"]
    for L0_size in L0_sizes:
        for file in files:
            command = "cat "
            for i in range(0, ROUNDS, ROUND_INCREMENT):
                command += f"results/{i}_{i+ROUND_INCREMENT}_{graph}_seed_{L0_size}_L0_{file}.txt "
            command += f"> complete/{graph}/0_{ROUNDS}_{graph}_seed_{L0_size}_L0_{file}.txt\n"
            bash_script+=command
    with open("cat_files.sh", 'w') as file:
        file.write(bash_script)
    os.chmod("cat_files.sh", 0o777)

def generate_cat_script_BiBFS():
    files = ["queries", "failed", "distances", "runtimes"]#, "distance_to_core",  "runtime_to_core",  "queries_accumulated",  "hops_in_core"]
    bash_script = "#!/bin/bash\n"
    for file in files:
        command = "cat "
        for i in range(0, ROUNDS, BiBFS_ROUND_INC):
            command += f"results/{i}_{i+BiBFS_ROUND_INC}_soc-twitter_BiBFS_{file}.txt "
        command += f"> complete/0_{ROUNDS}_soc-twitter_BiBFS_{file}.txt\n"
        bash_script+=command
    with open("cat_files_BiBFS.sh", 'w') as file:
        file.write(bash_script)
    os.chmod("cat_files_BiBFS.sh", 0o777)
def create_acc_queries_file_BiBFS():
    interval_size = 100
    visited_stem = "nodesTouched"
    query_stem = "queryTouched"
    visited = np.fromfile(f"bin/100_soc-twitter_BiBFS_{visited_stem}.bin", dtype='bool')
    query = np.fromfile(f"bin/100_soc-twitter_BiBFS_{query_stem}.bin", dtype=np.dtype('uint8')).astype('uint32')
    for i in range(200, ROUNDS+1, interval_size):
        cur_visited = np.fromfile(f"bin/{i}_soc-twitter_BiBFS_{visited_stem}.bin", dtype='bool')
        cur_query = np.fromfile(f"bin/{i}_soc-twitter_BiBFS_{query_stem}.bin", dtype=np.dtype('uint8')).astype('uint32')
        new_visited = cur_visited & np.logical_not(visited)
        # print(i - interval_size + cur_query[new_visited][0:25])
        # print("  ", cur_query[new_visited][0:25], i - interval_size)
        query[new_visited] = (i-interval_size) + cur_query[new_visited]
    del visited
    unique, counts = np.unique(query, return_counts=True)
    del query
    acc_queries = np.zeros(ROUNDS)
    acc_queries[0] = counts[np.where(unique == 0)[0][0]]
    for i in range(1, ROUNDS):
        acc_queries[i] = acc_queries[i-1]
        if(i in unique):
            ind = np.where(unique == i)[0][0]
            acc_queries[i] += counts[ind]
    print(acc_queries[0:50])
    np.savetxt("complete/0_10000_soc-twitter_BiBFS_queries_accumulated.txt", acc_queries )
def create_acc_queries_file():
    visited_stem = "nodesTouched"
    query_stem = "queryTouched"
    for L0_size in L0_sizes:
        visited = np.fromfile(f"bin/{ROUND_INCREMENT}_{graph}_seed_{L0_size}_L0_{visited_stem}.bin", dtype='bool')
        query = np.fromfile(f"bin/{ROUND_INCREMENT}_{graph}_seed_{L0_size}_L0_{query_stem}.bin", dtype=np.dtype('uint8')).astype('uint32')
        for i in range(ROUND_INCREMENT, ROUNDS+1, ROUND_INCREMENT):
            cur_visited = np.fromfile(f"bin/{i}_{graph}_seed_{L0_size}_L0_{visited_stem}.bin", dtype='bool')
            cur_query = np.fromfile(f"bin/{i}_{graph}_seed_{L0_size}_L0_{query_stem}.bin", dtype=np.dtype('uint8')).astype('uint32')
            new_visited = cur_visited & np.logical_not(visited)
            # print(i - interval_size + cur_query[new_visited][0:25])
            # print("  ", cur_query[new_visited][0:25], i - interval_size)
            query[new_visited] = (i-ROUND_INCREMENT) + cur_query[new_visited]
        del visited
        unique, counts = np.unique(query, return_counts=True)
        del query
        acc_queries = np.zeros(ROUNDS)
        acc_queries[0] = counts[np.where(unique == 0)[0][0]]
        for i in range(1, ROUNDS):
            acc_queries[i] = acc_queries[i-1]
            if(i in unique):
                ind = np.where(unique == i)[0][0]
                acc_queries[i] += counts[ind]
        print(acc_queries[0:50])
        np.savetxt(f"complete/{graph}/0_{ROUNDS}_{graph}_seed_{L0_size}_L0_queries_accumulated.txt", acc_queries )
def generate_BiBFS_scripts():
    for i in range(0, ROUNDS, BiBFS_ROUND_INC):
        file_path = f'scripts/runBiBFS_{i}.sh'
        script_content = f'''#!/bin/bash
#SBATCH -J {graph}_B{i}   # Job Name 
#SBATCH -p sched_mit_hill
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH -t 3:00:00   
#SBATCH -o ../output/%j_BiBFS_{i}.txt
module load gcc/11.2.0
cd ..
exe/sample BiBFS {graph} {i} {i+BiBFS_ROUND_INC}'''
        with open(file_path, 'w') as file:
            file.write(script_content)
        os.chmod(file_path, 0o777)
def generate_input_scripts():
    for i in range(0, ROUNDS, 50):
        file_path = f'scripts/generate_input_{i}.sh'
        script_content = f'''#!/bin/bash
#SBATCH -J soc-twitter_input_{i}             # Job Name 
#SBATCH -p sched_mit_hill
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH -t 3:00:00   
#SBATCH -o ../output/%j_input_{i}.txt
module load gcc/11.2.0
cd ..
exe/generateInput soc-twitter {i}'''

        with open(file_path, 'w') as file:
            file.write(script_content)
    os.chmod(file_path, 0o777)
def generate_BiBFS_input_scripts():
    for i in range(0, ROUNDS, 50):
        file_path = f'scripts/generate_BiBFS-input_{i}.sh'
        script_content = f'''#!/bin/bash
#SBATCH -J soc-twitter_input_{i}             # Job Name 
#SBATCH -p sched_mit_hill
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH -t 11:00:00   
#SBATCH -o ../output/%j_input_{i}.txt
module load gcc/11.2.0
cd ..
exe/sample BiBFS-input soc-twitter {i} {i+50}'''

        with open(file_path, 'w') as file:
            file.write(script_content)
    os.chmod(file_path, 0o777)
    
def create_submit_BiBFS_input_script():
    sbatch_str = "#!/bin/bash\n"
    sbatch_file = "scripts/sbatch_BiBFS-input.sh"
    for i in range(0, ROUNDS, 50):
        sbatch_str+= f"sbatch generate_BiBFS-input_{i}.sh\n"
    with open(sbatch_file, 'w') as file:
        file.write(sbatch_str)
    os.chmod(sbatch_file, 0o777)
def create_submit_input_script():
    sbatch_str = "#!/bin/bash\n"
    sbatch_file = "scripts/sbatch_input.sh"
    for i in range(0, ROUNDS, 50):
        sbatch_str+= f"sbatch generate_input_{i}.sh\n"
    with open(sbatch_file, 'w') as file:
        file.write(sbatch_str)
    os.chmod(sbatch_file, 0o777)
def create_submit_BiBFS_script():
    sbatch_str = "#!/bin/bash\n"
    sbatch_file = "scripts/submit_sbatch_BiBFS.sh"
    for i in range(0, ROUNDS, BiBFS_ROUND_INC):
        sbatch_str+= f"sbatch runBiBFS_{i}.sh\n"
    with open(sbatch_file, 'w') as file:
        file.write(sbatch_str)
    os.chmod(sbatch_file, 0o777)
def create_submit_L0_script():
    sbatch_str = ""
    sbatch_file = "scripts/submit_sbatch.sh"
    for L0_size in L0_sizes:
        for i in range(0, ROUNDS, ROUND_INCREMENT):
            sbatch_str+= f"sbatch runL0_{L0_size}_{i}.sh\n"
    with open(sbatch_file, 'w') as file:
        file.write(sbatch_str)
    os.chmod(sbatch_file, 0o777)
def print_methods():
    methods = ["generate_BiBFS_scripts","create_submit_BiBFS_script",  "create_submit_L0_script", "generate_L0_scripts"]
    for method in methods:
        print(method)
def run():
    generate_L0_scripts()
    generate_BiBFS_scripts()
    create_submit_BiBFS_script()
    create_submit_L0_script()
if __name__ == '__main__':
    globals()[sys.argv[1]]()
