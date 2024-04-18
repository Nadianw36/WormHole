import numpy as np
import io
import sys

if __name__ == "__main__":
    args = sys.argv[1:]
    graphname = args[0]
    L0_size = args[1]
    start_inquery = args[2]
    end_inquery = args[3]

    sanitized_COO_mapping = f"graphs/{graphname}_seed_{L0_size}_core-COO_mapping.txt"
    mapping = np.loadtxt(sanitized_COO_mapping)

    samples = ["random-L0", "random-L0-random-L1",
               "highest-deg-L0", "all-pairs-L0"]
    for sample in samples:

        filename = f"results/{graphname}/L0/BiBFS/{start_inquery}_{end_inquery}_{graphname}_seed_{L0_size}_L0-BiBFS_core-queries_{sample}.txt"
        print(filename)
        filename_output = f"inputs/MLL-inputs/{graphname}/{graphname}_seed_{L0_size}_{sample}.txt"

        queries = np.loadtxt(filename)
        remapped_queries = np.zeros_like(queries)

        for q in range(queries.shape[0]):
            query = queries[q]
            if (query[0] == -2 or query[1] == -2):
                remapped_queries[q][0] = -2
                remapped_queries[q][1] = -2
            elif (query[0] == -1 or query[1] == -1):
                remapped_queries[q][0] = -1
                remapped_queries[q][1] = -1
            else:
                remapped_queries[q][0] = np.where(mapping == query[0])[0][0]
                remapped_queries[q][1] = np.where(mapping == query[1])[0][0]

        np.savetxt(filename_output, remapped_queries, fmt='%i',)
