import numpy as np
import sys

if __name__ == "__main__":
    args = sys.argv[1:]
    graphname = args[0]
    GRAPHFILE = "graphs/" + graphname + ".txt"
    MAPPING = "graphs/" + graphname + "_mapping.txt"
    C00FILE = "graphs/" + graphname + "_sanitized.edges"
    edgelist = np.loadtxt(GRAPHFILE)
    print(edgelist.shape)
    inverse = np.zeros(int(edgelist[0, 0]))
    mappings = np.ones(int(edgelist[1:].max()) + 1)*-1
    sanitized = []
    sanitized.append([edgelist[0, 0], edgelist[0, 1]])
    counter = 0
    for i in range(1, len(edgelist)):
        e = edgelist[i]
        u = int(e[0])
        v = int(e[1])
        if mappings[u] == -1:
            mappings[u] = counter
            inverse[counter] = u
            counter += 1
        if mappings[v] == -1:
            mappings[v] = counter
            inverse[counter] = v
            counter += 1
        x = int(mappings[u])
        y = int(mappings[v])
        sanitized.append([x, y])

    coo = np.array(sanitized)
    coo_u = np.sort(coo[1:], axis=1)
    coo_u = np.unique(coo_u, axis=0)
    row_indices = np.where(coo_u[:, 0] == coo_u[:, 1])[0]
    if (len(row_indices)):
        coo_u = np.delete(coo_u, row_indices, axis=0)
    edges = len(coo_u)
    print(edges, np.max(coo_u)+1)
    vertices = np.max(coo_u) + 1
    coo = np.vstack(([vertices, edges], coo_u))
    np.savetxt(MAPPING, inverse, fmt='%i', delimiter='\t')
    np.savetxt(C00FILE, coo, fmt='%i', delimiter='\t')
