import os
import pickle
import sys
import networkx as nx


if __name__ == "__main__":
    pkl_path = os.path.abspath(sys.argv[1])
    with open(pkl_path, "rb") as f:
        G = pickle.load(f)

    print("Graph type:", type(G))
    print("Number of nodes:", G.number_of_nodes())
    print("Number of edges:", G.number_of_edges())
    # Count edges where the first letters of the nodes are different
    non_internal_edges = [
        (u, v) for u, v in G.edges()
        if u[0] != v[0]
    ]

    print("Non-internal edges:", len(non_internal_edges))
    print(non_internal_edges[:10])







