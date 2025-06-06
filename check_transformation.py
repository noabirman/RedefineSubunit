import networkx as nx
import os
import sys
import json

def validate_connected_components(transformations_path, subunit_names):
    G = nx.Graph()
    G.add_nodes_from(subunit_names)

    for fname in os.listdir(transformations_path):
        if fname.endswith(".json"):
            parts = fname.replace(".json", "").split("_")
            if len(parts) == 2:
                G.add_edge(parts[0], parts[1])

    if not nx.is_connected(G):
        print("Subunit transformation graph is not connected!")
        for component in nx.connected_components(G):
            print("Component:", component)

if __name__ == '__main__':
    if __name__ == '__main__':
        if len(sys.argv) == 3:
            transformations_path = os.path.abspath(sys.argv[1])
            subunits_info_path = os.path.abspath(sys.argv[2])

            with open(subunits_info_path, 'r') as f:
                subunits_info = json.load(f)

            subunit_names = [s["name"] for s in subunits_info.values()]
            validate_connected_components(transformations_path, subunit_names)
        else:
            print("usage: <script> transformations_path subunits_info")




