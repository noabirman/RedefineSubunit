import json
import os
import pickle
import sys

import create_graph_from_af_model

from vizualization_plots import show_subunit_circle_graph

uniprot_to_protein = {
    "P06239": "Lck",
    "Q04759": "PKCtheta",
    "P08575": "CD45",
    "P43403": "Zap70",
    "O43561": "LAT",
    "Q13094": "SLP76 (LCP2)",
    "O75791": "Gads",
    "P62993": "Grb2",
    "Q07889": "SOS1",
    "Q07890": "SOS2",
    "P22681": "Cbl",
    "P19174": "Plc-gamma1",
    "P41240": "Csk",
    "P05107": "LFA-1",
    "P15498": "VAV1",
    "P16333": "Nck",
    "Q08881": "ITK",
}

if __name__ == "__main__":
    # test the function
    folder = os.path.abspath(sys.argv[1])
    parent_dir = os.path.dirname(folder)
    chain_mapping_path = os.path.join(parent_dir, 'chain_id_mapping.json')

    subunits_info_path = os.path.join(folder, "subunits_info.json")

    with open(subunits_info_path, 'r') as f:
        subunits_info = json.load(f)


    with open(os.path.join(folder,"graph.pkl"), "rb") as f:
        graph = pickle.load(f)


    show_subunit_circle_graph(graph, subunits_info, folder,chain_mapping_path, uniprot_to_protein)

