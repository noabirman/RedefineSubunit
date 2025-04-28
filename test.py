import networkx as nx
from merge_graphs import merge_graphs, SubunitInfo, save_subunits_info

graph = nx.Graph()

# Add first high segment (indices 5-10)
graph.add_node(
    "A_1_high",
    data=SubunitInfo(
        name="A_1_high",
        chain_names=["A"],
        start=5,
        end=10,
        sequence="FGHIJK"  # corresponds to indices 5-10 in original sequence
    )
)

# Add second high segment (indices 15-18)
graph.add_node(
    "A_2_high",
    data=SubunitInfo(
        name="A_2_high",
        chain_names=["A"],
        start=15,
        end=18,
        sequence="PQRS"  # corresponds to indices 15-18 in original sequence
    )
)
save_subunits_info(graph, 'name_mapping.json', 'subunits_info.json','new_subunits.json')