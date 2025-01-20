from collections import defaultdict


def merge_graphs(graphs):
    # Helper to determine if two vertices overlap
    def overlap(v1, v2):
        return v1["chain"] == v2["chain"] and not (v1["end"] < v2["start"] or v2["end"] < v1["start"])

    # Step 1: Merge vertices
    merged_vertices = {}
    vertex_groups = defaultdict(list)

    # Create a mapping of vertices to groups by overlap
    for graph in graphs:
        for vertex in graph[0]:  # Graph vertices
            merged = False
            for group_key in list(vertex_groups.keys()):
                if any(overlap(vertex, v) for v in vertex_groups[group_key]):
                    vertex_groups[group_key].append(vertex)
                    merged = True
                    break
            if not merged:
                vertex_groups[vertex["name"]].append(vertex)

    # Create merged vertices
    for group_key, vertices in vertex_groups.items():
        names = "_".join(v["name"] for v in vertices)
        chain = vertices[0]["chain"]
        start = min(v["start"] for v in vertices)
        end = max(v["end"] for v in vertices)
        merged_vertices[names] = {"name": names, "chain": chain, "start": start, "end": end}

    # Step 2: Update edges
    merged_edges = set()
    vertex_mapping = {}
    for names, merged_vertex in merged_vertices.items():
        for original_vertex in vertex_groups[names.split("_")[0]]:
            vertex_mapping[original_vertex["name"]] = merged_vertex["name"]

    for graph in graphs:
        for edge in graph[1]:  # Graph edges
            v1, v2 = edge
            merged_edges.add((vertex_mapping[v1], vertex_mapping[v2]))

    # Step 3: Construct the merged graph
    merged_graph = (list(merged_vertices.values()), list(merged_edges))
    return merged_graph


# Example usage
g1 = (
    [
        {"name": "a1", "chain": "A", "start": 20, "end": 50},
        {"name": "b1", "chain": "B", "start": 30, "end": 100},
    ],
    [("a1", "b1")],
)

g2 = (
    [
        {"name": "b2", "chain": "B", "start": 10, "end": 80},
        {"name": "c1", "chain": "C", "start": 20, "end": 90},
    ],
    [("b2", "c1")],
)

g3 = (
    [
        {"name": "a2", "chain": "A", "start": 60, "end": 70},
        {"name": "c2", "chain": "C", "start": 80, "end": 150},
    ],
    [("a2", "c2")],
)

graphs = [g1, g2, g3]
merged_graph = merge_graphs(graphs)

# Output merged graph
print("Merged Vertices:", merged_graph[0])
print("Merged Edges:", merged_graph[1])
