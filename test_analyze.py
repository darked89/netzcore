import biana.utilities.graph_utilities as gu
g = gu.create_network_from_sif_file("/data/emre/toy_data/test_interactions_small.sif")
degrees = g.degree(with_labels=True)
node_to_values = gu.get_node_degree_related_values(g, ["v2","v3"])
for v in g.nodes():
    print v, degrees[v], node_to_values[v]
gu.create_R_analyze_network_script(g, ["v2","v3"])

