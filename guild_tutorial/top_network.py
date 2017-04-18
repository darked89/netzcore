#!/usr/bin/env python

import sys
sys.path.insert('/sbi/users/interchange/SRC/toolbox', 0)
import guild_utilities, network_utilities, stat_utilities 

def main():
    """
	Get nodes that are top scoring w.r.t. GUILD scores.
	Assumes that GUILD scores have been calculated already (i.e. python hello_world.py).
    """

    # Set the seed and network files
    data_dir = "../../DATA/guild_tutorial/"
    seed_file = data_dir + "seeds.txt"
    network_file = data_dir + "interactions.sif"
    enrichment_file = data_dir + "enrichment.txt"

    scoring_folder = data_dir + "test/"
    pvalue_file = scoring_folder + "output_scores.sif.netcombo.pval"
    subnetwork_file = scoring_folder + "subnetwork.sif"

    # Get GUILD scores
    node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)

    # Get top scoring, i.e. nodes that have p-value <= 0.05
    top_nodes = set()
    for node, vals in node_to_vals.iteritems():
	score, pval = vals
	if pval <= 0.05:
	    top_nodes.add(node)

    # Load interaction network
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)

    # Get subnetwork induced by top scoring nodes
    g_sub = network_utilities.get_subgraph(g, top_nodes)

    # Output subnetwork along with the inverted p-value scores (z-scores) calculated for edges
    f = open(subnetwork_file, 'w')
    for u, v in g_sub.edges():
	zscore_u = stat_utilities.convert_p_values_to_z_scores([node_to_vals[u][1]])[0]
	zscore_v = stat_utilities.convert_p_values_to_z_scores([node_to_vals[v][1]])[0]
	score = (zscore_u + zscore_v) / 2
	f.write("%s\t%f\t%s\n" % (u, score, v))
    f.close()
    return

if __name__ == "__main__":
    main()

