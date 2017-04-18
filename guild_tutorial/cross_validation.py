#!/usr/bin/env python

import sys
sys.path.insert('/sbi/users/interchange/SRC/toolbox', 0)
import guild_utilities, network_utilities 

def main():
    # Set the seed and network files
    data_dir = "../../DATA/guild_tutorial/"
    seed_file = data_dir + "seeds.txt"
    network_file = data_dir + "interactions.sif"
    scoring_folder = data_dir + "test/"
    executable_path = "../guild/scoreN"

    # Create input files for scoring
    guild_utilities.prepare_scoring(network_file, seed_file, scoring_folder, non_seed_score=0.01, seed_score=1.0, edge_score=1.0, n_sample=100, delim=" ")
    
    # Generate cross validation files
    node_scores_file = scoring_folder + "node_scores.sif"
    edge_scores_file = scoring_folder + "edge_scores_netshort.sif"

    # fill the code to get nodes, seed_to_score, edges and edge_to_score variables below
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)
    seeds = guild_utilities.get_nodes(seed_file)
    nodes = g.nodes()
    edges = g.edges()
    seed_to_score = dict([(node, 1) for node in seeds])
    edge_to_score = dict([((u,v), 1) for u,v in edges])

    guild_utilities.generate_cross_validation_node_score_files(nodes, seed_to_score, node_scores_file, xval = 3, default_score = 0.01, replicable = 123)

    guild_utilities.generate_cross_validation_edge_score_as_node_score_files(edges, seed_to_score, edge_to_score, edge_scores_file, xval = 3, default_score = 0.01, replicable = 123)

    # Run NetScore on these cross validation files
    guild_utilities.run_scoring(scoring_folder, executable_path, scoring_type="netscore", parameters={"n_iteration":2, "n_repetition":3}, qname=None, calculate_pvalue=True, xval=3)
    return

if __name__ == "__main__":
    main()

