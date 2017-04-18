#!/usr/bin/env python

import sys
sys.path.insert('/sbi/users/interchange/SRC/toolbox', 0)
import guild_utilities 

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

    # EXCERCISE - part A:
    # fill the code to get nodes, seed_to_score, edges and edge_to_score variables below
    # hint: have a look at the code of top_scoring and top_network 

    guild_utilities.generate_cross_validation_node_score_files(nodes, seed_to_score, node_scores_file, xval = 3, default_score = 0.01, replicable = 123)

    guild_utilities.generate_cross_validation_edge_score_as_node_score_files(edges, seed_to_score, edge_to_score, edge_scores_file, xval = 3, default_score = 0.01, replicable = 123)

    # EXCERCISE - part B:
    # modify to code below to run GUILD scoring, e.g., NetScore algorithm on these cross validation files
    # hint: run_scoring(scoring_folder, executable_path, scoring_type="netscore", parameters={"n_iteration":2, "n_repetition":3}, qname=None, calculate_pvalue=False, ...)
    return

if __name__ == "__main__":
    main()

