#!/usr/bin/env python

import sys
sys.path.insert('/sbi/users/interchange/SRC/toolbox', 0)
import guild_utilities

def main():
    # Set the seed and network files
    data_dir = "/sbi/users/interchange/DATA/guild_tutorial/"
    seed_file = data_dir + "seeds.txt"
    network_file = data_dir + "interactions.sif"
    scoring_folder = data_dir + "test/"
    executable_path = "/sbi/users/interchange/SRC/guild/scoreN"

    # Create input files for scoring
    guild_utilities.prepare_scoring(network_file, seed_file, scoring_folder, non_seed_score=0.01, seed_score=1.0, edge_score=1.0, n_sample=100, delim=" ")

    # Run GUILD and create output files, the case for Netcombo
    guild_utilities.run_scoring(scoring_folder, executable_path, scoring_type="netcombo")
    #run_scoring(scoring_folder, executable_path, scoring_type="netzcore", parameters={"n_iteration":5, "n_sample":100, "sampling_prefix":scoring_folder+"sampled_graph."}, qname=None)
    return

if __name__ == "__main__":
    main()

