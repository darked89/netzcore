#!/usr/bin/env python

import sys
sys.path.insert('/sbi/users/interchange/SRC/toolbox', 0)
import guild_utilities, functional_enrichment 

def main():
    """
	Get nodes that are top scoring w.r.t. GUILD scores.
	Assumes that GUILD scores have been calculated already (i.e. python hello_world.py).
    """

    # Set the seed and network files
    data_dir = "../../DATA/guild_tutorial/"
    seed_file = data_dir + "seeds.txt"
    network_file = data_dir + "interactions.sif"

    scoring_folder = data_dir + "test/"
    pvalue_file = scoring_folder + "output_scores.sif.netcombo.pval"
    enrichment_file = scoring_folder + "enrichment.txt"

    # Get GUILD scores
    node_to_vals = guild_utilities.get_values_from_pvalue_file(pvalue_file)

    # Get seeds
    seeds = guild_utilities.get_nodes(seed_file)

    # Get all nodes
    all_nodes = node_to_vals.keys()

    # Get top scoring, i.e. nodes that have p-value <= 0.05
    top_nodes = set()
    for node, vals in node_to_vals.iteritems():
	score, pval = vals
	if pval <= 0.05:
	    top_nodes.add(node)

    print "%d genes in network, %d in top scoring genes and %d of these top scoring are seeds" % (len(all_nodes), len(top_nodes), len(top_nodes & seeds))

    # Get functional enrichment for these nodes
    f = open(enrichment_file, 'w')
    functional_enrichment.check_functional_enrichment(list(top_nodes), all_nodes, "geneid", f.write, species = "Drosophila melanogaster") 
    f.close()
    return

if __name__ == "__main__":
    main()

