#########################################################################
# Main workflow for scoring 
#
# eg 25/06/2009
#########################################################################

from biana.utilities import biana_output_converter as biana_output_converter 
from biana.utilities import graph_utilities as network_utilities
import prepare_data
import score_network

#arff_file = "../data/aneurist.arff"
network_file = "/home/emre/arastirma/data/ppi_exp/human/human_interactome_biana/human.sif"
network_file_method_attribute = "/home/emre/arastirma/data/ppi_exp/human/human_interactome_biana/human_method_id.eda"
network_file_source_attribute = "/home/emre/arastirma/data/ppi_exp/human/human_interactome_biana/human_source.eda"
network_file_pubmed_attribute = "/home/emre/arastirma/data/ppi_exp/human/human_interactome_biana/human_pubmed.eda"
node_file = "/home/emre/arastirma/data/ppi_exp/human/human_interactome_biana/human_nodes.tsv"
node_file_score = "/home/emre/arastirma/data/aneurist/aneurist_Jul_05/aneursym_scores.txt"
network_file_filtered = "/home/emre/arastirma/data/ppi_exp/human/human_interactome_biana/human_no_tap.sif"
node_file_netzcore = "../data/input/node_scores.txt"
node_file_netzcore_prefix = "../data/input/node_scores"
edge_file_netzcore = "../data/input/edge_weights.txt"
edge_file_netzcore_combined = "../data/input/edge_weights_combined.txt"
edge_file_netzcore_relevance = "../data/input/edge_relevance_scores.txt"
edge_file_netzcore_relevance_prefix = "../data/input/edge_relevance_scores"
arff_file = "../data/arff/network_metrics.arff"
arff_file_prefix = "../data/arff/network_metrics"
result_file_prefix = "../data/output/fold"

def main():
    prepare(create_method_filtered_files = False, filter_hubs = False)
    #score()
    #test_timing_score_one_fold()
    return

def test_timing_score_one_fold():
    from time import clock
    t1 = clock()
    print score_network.test_run(edge_file_netzcore)
    t2 = clock()
    print "Time: ", t2-t1
    return

def score():
    #g = score_network.create_network_from_weight_and_score_files(edge_file_weights = edge_file_netzcore, edge_file_scores = edge_file_netzcore_relevance)
    #node_to_score = score_network.score_by_shortest_paths(g) 
    #scores = node_to_score.values()
    #scores.sort()
    #print scores[-20:]
    score_network.run_and_assess_performance_of_folds(10, edge_file_netzcore, edge_file_netzcore_relevance_prefix, node_file_netzcore, node_file_netzcore_prefix, result_file_prefix)
    return

def prepare(create_method_filtered_files = False, filter_hubs = False):

    #prepare_data.sample_network_preserving_topology(edge_file_netzcore_relevance[:-4] + ".sif", 100, "../data/sampled_graphs/sampled_graph.txt.")
    prepare_data.sample_network_preserving_topology("../data/toy_data/test_interactions_small.sif", 4, "../data/sampled_graphs/sampled_graph.txt.")

    return

    # Filter by detection method (non-tap interactions)
    if create_method_filtered_files:
	prepare_data.create_method_filtered_network_files(network_file_method_attribute)

    # Load network
    #create_network_from_sif_file(network_file)
    g = network_utilities.create_network_from_sif_file(network_file_filtered)

    network_utilities.analyze_network(g)

    # Filter network by degree and get largest connected component
    if filter_hubs:
	g = network_utilities.filter_network(g, 90) 
    else:
	g = network_utilities.filter_network(g, 175)

    # Get degrees of highly connected nodes
    network_utilities.analyze_network(g)

    # Check how many of the seed genes we cover in the network
    seeds = prepare_data.check_seed_coverage_and_get_seed_nodes(g, node_file, node_file_score)

    print seeds

    return

    # Get scores of each node
    node_to_score, seeds_all = prepare_data.map_scores_to_biana_nodes(node_file, node_file_score)

    # Assign node scores either as node relevance score or as edge relevance score
    prepare_data.assign_node_scores(g = g, node_to_score = node_to_score, out_file = node_file_netzcore, as_edge_relevance_score = False, seeds_ignored=None)
    prepare_data.assign_node_scores(g = g, node_to_score = node_to_score, out_file = edge_file_netzcore_relevance, as_edge_relevance_score = True, seeds_ignored=None)
    
    return
    
    # Create arff file for all network
    network_utilities.create_arff_file_with_network_metrics(g, node_to_score, seeds, arff_file)

    # Generate cross validation node files
    prepare_data.generate_cross_validation_node_files(g = g, node_to_score = node_to_score, seeds = seeds, node_file_netzcore_prefix = node_file_netzcore_prefix, edge_file_netzcore_relevance_prefix =  edge_file_netzcore_relevance_prefix, arff_file_prefix = arff_file_prefix)

    # Assign edge scores
    prepare_data.assign_edge_reliability_scores(g = g, network_file_method_attribute = network_file_method_attribute, network_file_source_attribute = network_file_source_attribute, network_file_pubmed_attribute = network_file_pubmed_attribute, out_file = edge_file_netzcore)

    #create_arff_file(node_score_file, network_file, arff_file)

    #prepare_data.generate_cross_validation_test_test_arff_files(10, arff_file, arff_file+".id_score_filtered.arff", arff_file_prefix)
    return

if __name__ == "__main__":
    main()

