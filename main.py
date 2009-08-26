#########################################################################
# Main workflow for scoring 
#
# eg 25/06/2009
#########################################################################

from biana.utilities import biana_output_converter as biana_output_converter 
from biana.utilities import graph_utilities as network_utilities
import prepare_data
import score_network
import analyze_results
import os, os.path

# Scoring related parameters 
DEFAULT_NON_SEED_SCORE = 0

# Project directory structure
data_dir = ".." + os.sep + "data" + os.sep
input_dir = data_dir + os.sep + "input" + os.sep
output_dir = data_dir + os.sep + "output" + os.sep
sampling_dir = data_dir + "sampled_grahps" + os.sep

# BIANA node & network files
node_file_prefix = data_dir + "human_interactome_biana" + os.sep + "human_nodes"
network_file_prefix = data_dir + "human_interactome_biana" + os.sep + "human_network"
network_file_filtered_by_method = network_file_prefix + "_no_tap.sif"
network_file_filtered_by_degree = network_file_filtered_by_method[:-4] + "no_hub.sif"

# Disease association files
aneurysm_scores_file = data_dir + "aneurist" + os.sep + "aneurysm_associated_genes.txt"
aneurysm_scores_all_equal_file = data_dir + "aneurist" + os.sep + "aneurysm_associated_genes_all_equal.txt"
association_scores_file = aneurysm_scores_all_equal_file

# Node & edge score files
node_scores_file = input_dir + "human_node_scores.sif"
edge_scores_file = input_dir + "human_filtered.sif"
sampled_file_prefix = sampling_dir + "sampled_graph"

# Filtered files
network_file_method_attribute = network_file_prefix + "_method_id.eda"
network_file_source_attribute = network_file_prefix + "_source.eda"
network_file_pubmed_attribute = network_file_prefix + "_pubmed.eda"

# Old scoring files
node_file_netzcore = "../data/input/node_scores.txt"
node_file_netzcore_prefix = "../data/input/node_scores"
edge_file_netzcore = "../data/input/edge_weights.txt"
edge_file_netzcore_combined = "../data/input/edge_weights_combined.txt"
edge_file_netzcore_relevance = "../data/input/edge_relevance_scores.txt"
edge_file_netzcore_relevance_prefix = "../data/input/edge_relevance_scores"
arff_file = "../data/arff/network_metrics.arff"
arff_file_prefix = "../data/arff/network_metrics"
result_file_prefix = "../data/output/fold"

# Resulting score files
netshort_results = output_dir + "netshort_onlyEdgeRelevance_noNodeInitialScoreAccumulation.txt"
netrank_results = output_dir + "netrank_onlyEdgeRelevance_noNodeInitialScoreAccumulation.txt"
netscore_results_1 = output_dir + "netscore_onlyEdgeRelevance_noNodeInitialScoreAccumulation_i1.txt"
netscore_results_3 = output_dir + "netscore_onlyEdgeRelevance_noNodeInitialScoreAccumulation_i3.txt"
netscore_results_no_1 = output_dir + "netscore_noEdgeRelevance_noNodeInitialScoreAccumulation_i1.txt"
netscore_results_no_3 = output_dir + "netscore_noEdgeRelevance_noNodeInitialScoreAccumulation_i3.txt"
netzcore_results_1 = output_dir + "netzcore_onlyEdgeRelevance_noNodeInitialScoreAccumulation_i1.txt"
netzcore_results_3 = output_dir + "netzcore_onlyEdgeRelevance_noNodeInitialScoreAccumulation_i3.txt"

def main():
    prepare(filter_hubs = False, sample_graphs = False)
    #analyze()
    #score()
    #test_timing_score_one_fold()
    return


def analyze():
    result_files = [ netshort_results, netrank_results, netscore_results_1, netscore_results_3,  netscore_results_no_1, netscore_results_no_3, netzcore_results_1, netzcore_results_3 ]
    for percentage in (10, 25, 50):
	print "---- %s:" % percentage
	for f in result_files:
	    print f
	    print analyze_results.calculate_seed_coverage_at_given_percentage(node_file_netzcore[:-3]+"sif", f, percentage)
    return

def score():
    #g = score_network.create_network_from_weight_and_score_files(edge_file_weights = edge_file_netzcore, edge_file_scores = edge_file_netzcore_relevance)
    #node_to_score = score_network.score_by_shortest_paths(g) 
    #scores = node_to_score.values()
    #scores.sort()
    #print scores[-20:]
    score_network.run_and_assess_performance_of_folds(10, edge_file_netzcore, edge_file_netzcore_relevance_prefix, node_file_netzcore, node_file_netzcore_prefix, result_file_prefix)
    return

def prepare(filter_hubs = False, sample_graphs = False):
    """
	Creates necessary files for scoring
    """

    # Create PPI network using BIANA
    if not (os.path.exists(node_file_prefix+".tsv") and os.path.exists(network_file_prefix+".sif")):
	print "Creating BIANA Human interactome"
	prepare_data.create_human_interactome_using_biana(node_file_prefix = node_file_prefix, network_files_prefix = network_file_prefix, network_type="functional", load_from_saved_session=True) #"experimental") 

    # Filter by detection method (non-tap interactions)
    if not os.path.exists(network_file_filtered_by_method): 
	print "Creating method filtered files"
	prepare_data.create_method_filtered_network_files(network_file_method_attribute = network_file_prefix + "_method_id.eda", network_file_method_filtered_prefix = network_file_prefix)

	# Filter network by degree and get largest connected component
	if filter_hubs:
	    print "Filtering by degree"
	    #prepare_data.create_degree_filtered_network_file(network_file_filtered_by_method, network_file_filtered_by_degree, 90)
	    prepare_data.create_degree_filtered_network_file(network_file_filtered_by_method, network_file_filtered_by_degree, 175)
	    network_file_filtered  =  network_file_filtered_by_degree
	else: 
	    print "No degree filtering, using method filtered file!"
	    network_file_filtered = network_file_filtered_by_method

    # Create node scores file and check how many of the seed genes we cover in the network
    if not os.path.exists(node_scores_file): 
	print "Creating node scores file"
	prepare_data.create_node_scores_file(network_file = network_file_filtered, network_file_identifier_type = "user entity id", node_biana_file = node_file_prefix + ".tsv", association_scores_file = association_scores_file, association_scores_file_identifier_type = "genesymbol", node_scores_file = node_scores_file, ignored_seed_nodes = None, default_non_seed_score = DEFAULT_NON_SEED_SCORE)

    return

    if not os.path.exists(sampling_dir): 
	os.mkdir(sampling_dir)
	
    if sample_graphs:
	print "Creating sampled networks"
	prepare_data.sample_network_preserving_topology(edge_file_netzcore_relevance[:-4] + ".sif", 100, sampled_file_prefix + ".sif.")
    #prepare_data.sample_network_preserving_topology("../data/toy_data/test_interactions_small.sif", 4, "../data/sampled_graphs/sampled_graph.sif.")

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

    #create_arff_file(node_score_file, network_file_filtered, arff_file)

    #prepare_data.generate_cross_validation_test_test_arff_files(10, arff_file, arff_file+".id_score_filtered.arff", arff_file_prefix)
    return


def test_timing_score_one_fold():
    from time import clock
    t1 = clock()
    print score_network.test_run(edge_file_netzcore)
    t2 = clock()
    print "Time: ", t2-t1
    return


if __name__ == "__main__":
    main()

