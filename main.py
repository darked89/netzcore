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
MODE = "all" # prepare, score, analyze

PPI = "biana" 
#PPI = "goh" 
#PPI = "rhodes"

N_REPEATITION = 3
N_ITERATION = 3

DEFAULT_NON_SEED_SCORE = 0.01
ALLOWED_MAX_DEGREE = 100000 #175 #90
N_SAMPLE_GRAPH = 100
N_X_VAL = 5

# Project directory structure
data_dir = ".." + os.sep + "data" + os.sep

# BIANA node & network files
biana_node_file_prefix = data_dir + "human_interactome_biana" + os.sep + "human_nodes"
biana_network_file_prefix = data_dir + "human_interactome_biana" + os.sep + "human_network"
biana_network_file_filtered_by_method = biana_network_file_prefix + "_no_tap.sif"
biana_network_file_filtered_by_degree = biana_network_file_filtered_by_method[:-4] + "_degree_filtered.sif"

# Edge attribute files
biana_network_file_method_attribute = biana_network_file_prefix + "_method_id.eda"
biana_network_file_source_attribute = biana_network_file_prefix + "_source.eda"
biana_network_file_pubmed_attribute = biana_network_file_prefix + "_pubmed.eda"

# PPIs from existing studies
goh_network_file = data_dir + "goh07_human_ppi" + os.sep + "ppi.sif"
rhodes_network_file = data_dir + "rhodes05_human_probabilistic_ppi" + os.sep + "ppi.sif"
goh_network_file_filtered_by_degree = goh_network_file[:-4] + "_degree_filtered.sif"
rhodes_network_file_filtered_by_degree = rhodes_network_file[:-4] + "_degree_filtered.sif"

# Gene info file 
gene_info_file = data_dir + "gene_info" + os.sep + "genes.tsv"

# Disease association files
aneurysm_scores_file = data_dir + "aneurist" + os.sep + "aneurysm_associated_genes.txt"
aneurysm_scores_all_equal_file = data_dir + "aneurist" + os.sep + "aneurysm_associated_genes_all_equal.txt"


# Network specific

# BIANA ppi
if PPI == "biana":
    node_description_file = biana_node_file_prefix + ".tsv"
    network_file = biana_network_file_filtered_by_method
    network_file_filtered = biana_network_file_filtered_by_degree 
    input_dir = data_dir + "input_biana" + os.sep
    output_dir = data_dir + "output_biana" + os.sep
    network_file_identifier_type = "user entity id"
# Goh ppi
elif PPI == "goh":
    node_description_file = gene_info_file 
    network_file = goh_network_file
    network_file_filtered = goh_network_file_filtered_by_degree 
    input_dir = data_dir + "input_goh" + os.sep
    output_dir = data_dir + "output_goh" + os.sep
    network_file_identifier_type = "geneid"
# Rhodes ppi
elif PPI == "rhodes":
    node_description_file = gene_info_file 
    network_file = rhodes_network_file
    network_file_filtered = rhodes_network_file_filtered_by_degree 
    input_dir = data_dir + "input_rhodes" + os.sep
    output_dir = data_dir + "output_rhodes" + os.sep
    network_file_identifier_type = "geneid"
else:
    raise ValueError("Unrecognized ppi!")

# Common to all ppi networks
association_scores_file = aneurysm_scores_all_equal_file
association_scores_file_identifier_type = "genesymbol"

node_scores_file = input_dir + "node_scores.sif"
edge_scores_file = input_dir + "edge_scores.sif"
sampling_dir = input_dir + "sampled_graphs" + os.sep
output_scores_file = output_dir + "node_scores_r%dn%d.sif" % (N_REPEATITION, N_ITERATION)
log_file = output_dir + "log_r%dn%d.txt" % (N_REPEATITION, N_ITERATION)
sampled_file_prefix = sampling_dir + "sampled_graph"

def main():
    if MODE == "prepare":
	prepare()
    elif MODE == "score":
	score()
    elif MODE == "analyze":
	analyze()
    elif MODE == "all":
	prepare()
	score()
	analyze()
    else:
	raise ValueError("Unrecognized mode!")
    return


def prepare():
    """
	Creates necessary files for scoring
    """
    if not os.path.exists(input_dir): 
	os.mkdir(input_dir)
    if not os.path.exists(sampling_dir): 
	os.mkdir(sampling_dir)
    if not os.path.exists(output_dir): 
	os.mkdir(output_dir)

    if PPI == "biana":
	create_biana_network()

    prepare_scoring_files()

    return

def score():
    #score_original()
    score_xval()
    return

def analyze():
    #analyze_original()
    analyze_xval_percentage()
    analyze_xval()
    return


def score_xval():
    score_command = "scoreNetwork/scoreN -s s -n %s.%d -e %s -o %s.ns.%d -r %d -i %d &> %s.ns.%d" 
    for k in range(1, N_X_VAL+1):
	print score_command % (node_scores_file, k, edge_scores_file, output_scores_file, k, N_REPEATITION, N_ITERATION, log_file, k)
	os.system( score_command % (node_scores_file, k, edge_scores_file, output_scores_file, k, N_REPEATITION, N_ITERATION, log_file, k) )
    return


def score_original():
    netscore_command = "scoreNetwork/scoreN -s s -n %s -e %s -o %s.ns -r %d -i %d &> %s.ns" % (node_scores_file, edge_scores_file, output_scores_file, N_REPEATITION, N_ITERATION, log_file)
    netzcore_command = "scoreNetwork/scoreN -s z -n %s -e %s -o %s.nz -i %d -x %d -d %s &> %s.nz" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, log_file)
    if not os.path.exists(output_scores_file+".ns"):
	print netscore_command
	os.system(netscore_command)
    if not os.path.exists(output_scores_file+".nz"):
	pass
	#print netzcore_command
	#os.system(netzcore_command)
    return


def analyze_xval():
    for tScore in [ i*0.1 for i in xrange(0, 10,1)]:
	print "t:", tScore
	nTP_sum, nFP_sum, nFN_sum, nTN_sum = 0.0, 0.0, 0.0, 0.0
	for k in range(1, N_X_VAL+1):
	    print output_scores_file+".ns.%d"%k, node_scores_file+".%d.test"%k 
	    nTP, nFP, nFN, nTN = analyze_results.calculate_performance_metric_counts(file_result = output_scores_file+".ns.%d"%k, file_seed_test_scores = node_scores_file+".%d.test"%k, file_seed_scores = node_scores_file, score_threshold = tScore, n_random_negative_folds = 10)
	    (acc, sens, spec, ppv) = analyze_results.calculatePerformance(nTP, nFP, nFN, nTN)
	    print "A:", acc, "S:", sens, "P:", ppv
	    nTP_sum += nTP
	    nFP_sum += nFP
	    nFN_sum += nFN
	    nTN_sum += nTN
	# Calculate Sensitivity (TP/T[TP+FN]) (aka TPR) & Calculate PPV (TP/P[TP+FP]) (use randomly selected non-seed nodes as negatives [average w.r.t. n different random selection results])
	(acc, sens, spec, ppv) = analyze_results.calculatePerformance(nTP, nFP, nFN, nTN)
	print "Avg. A:", acc, "S:", sens, "P:", ppv
	# Calculate 1-Specificity (1-TN/N[FP+TN] = FP/N) (aka FPR) 
	fpr = 1-spec
    return


def analyze_xval_percentage():
    for percentage in (10, 25, 50):
	print "---- %s:" % percentage
	for k in range(1, N_X_VAL+1):
	    print output_scores_file+".ns.%d"%k, node_scores_file+".%d.test"%k 
	    print analyze_results.calculate_seed_coverage_at_given_percentage(output_scores_file+".ns.%d"%k, node_scores_file+".%d.test"%k, percentage, DEFAULT_NON_SEED_SCORE)
    return


def analyze_original():
    #result_files = [ output_scores_file+".ns", output_scores_file+".nz" ]
    result_files = [ output_scores_file+".ns" ]
    for percentage in (10, 25, 50):
	print "---- %s:" % percentage
	for result_file in result_files:
	    print result_file
	    print analyze_results.calculate_seed_coverage_at_given_percentage(result_file, node_scores_file, percentage, DEFAULT_NON_SEED_SCORE)
    return


def create_biana_network():
    # Create PPI network using BIANA
    if not (os.path.exists(biana_node_file_prefix+".tsv") and os.path.exists(biana_network_file_prefix+".sif")):
	print "Creating BIANA Human interactome", biana_node_file_prefix+".tsv", biana_network_file_prefix+".sif"
	prepare_data.create_human_interactome_using_biana(node_file_prefix = biana_node_file_prefix, network_files_prefix = biana_network_file_prefix, network_type="functional", load_from_saved_session=False) #"experimental")  

    # Filter by detection method (non-tap interactions)
    if not os.path.exists(biana_network_file_filtered_by_method): 
	print "Creating method filtered files",  biana_network_file_prefix + ".sif", "->", biana_network_file_filtered_by_method
	prepare_data.create_method_filtered_network_files(network_file_method_attribute = biana_network_file_prefix + "_method_id.eda", network_file_method_filtered_prefix = biana_network_file_prefix)
    return


def prepare_scoring_files_from_network_and_association_files(network_file, network_file_filtered, network_file_identifier_type, node_description_file, association_scores_file, association_scores_file_identifier_type, node_scores_file, edge_scores_file, sampled_file_prefix):
    # Filter network by degree and get largest connected component
    if not os.path.exists(network_file_filtered): 
	print "Filtering by degree", network_file, "->", network_file_filtered
	prepare_data.create_degree_filtered_network_file(network_file, network_file_filtered, ALLOWED_MAX_DEGREE)

    # Create node scores file and check how many of the seed genes we cover in the network
    if not os.path.exists(node_scores_file): 
	print "Creating node scores file(s)", node_scores_file
	#prepare_data.create_node_scores_file(network_file = network_file_filtered, network_file_identifier_type = network_file_identifier_type, node_description_file = node_description_file, association_scores_file = association_scores_file, association_scores_file_identifier_type = association_scores_file_identifier_type, node_scores_file = node_scores_file, ignored_seed_nodes = None, default_non_seed_score = DEFAULT_NON_SEED_SCORE)
	nodes = prepare_data.get_nodes_in_network(network_file =  network_file_filtered)
	node_to_score = prepare_data.get_node_association_score_mapping(network_file = network_file_filtered, network_file_identifier_type = network_file_identifier_type, node_description_file = node_description_file, association_scores_file = association_scores_file, association_scores_file_identifier_type = association_scores_file_identifier_type)
	prepare_data.create_node_scores_file(nodes = nodes, node_to_score = node_to_score, node_scores_file = node_scores_file, ignored_nodes = None, default_score = DEFAULT_NON_SEED_SCORE)
	prepare_data.generate_cross_validation_node_score_files(nodes = nodes, node_to_score = node_to_score, node_scores_file = node_scores_file, xval = 5, default_score = DEFAULT_NON_SEED_SCORE)
    if not os.path.exists(edge_scores_file): 
	print "Creating edge scores file", edge_scores_file
	prepare_data.create_edge_scores_file(network_file = network_file_filtered, edge_scores_file = edge_scores_file)
    if not os.path.exists(sampled_file_prefix + ".sif.1"): 
	print "Creating sampled networks"
	prepare_data.sample_network_preserving_topology(edge_scores_file, N_SAMPLE_GRAPH, sampled_file_prefix + ".sif.")
    return


def prepare_scoring_files():
    prepare_scoring_files_from_network_and_association_files(network_file, network_file_filtered, network_file_identifier_type, node_description_file, association_scores_file, association_scores_file_identifier_type, node_scores_file, edge_scores_file, sampled_file_prefix)
    return


def old_prepare():
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
    prepare_data.assign_edge_reliability_scores(g = g, network_file_method_attribute = biana_network_file_method_attribute, network_file_source_attribute = biana_network_file_source_attribute, network_file_pubmed_attribute = biana_network_file_pubmed_attribute, out_file = edge_file_netzcore)

    #create_arff_file(node_score_file, network_file_filtered, arff_file)

    #prepare_data.generate_cross_validation_test_test_arff_files(10, arff_file, arff_file+".id_score_filtered.arff", arff_file_prefix)
    return


def old_analyze():
    result_files = [ netshort_results, netrank_results, netscore_results_1, netscore_results_3,  netscore_results_no_1, netscore_results_no_3, netzcore_results_1, netzcore_results_3 ]
    for percentage in (10, 25, 50):
	print "---- %s:" % percentage
	for f in result_files:
	    print f
	    print analyze_results.calculate_seed_coverage_at_given_percentage(f, node_file_netzcore[:-3]+"sif", percentage)
    return


def old_score():
    #g = score_network.create_network_from_weight_and_score_files(edge_file_weights = edge_file_netzcore, edge_file_scores = edge_file_netzcore_relevance)
    #node_to_score = score_network.score_by_shortest_paths(g) 
    #scores = node_to_score.values()
    #scores.sort()
    #print scores[-20:]
    score_network.run_and_assess_performance_of_folds(10, edge_file_netzcore, edge_file_netzcore_relevance_prefix, node_file_netzcore, node_file_netzcore_prefix, result_file_prefix)
    return


def test_timing_score_one_fold():
    from time import clock
    t1 = clock()
    print score_network.test_run(edge_file_netzcore)
    t2 = clock()
    print "Time: ", t2-t1
    return

# Node & edge score files
#node_scores_file = input_dir + "human_node_scores.sif"
#edge_scores_file = input_dir + "human_edge_scores.sif"
#goh_node_scores_file = input_dir + "goh_node_scores.sif"
#goh_edge_scores_file = input_dir + "goh_edge_scores.sif"
#rhodes_node_scores_file = input_dir + "rhodes_node_scores.sif"
#rhodes_edge_scores_file = input_dir + "rhodes_edge_scores.sif"

# Sampling related
#sampled_file_prefix = sampling_dir + "sampled_graph"
#goh_sampled_file_prefix = sampling_dir[:-1] + "_goh" + os.sep + "sampled_graph"
#rhodes_sampled_file_prefix = sampling_dir[:-1] + "_rhodes" + os.sep + "sampled_graph"


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


if __name__ == "__main__":
    main()

