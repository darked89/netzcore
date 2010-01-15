#!/usr/bin/env python

#########################################################################
# Main workflow for scoring 
#
# eg 25/06/2009
#########################################################################

import prepare_data
import analyze_results
import os, os.path
from string import Template

# obsolete
import score_network # obsolete

#only_print_command = True
only_print_command = False

# Scoring related parameters 
MODE = "all" # prepare, score, analyze

#PPI = "biana" # biana output network as it is (do not forget to revert _degree_filtered.sif.original to _degree_filtered.sif, this one is _degree_filtere_disconnected_only.sif)
#PPI = "biana_no_reliability" # only largest connected component (lcc)
#PPI = "biana_reliability" # reliability filtered lcc
PPI = "biana_no_tap_no_reliability" # tap filtered lcc
#PPI = "biana_no_tap_no_reliability_1e-5" # tap filtered lcc with non seed scores of 1e-5
#PPI = "biana_no_tap_reliability" # tap & reliability filtered lcc
#PPI = "biana_no_tap_relevance" # tap filtered & string edge score assigned lcc
#PPI = "biana_no_tap_exp_db_relevance" tap filtered & string exp & db edge score assigned lcc
#PPI = "biana_no_tap_corelevance" # tap filtered & string co-exp score assigned lcc
#PPI = "biana_no_tap_reliability_relevance" # tap & reliability filtered & string edge score assigned lcc
#PPI = "goh" 
#PPI = "rhodes"
#PPI = "ori_coexpression_1e-2"
#PPI = "ori_network" # these are all in 1e-5
#PPI = "ori_coexpression" 
#PPI = "ori_coexpression_colocalization"
#PPI = "ori_colocalization"
#PPI = "ori_coexpression_colocalization_1e-2"

ASSOCIATION = "aneurysm"
#ASSOCIATION = "breast_cancer"
#ASSOCIATION = "apoptosis"

#SCORING = "ns" #"netscore"
SCORING = "nz" #"netzcore"
#SCORING = "nh" #"netzscore"
#SCORING = "n1" #"netz1score"
#SCORING = "nd" #"netshort"
#SCORING = "nw" #"netween"
#SCORING = "nr" #"netrank"
#SCORING = "nx" #"netrandom"
#SCORING = "nb" #"netZscore" (cortesy of baldo)
#SCORING = "ff" #"FunctionalFlow" 

N_REPETITION = 1
N_ITERATION = 5

DEFAULT_NON_SEED_SCORE = 0.01 # Default score for non-seed nodes
ALLOWED_MAX_DEGREE = 100000 #175 #90 # Max degree allowed in the graph filtering
N_SAMPLE_GRAPH = 100 # Number of random graphs to be generated
N_X_VAL = 5 # Number of cross validation folds
N_RANDOM_NEGATIVE_FOLDS = None #10 # Number of non-seed scores to be averaged for negative score calculation
REPLICABLE = True #False # Assign a predefined seed in randomization for N_RANDOM_NEGATIVE_FOLD generation

# Data directory of the project
data_dir = ".." + os.sep + "data" + os.sep
data_dir = os.path.abspath(data_dir) + os.sep

# BIANA node & network files
biana_node_file_prefix = data_dir + "human_interactome_biana" + os.sep + "human_nodes"
biana_network_file_prefix = data_dir + "human_interactome_biana" + os.sep + "human_network"

# PPIs from existing studies
goh_network_file = data_dir + "goh07_human_ppi" + os.sep + "ppi.sif"
rhodes_network_file = data_dir + "rhodes05_human_probabilistic_ppi" + os.sep + "ppi.sif"
goh_network_file_filtered_by_degree = goh_network_file[:-4] + "_degree_filtered.sif"
rhodes_network_file_filtered_by_degree = rhodes_network_file[:-4] + "_degree_filtered.sif"

# Gene info file 
gene_info_file = data_dir + "gene_info" + os.sep + "genes.tsv"

# Disease association files (Association data to be used)
association_scores_validation_file = None
if ASSOCIATION == "aneurysm":
    association_dir = data_dir + "aneurist" + os.sep
    #aneurysm_scores_file = association_dir + "aneurysm_associated_genes.txt"
    association_scores_file = association_dir + "aneurysm_associated_genes_all_equal.txt"
    association_scores_file_identifier_type = "genesymbol"
    association_scores_validation_file = association_dir + "aneurysm_new_9.txt" 
elif ASSOCIATION == "breast_cancer":
    association_dir = data_dir + "osiris" + os.sep
    association_scores_file = association_dir + "breast_cancer_genes_all_equal.txt"
    association_scores_file_identifier_type = "genesymbol"
#elif ASSOCIATION == "apoptosis":
#    association_dir = data_dir + "apoptosis" + os.sep
#    association_scores_file = apoptosis_scores_all_equal_file
#    association_scores_file_identifier_type = "uniprotaccession"
else:
    raise ValueError("Unrecognized association!")


# Network specific
interaction_relevance_file = None
interaction_relevance_file2 = None
# BIANA ppi
if PPI.startswith("biana"):
    node_description_file = biana_node_file_prefix + ".tsv"
    network_file_identifier_type = "user entity id"
    if PPI == "biana":
	biana_network_file_filtered_by_method = biana_network_file_prefix + ".sif"
	biana_network_file_filtered_by_reliability = biana_network_file_filtered_by_method[:-4] + "_reliability_filtered.sif"
	network_file = biana_network_file_filtered_by_method 
    elif PPI == "biana_no_reliability":
	biana_network_file_filtered_by_method = biana_network_file_prefix + ".sif"
	biana_network_file_filtered_by_reliability = biana_network_file_filtered_by_method[:-4] + "_reliability_filtered.sif"
	network_file = biana_network_file_filtered_by_method  
    elif PPI == "biana_reliability":
	biana_network_file_filtered_by_method = biana_network_file_prefix + ".sif"
	biana_network_file_filtered_by_reliability = biana_network_file_filtered_by_method[:-4] + "_reliability_filtered.sif"
	network_file = biana_network_file_filtered_by_reliability 
    elif PPI.startswith("biana_no_tap"):
	biana_network_file_filtered_by_method = biana_network_file_prefix + "_no_tap.sif"
	biana_network_file_filtered_by_reliability = biana_network_file_filtered_by_method[:-4] + "_reliability_filtered.sif"
	if PPI == "biana_no_tap_no_reliability":
	    network_file = biana_network_file_filtered_by_method # Using only non-tap interactions
	elif PPI == "biana_no_tap_no_reliability_1e-5":
	    network_file = biana_network_file_filtered_by_method # Using only non-tap interactions
	    DEFAULT_NON_SEED_SCORE = 0.00001 
	elif PPI == "biana_no_tap_reliability":
	    network_file = biana_network_file_filtered_by_reliability # Using reliability filtered non-tap interactions
	elif PPI == "biana_no_tap_relevance":
	    network_file = biana_network_file_filtered_by_method 
	    interaction_relevance_file = biana_network_file_prefix + "_stringscore.eda"
	elif PPI == "biana_no_tap_exp_db_relevance":
	    network_file = biana_network_file_filtered_by_method 
	    interaction_relevance_file = biana_network_file_prefix + "_stringscore_experimental.eda"
	    interaction_relevance_file2 = biana_network_file_prefix + "_stringscore_db.eda"
	elif PPI == "biana_no_tap_corelevance":
	    network_file = biana_network_file_filtered_by_method 
	    interaction_relevance_file = biana_network_file_prefix + "_stringscore_coexpression.eda"
	elif PPI == "biana_no_tap_reliability_relevance":
	    network_file = biana_network_file_filtered_by_reliability 
	    interaction_relevance_file = biana_network_file_prefix + "_stringscore.eda"
	else:
	    raise ValueError("Unrecognized ppi!")
    else:
	raise ValueError("Unrecognized ppi!")
    network_file_filtered = network_file[:-4] + "_degree_filtered.sif" # Using only the largest strongly connected component
# Goh ppi
elif PPI == "goh":
    node_description_file = gene_info_file 
    network_file_identifier_type = "geneid"
    network_file = goh_network_file
    network_file_filtered = goh_network_file_filtered_by_degree 
# Rhodes ppi
elif PPI == "rhodes":
    node_description_file = gene_info_file 
    network_file_identifier_type = "geneid"
    network_file = rhodes_network_file
    network_file_filtered = rhodes_network_file_filtered_by_degree 
elif PPI.startswith("ori"):
    network_base_dir = data_dir + "human_interactome_ori" + os.sep
    node_description_file = biana_node_file_prefix + ".tsv"
    network_file_identifier_type = "user entity id"
    if PPI == "ori_network":
	network_dir = network_base_dir + "network" + os.sep
	node_file = network_dir + "aneurist.noa"
	network_file = network_dir + "aneurist.sif"
	network_file_filtered = network_file
	interaction_relevance_file = network_dir + "aneurist.eda"
	DEFAULT_NON_SEED_SCORE = 0.00001 
    elif PPI == "ori_coexpression_1e-2":
	network_dir = network_base_dir + "coexpression" + os.sep
	node_file = network_dir + "aneurist.noa"
	network_file = network_dir + "aneurist.sif"
	network_file_filtered = network_file
	interaction_relevance_file = network_dir + "aneurist.eda"
    elif PPI == "ori_coexpression":
	network_dir = network_base_dir + "coexpression" + os.sep
	node_file = network_dir + "aneurist.noa"
	network_file = network_dir + "aneurist.sif"
	network_file_filtered = network_file
	interaction_relevance_file = network_dir + "aneurist.eda"
	DEFAULT_NON_SEED_SCORE = 0.00001 
    elif PPI == "ori_coexpression_colocalization_1e-2":
	network_dir = network_base_dir + "coexpression_colocalization" + os.sep
	node_file = network_dir + "aneurist.noa"
	network_file = network_dir + "aneurist.sif"
	network_file_filtered = network_file
	interaction_relevance_file = network_dir + "aneurist.eda"
    elif PPI == "ori_coexpression_colocalization":
	network_dir = network_base_dir + "coexpression_colocalization" + os.sep
	node_file = network_dir + "aneurist.noa"
	network_file = network_dir + "aneurist.sif"
	network_file_filtered = network_file
	interaction_relevance_file = network_dir + "aneurist.eda"
	DEFAULT_NON_SEED_SCORE = 0.00001 
    elif PPI == "ori_colocalization":
	network_dir = network_base_dir + "colocalization" + os.sep
	node_file = network_dir + "aneurist.noa"
	network_file = network_dir + "aneurist.sif"
	network_file_filtered = network_file
	interaction_relevance_file = network_dir + "aneurist.eda"
	DEFAULT_NON_SEED_SCORE = 0.00001 
else:
    raise ValueError("Unrecognized ppi!")


# Human readable title for the run
title = "%s - %s - %s" % (PPI, ASSOCIATION, SCORING)

# Project directory structure
input_base_dir = data_dir + "input" + os.sep
input_base_dir_network = input_base_dir + PPI + os.sep
input_base_dir_association = input_base_dir_network + ASSOCIATION + os.sep
input_dir = input_base_dir_association
sampling_dir = input_base_dir_network + "sampled_graphs" + os.sep
output_base_dir = data_dir + "output" + os.sep
output_base_dir_network = output_base_dir + PPI + os.sep
output_base_dir_association = output_base_dir_network + ASSOCIATION + os.sep
output_dir = output_base_dir_association + SCORING + os.sep

# Create directory hirerarchy
if not os.path.exists(input_base_dir): 
    os.mkdir(input_base_dir)
if not os.path.exists(input_base_dir_network): 
    os.mkdir(input_base_dir_network)
if not os.path.exists(input_base_dir_association): 
    os.mkdir(input_base_dir_association)
if not os.path.exists(input_dir): 
    os.mkdir(input_dir)
if not os.path.exists(sampling_dir): 
    os.mkdir(sampling_dir)
if not os.path.exists(output_base_dir): 
    os.mkdir(output_base_dir)
if not os.path.exists(output_base_dir_network): 
    os.mkdir(output_base_dir_network)
if not os.path.exists(output_base_dir_association): 
    os.mkdir(output_base_dir_association)
if not os.path.exists(output_dir): 
    os.mkdir(output_dir)

if SCORING == "ns" or SCORING == "nh" or SCORING == "n1":
    title += " - r%d i%d" % (N_REPETITION, N_ITERATION)
    output_dir = output_dir + "r%di%d" % (N_REPETITION, N_ITERATION) + os.sep
    if not os.path.exists(output_dir): 
	os.mkdir(output_dir)
elif SCORING == "nz" or SCORING == "ff" or SCORING == "nb":
    title += " - i%d" % N_ITERATION
    output_dir = output_dir + "i%d" % N_ITERATION + os.sep
    if not os.path.exists(output_dir): 
	os.mkdir(output_dir)

# Input/Output score file
seed_scores_file = input_dir + "seed_scores.sif"
node_scores_file = input_dir + "node_scores.sif"
#edge_scores_file = input_dir + "edge_scores.sif"
edge_scores_file = input_base_dir_network + "edge_scores.sif"
edge_scores_as_node_scores_file = input_dir + "edge_scores_as_node_scores.sif"
reliability_filtered_edge_scores_file = input_base_dir_network + "reliability_filtered_edge_scores.sif"
output_scores_file = output_dir + "node_scores.sif" # "_r%dn%d.sif" % (N_REPETITION, N_ITERATION)
score_log_file = output_dir + "log.txt" # "_r%dn%d.txt" % (N_REPETITION, N_ITERATION)
sampled_file_prefix = sampling_dir + "sampled_graph"

# Log (README) files 
log_file = output_dir + "README"
input_log_file = input_dir + "README"

# Analysis files
predictions_file = output_dir + "predictions.txt"
labels_file = output_dir + "labels.txt"
r_script_file = output_dir + "results.r"
tex_script_file = output_dir + "results.tex"

score_xval_commands = { "ns": Template("scoreNetwork/scoreN -s s -n %s.$fold -e %s -o %s.$fold -r %d -i %d &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, score_log_file)),
			"nz": Template("scoreNetwork/scoreN -s z -n %s.$fold -e %s -o %s.$fold -i %d -x %d -d %s &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file)),
			"nh": Template("scoreNetwork/scoreN -s h -n %s.$fold -e %s -o %s.$fold -r %d -i %d -x %d -d %s &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file)),
			"n1": Template("scoreNetwork/scoreN -s 1 -n %s.$fold -e %s -o %s.$fold -r %d -i %d -x %d -d %s &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file)),
			"nd": Template("scoreNetwork/scoreN -s d -n %s.$fold -e %s.$fold -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_as_node_scores_file, output_scores_file, score_log_file)),
			"nw": Template("scoreNetwork/scoreN -s w -n %s.$fold -e %s.$fold -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_as_node_scores_file, output_scores_file, score_log_file)),
			"nr": Template("scoreNetwork/scoreN -s r -n %s.$fold -e %s -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file)),
			"nx": Template("scoreNetwork/scoreN -s x -n %s.$fold -e %s -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file)),
			"nb": Template("./netscore -c %s.$fold -i %s -o %s.$fold -t 0 -z 0 -nr 100 -r 1 -zp 0 -n %d -nd 2 -mx 1 -ms 3 -mn 0 -dn 2 -de 2 -mxe 0 -mne 0.00000001 -mnd 0.0000001 -mnde 0.0000001 -mnst 20 -mnste 20 -dxi 1 -dxn 0 -dxe 0 -e 0.0000001 &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, score_log_file)),
			"ff": Template("./fFlow %s.$fold %s %s.$fold %d &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, score_log_file)),
		      }

score_commands = { "ns": "scoreNetwork/scoreN -s s -n %s -e %s -o %s -r %d -i %d &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, score_log_file),
		   "nz": "scoreNetwork/scoreN -s z -n %s -e %s -o %s -i %d -x %d -d %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file), 
		   "nh": "scoreNetwork/scoreN -s h -n %s -e %s -o %s -r %d -i %d -x %d -d %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file), 
		   "n1": "scoreNetwork/scoreN -s 1 -n %s -e %s -o %s -r %d -i %d -x %d -d %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file), 
		   "nd": "scoreNetwork/scoreN -s d -n %s -e %s -o %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file),
		   "nw": "scoreNetwork/scoreN -s w -n %s -e %s -o %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file),
		   "nr": "scoreNetwork/scoreN -s r -n %s -e %s -o %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file),
		   "nx": "scoreNetwork/scoreN -s x -n %s -e %s -o %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file),
		   "nb": "./netscore -c %s -i %s -o %s -t 0 -z 0 -nr 100 -r 1 -zp 0 -n %d -nd 2 -mx 1 -ms 3 -mn 0 -dn 2 -de 2 -mxe 0 -mne 0.00000001 -mnd 0.0000001 -mnde 0.0000001 -mnst 20 -mnste 20 -dxi 1 -dxn 0 -dxe 0 -e 0.0000001 &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, score_log_file),
		   "ff": "./fFlow %s %s %s %d &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, score_log_file),
		 }


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
    # Create PPI network if necessary
    if PPI.startswith("biana"):
	create_biana_network_files(biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, biana_network_file_filtered_by_reliability)

    # Filter network by degree and get largest connected component
    if not os.path.exists(network_file_filtered): 
	print "Filtering by degree", network_file, "->", network_file_filtered
	prepare_data.analyze_network(network_file)
	prepare_data.create_degree_filtered_network_file(network_file, network_file_filtered, ALLOWED_MAX_DEGREE)
    prepare_data.analyze_network(network_file_filtered, out_file = input_log_file)

    # Get node to association mapping
    if PPI.startswith("ori"):
	os.system("awk '{ if($2 == 1) print $1, $2}' %s > %s" % (node_file, seed_scores_file))

    seed_to_score = None
    if not os.path.exists(seed_scores_file): 
	seed_to_score = prepare_data.get_node_association_score_mapping(network_file = network_file_filtered, network_file_identifier_type = network_file_identifier_type, node_description_file = node_description_file, association_scores_file = association_scores_file, association_scores_file_identifier_type = association_scores_file_identifier_type, log_file = input_log_file)
	# Create initial data analysis file
	#if not os.path.exists(input_dir + "analyze_network.r"):
	prepare_data.create_R_analyze_network_script(network_file_filtered, seeds=seed_to_score.keys(), out_path=input_dir, title = PPI + "_" + ASSOCIATION)

    # Prepare scoring files
    prepare_scoring_files(network_file_filtered, seed_to_score, node_scores_file, edge_scores_file, sampled_file_prefix)
    return


def score():
    """
	Runs or prints commands to run scoring method on the input files
    """
    score_original()
    score_xval()
    return


def analyze():
    """
	Does cross validation and percentage analysis on the output files
    """
    analyze_original()
    analyze_xval_percentage()
    analyze_xval() 
    return


def generate_score_xval_command(k):
    #return score_xval_command % (node_scores_file, k, edge_scores_file, output_scores_file, k, N_REPETITION, N_ITERATION, score_log_file, k)
    return score_xval_commands[SCORING].substitute(fold = '%d' % k)


def score_xval():
    if not os.path.exists(output_scores_file + ".1"):
	f = open(log_file, "a")
	for k in range(1, N_X_VAL+1):
	    if only_print_command:
		print generate_score_xval_command(k)
	    else:
		f.write("%s\n" % generate_score_xval_command(k))
		os.system( generate_score_xval_command(k) )
	f.close()
    return


def score_original():
    if not os.path.exists(output_scores_file):
	f = open(log_file, "a")
	#print score_commands[SCORING]
	if only_print_command:
	    print score_commands[SCORING]
	else:
	    f.write("%s\n" % score_commands[SCORING])
	    os.system(score_commands[SCORING])
	f.close()
    return


def analyze_xval():
    if os.path.exists(r_script_file):
	return
    list_node_scores_and_labels = []
    for k in range(1, N_X_VAL+1):
	node_validation_data = analyze_results.get_validation_node_scores_and_labels(file_result = output_scores_file+".%d"%k, file_seed_test_scores = node_scores_file+".%d.test"%k, file_node_scores = node_scores_file, n_random_negative_folds = N_RANDOM_NEGATIVE_FOLDS, default_score = DEFAULT_NON_SEED_SCORE, replicable = REPLICABLE)
	list_node_scores_and_labels.append(node_validation_data)
    analyze_results.create_ROCR_files(list_node_scores_and_labels, predictions_file, labels_file)
    analyze_results.create_R_script(r_script_file, output_dir, title) # os.path.basename(output_dir))
    os.system("R CMD BATCH %s" % (r_script_file))
    os.system("convert %sperformance.eps %sperformance.jpg" % (output_dir, output_dir))
    analyze_results.create_tex_script(tex_script_file, output_dir, title)

    #for tScore in [ i*0.01 for i in xrange(0, 100, 5)]:
    f = open(log_file, "a")
    for tScore in [ 0.5 ]:
	#print "---- t:", tScore
	f.write("---- t: %s\n" % tScore)
	nTP_sum, nFP_sum, nFN_sum, nTN_sum = 0.0, 0.0, 0.0, 0.0
	for k in range(1, N_X_VAL+1):
	    ##print output_scores_file+".ns.%d"%k, node_scores_file+".%d.test"%k 
	    nTP, nFP, nFN, nTN = analyze_results.calculate_performance_metric_counts(file_result = output_scores_file+".%d" % k, file_seed_test_scores = node_scores_file+".%d.test"%k, file_node_scores = node_scores_file, score_threshold = tScore, n_random_negative_folds = N_RANDOM_NEGATIVE_FOLDS, default_score = DEFAULT_NON_SEED_SCORE, replicable = REPLICABLE)
	    (acc, sens, spec, ppv) = analyze_results.calculatePerformance(nTP, nFP, nFN, nTN)
	    ##print "A:", acc, "S:", sens, "P:", ppv
	    nTP_sum += nTP
	    nFP_sum += nFP
	    nFN_sum += nFN
	    nTN_sum += nTN
	# Calculate Sensitivity (TP/T[TP+FN]) (aka TPR) & Calculate PPV (TP/P[TP+FP]) (use randomly selected non-seed nodes as negatives [average w.r.t. n different random selection results])
	nTP = nTP_sum / N_X_VAL
	nFP = nFP_sum / N_X_VAL
	nFN = nFN_sum / N_X_VAL
	nTN = nTN_sum / N_X_VAL
	(acc, sens, spec, ppv) = analyze_results.calculatePerformance(nTP, nFP, nFN, nTN)
	#print "Avg. A:", acc, "S:", sens, "P:", ppv
	f.write("Avg. A: %s S: %s P: %s\n" % (acc, sens, ppv))
	# Calculate 1-Specificity (1-TN/N[FP+TN] = FP/N) (aka FPR) 
	fpr = 1-spec
    f.close()
    return


def analyze_xval_percentage():
    f = open(log_file, "a")
    for percentage in (10, 25, 50):
	#print "---- %s%%:" % percentage
	f.write("---- %s%%:\n" % percentage)
	n_seed_sum = 0.0
	for k in range(1, N_X_VAL+1):
	    ##print output_scores_file+".%s.%d" % (SCORING, k), node_scores_file+".%d.test"%k 
	    n_seed, n, i = analyze_results.calculate_seed_coverage_at_given_percentage(output_scores_file + ".%d" % k, node_scores_file+".%d.test" % k, percentage, DEFAULT_NON_SEED_SCORE)
	    ##print n_seed, n, i
	    if i > n+1:
		#print "Warning: Real coverage percentage is disputed due to equal scores!", i, n
		f.write("Warning: Real coverage percentage is disputed due to equal scores! %i %i\n" % (i, n))
	    n_seed_sum += n_seed
	n_seed = n_seed_sum / N_X_VAL
	#print "Avg:", n_seed, "over", n
	f.write("Avg: %i over %i\n" % (n_seed, n))
    f.close()
    return


def analyze_original():
    if not os.path.exists(output_scores_file):
	raise Exception("Output score file does not exist!")
    f = open(log_file, "a")
    for percentage in (10, 25, 50):
	#print "---- %s:" % percentage
	f.write("---- %s:\n" % percentage)
	f.write("%s\n" % output_scores_file)
	f.write("%s\n" % str(analyze_results.calculate_seed_coverage_at_given_percentage(output_scores_file, node_scores_file, percentage, DEFAULT_NON_SEED_SCORE)))
	if not os.path.exists(output_scores_file+"."+association_scores_file_identifier_type):
	    if PPI.startswith("biana"):
		prepare_data.convert_file_using_new_id_mapping(output_scores_file, node_description_file, network_file_identifier_type, "geneid", output_scores_file+".geneid")
		prepare_data.convert_file_using_new_id_mapping(output_scores_file+".geneid", gene_info_file, "geneid", association_scores_file_identifier_type, output_scores_file+"."+association_scores_file_identifier_type)
	    elif PPI.startswith("ori"):
		pass
	    else:
		prepare_data.convert_file_using_new_id_mapping(output_scores_file, node_description_file, network_file_identifier_type, association_scores_file_identifier_type, output_scores_file+"."+association_scores_file_identifier_type)
	if association_scores_validation_file is not None:
	    f.write("Validation seed coverage:")
	    if PPI.startswith("ori"):
		pass
	    else:
		f.write("%s\n" % str(analyze_results.calculate_seed_coverage_at_given_percentage(output_scores_file+"."+association_scores_file_identifier_type, association_scores_validation_file, percentage, DEFAULT_NON_SEED_SCORE)))
    f.close()
    return


def create_biana_network_files(biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, network_file_filtered_by_reliability):
    # Create PPI network using BIANA
    if not (os.path.exists(biana_node_file_prefix+".tsv") and os.path.exists(biana_network_file_prefix+".sif")):
	print "Creating BIANA Human interactome", biana_node_file_prefix+".tsv", biana_network_file_prefix+".sif"
	prepare_data.create_human_interactome_using_biana(node_file_prefix = biana_node_file_prefix, network_files_prefix = biana_network_file_prefix, network_type="functional", load_from_saved_session=False) #"experimental")  
    # Filter by detection method (non-tap interactions)
    if not os.path.exists(biana_network_file_filtered_by_method): 
	print "Creating method filtered files", biana_network_file_prefix + ".sif", "->", biana_network_file_filtered_by_method
	prepare_data.create_method_filtered_network_files(network_file_prefix = biana_network_file_prefix)
    # Filter by reliability based on number of detection methods & sources & pubmeds 
    if not os.path.exists(network_file_filtered_by_reliability): 
	print "Creating reliability filtered files", biana_network_file_filtered_by_method, "->", network_file_filtered_by_reliability
	prepare_data.create_edge_reliability_filtered_network_file(network_file = biana_network_file_filtered_by_method, network_file_prefix = biana_network_file_prefix, out_file = network_file_filtered_by_reliability)
    return


def prepare_scoring_files(network_file_filtered, seed_to_score, node_scores_file, edge_scores_file, sampled_file_prefix):
    # Create node scores file and check how many of the seed genes we cover in the network
    # Create seed scores file from association score to seed node mapping
    if not os.path.exists(seed_scores_file): 
	print "Creating seed scores file", seed_scores_file
	seeds = seed_to_score.keys()
	prepare_data.create_node_scores_file(nodes = seeds, node_to_score = seed_to_score, node_scores_file = seed_scores_file, ignored_nodes = None, default_score = DEFAULT_NON_SEED_SCORE)
    # Create node scores (original + xval) files using previously created seed score file (all non-seeds will have DEFAULT_NON_SEED_SCORE)
    if not os.path.exists(node_scores_file): 
	print "Creating node score files", node_scores_file
	nodes = prepare_data.get_nodes_in_network(network_file = network_file_filtered)
	seed_to_score = prepare_data.get_node_to_score_from_node_scores_file(seed_scores_file)
	prepare_data.create_node_scores_file(nodes = nodes, node_to_score = seed_to_score, node_scores_file = node_scores_file, ignored_nodes = None, default_score = DEFAULT_NON_SEED_SCORE)
	prepare_data.generate_cross_validation_node_score_files(nodes = nodes, seed_to_score = seed_to_score, node_scores_file = node_scores_file, xval = N_X_VAL, default_score = DEFAULT_NON_SEED_SCORE)
    # Create edge scores (original + xval) files as well as node scores as edge scores files
    if not os.path.exists(edge_scores_file): 
	print "Creating edge score files", edge_scores_file
	edges = prepare_data.get_edges_in_network(network_file = network_file_filtered)
	if all([interaction_relevance_file, interaction_relevance_file2]):
	    edge_to_score = prepare_data.get_edge_to_score_from_sif_attribute_file(interaction_relevance_file)
	    edge_to_score = dict([(e, (sum([float(i) for i in v])/len(v))/1000 + 0.001) for e,v in edge_to_score.iteritems()])
	    edge_to_score2 = prepare_data.get_edge_to_score_from_sif_attribute_file(interaction_relevance_file2)
	    edge_to_score2 = dict([(e, (sum([float(i) for i in v])/len(v))/1000) for e,v in edge_to_score2.iteritems()])
	    edge_to_score = dict([(e, v+edge_to_score2[e]) for e,v in edge_to_score.iteritems()])
	elif interaction_relevance_file is not None:
	    edge_to_score = prepare_data.get_edge_to_score_from_sif_attribute_file(interaction_relevance_file)
	    # Normalizing string score (0.001 is added to avoid 0 score edges)
	    edge_to_score = dict([(e, (sum([float(i) for i in v])/len(v))/1000 + 0.001) for e,v in edge_to_score.iteritems()])
	else:
	    edge_to_score = dict([(e, 1) for e in edges])
	prepare_data.create_edge_scores_file(network_file = network_file_filtered, edge_scores_file = edge_scores_file, edge_to_score = edge_to_score, default_score=DEFAULT_NON_SEED_SCORE)
	#prepare_data.old_create_edge_scores_as_node_scores_file(network_file = network_file_filtered, edge_scores_file = edge_scores_as_node_scores_file, ignored_nodes = None, default_score = DEFAULT_NON_SEED_SCORE)
	seed_to_score = prepare_data.get_node_to_score_from_node_scores_file(seed_scores_file)
	prepare_data.create_edge_scores_as_node_scores_file(edges = edges, node_to_score = seed_to_score, edge_scores_file = edge_scores_as_node_scores_file, ignored_nodes = None, default_score = DEFAULT_NON_SEED_SCORE)
	prepare_data.generate_cross_validation_edge_score_as_node_score_files(edges = edges, seed_to_score = seed_to_score, edge_scores_file = edge_scores_as_node_scores_file, xval = N_X_VAL, default_score = DEFAULT_NON_SEED_SCORE)
    # Create random network files
    if not os.path.exists(sampled_file_prefix + ".sif.1"): 
	print "Creating sampled networks"
	prepare_data.sample_network_preserving_topology(edge_scores_file, N_SAMPLE_GRAPH, sampled_file_prefix + ".sif.")
    return


########################### Oldies & Goldies #########################
######################################################################

def old_score_xval_random():
    if not os.path.exists(output_scores_file+"random.1"):
	for k in range(1, N_X_VAL+1):
	    score_network.score_by_random_model(node_scores_file + ".%d" % k, edge_scores_file, output_scores_file + ".random.%d" % k)
    return


def old_prepare():
    # Get scores of each node
    node_to_score, seeds_all = prepare_data.map_scores_to_biana_nodes(node_file, node_file_score)

    # Assign node scores either as node relevance score or as edge relevance score
    prepare_data.assign_node_scores(g = g, node_to_score = node_to_score, out_file = node_file_netzcore, as_edge_relevance_score = False, seeds_ignored=None)
    prepare_data.assign_node_scores(g = g, node_to_score = node_to_score, out_file = edge_file_netzcore_relevance, as_edge_relevance_score = True, seeds_ignored=None)
    
    return
    
    # Create arff file for all network
    #network_utilities.create_arff_file_with_network_metrics(g, node_to_score, seeds, arff_file)

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

