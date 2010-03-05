#!/usr/bin/env python

#########################################################################
# Main workflow for scoring 
#
# eg 25/06/2009
#########################################################################

import prepare_data
import analyze_results
import calculate_mean_and_sigma
import os, os.path
from string import Template

#only_print_command = True
only_print_command = False
use_cluster = True
#use_cluster = False

N_LINKER_THRESHOLD = 2
DEFAULT_SEED_SCORE = 1.0 # Default score for seed nodes, used when no score given in assoc file
DEFAULT_NON_SEED_SCORE = 0.01 # Default score for non-seed nodes
ALLOWED_MAX_DEGREE = 100000 #175 #90 # Max degree allowed in the graph filtering
N_SAMPLE_GRAPH = 100 # Number of random graphs to be generated
N_X_VAL = 5 #182 # Number of cross validation folds
N_RANDOM_NEGATIVE_FOLDS = None #0 #None #10 # Number of non-seed scores to be averaged for negative score calculation, 
			    # If 0 all non seeds are included as they are, If None all non seeds are included averaged to scale with the size of test seeds
REPLICABLE = True #False # Assign a predefined seed in randomization for N_RANDOM_NEGATIVE_FOLD generation

# Directory of the project
base_dir = ".."
base_dir = os.path.abspath(base_dir) + os.sep
data_dir = base_dir + "data" + os.sep
data_dir = os.path.abspath(data_dir) + os.sep
src_dir = base_dir + "src" + os.sep

# BIANA node & network files
biana_node_file_prefix = data_dir + "human_interactome_biana" + os.sep + "human_nodes"
biana_network_file_prefix = data_dir + "human_interactome_biana" + os.sep + "human_network"

# PPIs from existing studies
goh_network_file = data_dir + "goh07_human_ppi" + os.sep + "ppi.sif"
rhodes_network_file = data_dir + "rhodes05_human_probabilistic_ppi" + os.sep + "ppi.sif"
goh_network_file_filtered_by_degree = goh_network_file[:-4] + "_degree_filtered.sif"
rhodes_network_file_filtered_by_degree = rhodes_network_file[:-4] + "_degree_filtered.sif"
rhodes_interaction_relevance_file = rhodes_network_file[:-4] + ".eda"

# Gene info file 
gene_info_file = data_dir + "gene_info" + os.sep + "genes.tsv"

goh_phenotypes = ["goh_developmental", "goh_connective_tissue", "goh_ear,nose,throat", "goh_endocrine", "goh_psychiatric", "goh_immunological", "goh_neurological", "goh_respiratory", "goh_multiple", "goh_renal", "goh_skeletal", "goh_bone", "goh_dermatological", "goh_cancer", "goh_ophthamological", "goh_metabolic", "goh_nutritional", "goh_muscular", "goh_hematological", "goh_gastrointestinal", "goh_cardiovascular"] #,"goh_unclassified"] 

omim_phenotypes = ["alzheimer", "breast cancer", "diabetes", "insulin", "anemia", "aneurysm", "myopathy", "neuropathy", "obesity", "parkinson disease", "prostate cancer", "hypertension", "leukemia", "lung cancer", "autism", "asthma", "ataxia", "epilepsy", "schizophrenia", "cardiomyopathy", "cataract", "spastic paraplegia", "lymphoma", "mental retardation", "systemic lupus erythematosus"]
omim_phenotypes = [ "omim_" + "_".join(p.split()) for p in omim_phenotypes ]

def main():
    MODE = "prepare" # prepare, score, analyze
    ignore_experiment_failures = False

    ppis = ["goh", "entrez", "biana_no_tap_no_reliability", "biana_no_reliability", "biana_no_tap_relevance"] #["goh"] #["piana_joan_exp", "piana_joan_all"] #["david"] #["goh", "biana_no_tap_no_reliability", "biana_no_reliability", "biana_no_tap_relevance"]
    phenotypes = omim_phenotypes + goh_phenotypes #["apoptosis_joan"] #["alzheimer_david_CpOGU", "alzheimer_david_CpOIN", "alzheimer_david_RpOGU", "alzheimer_david_RpOIN"] #["aneurysm", "breast_cancer"]
    scoring_parameters = []
    #scoring_parameters += [("ff", 1, 5), ("nz", 1, 5), ("ns", 3, 2)] 
    scoring_parameters += [("nr", 1, 1)]
    #scoring_parameters += [("nd", 1, 1)]
    #scoring_parameters += [("nw",1, 1)]
    #scoring_parameters += [("nx", 1, 1)]
    #scoring_parameters += [("ff", 1, i) for i in xrange(4,7)]
    #scoring_parameters += [("nz", 1, i) for i in xrange(4,7)]
    #scoring_parameters += [("ns", r, i) for r in xrange(1,4) for i in xrange(1,4)]
    ##scoring_parameters += [("ff", 1, i) for i in xrange(1,9)]
    ##scoring_parameters += [("nz", 1, i) for i in xrange(1,9)]
    ##scoring_parameters += [("ns", r, i) for r in xrange(1,9) for i in xrange(1,5)]
    ##scoring_parameters += [("ns", r, i) for r in xrange(4,9) for i in xrange(1,3)]
    ##scoring_parameters += [("nh", r, i) for r in (1,2,3) for i in xrange(1,5)]
    ##scoring_parameters += [("n1", r, i) for r in (1,2,3) for i in xrange(1,5)]

    experiments = []
    for ppi in ppis:
	for phenotype in phenotypes:
	    for parameters in scoring_parameters:
		experiments.append((ppi, phenotype, parameters[0], parameters[1], parameters[2]))

    for experiment in experiments:
	PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION = experiment
	print "Running experiment:", experiment
	if ignore_experiment_failures:
	    try:
		run_experiment(MODE, PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION)
	    except:
		print "!Problem!"
	else:
	    run_experiment(MODE, PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION)
    return

# Scoring related parameters 
#PPI = "biana" # biana output network as it is (do not forget to revert _degree_filtered.sif.original to _degree_filtered.sif, this one is _degree_filtere_disconnected_only.sif)
#PPI = "biana_no_reliability" # only largest connected component (lcc)
#PPI = "biana_reliability" # reliability filtered lcc
#PPI = "biana_no_tap_no_reliability" # tap filtered lcc
#PPI = "biana_no_tap_no_reliability_1e-5" # tap filtered lcc with non seed scores of 1e-5
#PPI = "biana_no_tap_reliability" # tap & reliability filtered lcc
#PPI = "biana_no_tap_relevance" # tap filtered & string edge score assigned lcc # manually assigned +1 to all edge scores to reduce max/min edge score ratio
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
#PPI = "david_CpOGU"
#PPI = "david_CpOIN"
#PPI = "david_RpOGU"
#PPI = "david_RpOIN"

#ASSOCIATION = "aneurysm"
#ASSOCIATION = "breast_cancer"
#ASSOCIATION = "apoptosis"
#ASSOCIATION = "alzheimer"

#SCORING = "ns" #"netscore"
#SCORING = "nz" #"netzcore"
#SCORING = "nh" #"netzscore"
#SCORING = "n1" #"netz1score"
#SCORING = "nd" #"netshort"
#SCORING = "nw" #"netween"
#SCORING = "nl" #"netlink"
#SCORING = "nr" #"netrank"
#SCORING = "nx" #"netrandom"
#SCORING = "nb" #"netZscore" (cortesy of baldo)
#SCORING = "ff" #"FunctionalFlow" 

#N_REPETITION = 2
#N_ITERATION = 2

#def run_experiment(MODE, PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION, biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, biana_network_file_filtered_by_reliability, network_file, network_file_filtered, input_log_file, node_file, seed_scores_file, network_file_identifier_type, node_description_file, association_scores_file, association_scores_file_identifier_type, input_dir, node_scores_file, edge_scores_file, sampled_file_prefix, output_scores_file, log_file):
def run_experiment(MODE, PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION):

    # Create experiment parameters
    global DEFAULT_NON_SEED_SCORE 

    # Disease association files (Association data to be used)
    association_scores_validation_file = None
    if ASSOCIATION == "aneurysm":
	association_dir = data_dir + "aneurysm" + os.sep
	#aneurysm_scores_file = association_dir + "aneurysm_associated_genes.txt"
	association_scores_file = association_dir + "aneurysm_associated_genes_all_equal.txt"
	association_scores_file_identifier_type = "genesymbol"
	association_scores_validation_file = association_dir + "aneurysm_new_9.txt" 
    elif ASSOCIATION == "breast_cancer":
	association_dir = data_dir + "breast_cancer" + os.sep
	#association_scores_file = association_dir + "breast_cancer_gene_names_all_equal.txt"
	#association_scores_file_identifier_type = "genesymbol"
	association_scores_file = association_dir + "breast_cancer_gene_ids_all_equal.txt"
	association_scores_file_identifier_type = "geneid"
    elif ASSOCIATION.startswith("alzheimer_david"):
	association_dir = data_dir + "alzheimer_david" + os.sep
	association_scores_file = association_dir + "alzheimer_" + ASSOCIATION[-5:] + "_seed.list"
	association_scores_file_identifier_type = "uniprotentry"
    elif ASSOCIATION == "apoptosis_joan":
        association_dir = data_dir + "apoptosis_joan" + os.sep
        association_scores_file = association_dir + "Apoptosis_uniProt_nodes.list"
        association_scores_file_identifier_type = "uniprotentry"
    elif ASSOCIATION == "alzheimer":
	association_dir = data_dir + "alzheimer" + os.sep
	association_scores_file = None
	association_scores_file_identifier_type = "genesymbol"
    elif ASSOCIATION.startswith("goh_"):
	association_dir = data_dir + "goh_disease_data" + os.sep
	association_scores_file = association_dir + ASSOCIATION + ".txt"
	association_scores_file_identifier_type = "genesymbol"
    elif ASSOCIATION.startswith("omim_"):
	association_dir = data_dir + "omim" + os.sep
	association_scores_file = association_dir + ASSOCIATION + ".txt"
	association_scores_file_identifier_type = "genesymbol"
    else:
	raise ValueError("Unrecognized association!")


    # Network specific
    interaction_relevance_file = None
    interaction_relevance_file2 = None
    biana_network_file_filtered_by_method = None
    biana_network_file_filtered_by_reliability = None
    node_file = None
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
    # Entrez ppi
    elif PPI == "entrez":
	node_description_file = gene_info_file 
	network_file_identifier_type = "geneid"
	network_file = data_dir + "entrez10_human_ppi" + os.sep + "ppi.sif"
	network_file_filtered = network_file[:-4] + "_degree_filtered.sif"
    # Rhodes ppi
    elif PPI == "rhodes":
	node_description_file = gene_info_file 
	network_file_identifier_type = "geneid"
	network_file = rhodes_network_file
	network_file_filtered = rhodes_network_file_filtered_by_degree 
	#interaction_relevance_file = rhodes_interaction_relevance_file # need to rescale / cluster scores because max/min >= 10000
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
    elif PPI == "david": #elif PPI.startswith("david"):
	network_dir = data_dir + "human_interactome_david" + os.sep
	network_file_identifier_type = "user entity id"
	network_file = network_dir + "human_network.sif"
	node_description_file = network_dir + "human_nodes.tsv"
	#if PPI == "david_CpOGU":
	    #association_scores_file = association_dir + "alzheimer_CpOGU_seed.list"
	    #node_file = network_dir + "alzheimer_CpOGU_seed.list"
	    #network_file = network_dir + "alzheimer_CpOGU_network.sif"
	    #node_description_file = network_dir + "alzheimer_CpOGU_network_all.tab"
	    #network_file_filtered = network_file
	#elif PPI == "david_CpOIN":
	    #association_scores_file = association_dir + "alzheimer_CpOIN_seed.list"
	    #node_file = network_dir + "alzheimer_CpOIN_seed.list"
	    #network_file = network_dir + "alzheimer_CpOIN_network.sif"
	    #node_description_file = network_dir + "alzheimer_CpOIN_network_all.tab"
	    #network_file_filtered = network_file
	#elif PPI == "david_RpOGU":
	    #association_scores_file = association_dir + "alzheimer_RpOGU_seed.list"
	    #node_file = network_dir + "alzheimer_RpOGU_seed.list"
	    #network_file = network_dir + "alzheimer_RpOGU_network.sif"
	    #node_description_file = network_dir + "alzheimer_RpOGU_network_all.tab"
	    #network_file_filtered = network_file
	#elif PPI == "david_RpOIN":
	    #association_scores_file = association_dir + "alzheimer_RpOIN_seed.list"
	    #node_file = network_dir + "alzheimer_RpOIN_seed.list"
	    #network_file = network_dir + "alzheimer_RpOIN_network.sif"
	    #node_description_file = network_dir + "alzheimer_RpOIN_network_all.tab"
	    #network_file_filtered = network_file
	network_file_filtered = network_file[:-4] + "_degree_filtered.sif" # Using only the largest strongly connected component
    elif PPI.startswith("piana_joan"):
	network_dir = data_dir + "apoptosis_joan" + os.sep
	node_description_file = None #biana_node_file_prefix + ".tsv"
	network_file_identifier_type = "uniprotentry"
	node_file = association_scores_file  # association_dir + "Apoptosis_uniProt_nodes.list" 
	if PPI == "piana_joan_exp":
	    network_file = network_dir + "piana_exp_human.sif"
	elif PPI == "piana_joan_all":
	    network_file = network_dir + "piana_all_human.sif"
	network_file_filtered = network_file
	interaction_relevance_file = None 
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
    elif SCORING == "nl":
	title += " - t%f" % N_LINKER_THRESHOLD
	output_dir = output_dir + "t%f" % N_LINKER_THRESHOLD + os.sep
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
    job_file = output_dir + "job.sh"
    output_log_file = output_base_dir_association + "README"

    # Analysis files
    predictions_file = output_dir + "predictions.txt"
    labels_file = output_dir + "labels.txt"
    r_script_file = output_dir + "results.r"
    tex_script_file = output_dir + "results.tex"

    score_xval_commands = { "ns": Template(src_dir + "scoreNetwork/scoreN -s s -n %s.$fold -e %s -o %s.$fold -r %d -i %d &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, score_log_file)),
			    "nz": Template(src_dir + "scoreNetwork/scoreN -s z -n %s.$fold -e %s -o %s.$fold -i %d -x %d -d %s &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file)),
			    "nh": Template(src_dir + "scoreNetwork/scoreN -s h -n %s.$fold -e %s -o %s.$fold -r %d -i %d -x %d -d %s &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file)),
			    "n1": Template(src_dir + "scoreNetwork/scoreN -s 1 -n %s.$fold -e %s -o %s.$fold -r %d -i %d -x %d -d %s &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file)),
			    "nd": Template(src_dir + "scoreNetwork/scoreN -s d -n %s.$fold -e %s.$fold -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_as_node_scores_file, output_scores_file, score_log_file)),
			    "nw": Template(src_dir + "scoreNetwork/scoreN -s w -n %s.$fold -e %s -o %s.$fold -t %f &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, DEFAULT_NON_SEED_SCORE, score_log_file)),
			    "nl": Template(src_dir + "scoreNetwork/scoreN -s l -n %s.$fold -e %s -o %s.$fold -t %f &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_LINKER_THRESHOLD, score_log_file)),
			    "nr": Template(src_dir + "scoreNetwork/scoreN -s r -n %s.$fold -e %s -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file)),
			    "nr": Template(src_dir + "scoreNetwork/scoreN -s r -n %s.$fold -e %s -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file)),
			    "nx": Template(src_dir + "scoreNetwork/scoreN -s x -n %s.$fold -e %s -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file)),
			    "nb": Template(src_dir + "./netscore -c %s.$fold -i %s -o %s.$fold -t 0 -z 0 -nr 100 -r 1 -zp 0 -n %d -nd 2 -mx 1 -ms 3 -mn 0 -dn 2 -de 2 -mxe 0 -mne 0.00000001 -mnd 0.0000001 -mnde 0.0000001 -mnst 20 -mnste 20 -dxi 1 -dxn 0 -dxe 0 -e 0.0000001 &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, score_log_file)),
			    "ff": Template(src_dir + "./fFlow %s.$fold %s %s.$fold %d %f &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, DEFAULT_NON_SEED_SCORE, score_log_file)),
			  }

    score_commands = { "ns": src_dir + "scoreNetwork/scoreN -s s -n %s -e %s -o %s -r %d -i %d &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, score_log_file),
		       "nz": src_dir + "scoreNetwork/scoreN -s z -n %s -e %s -o %s -i %d -x %d -d %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file), 
		       "nh": src_dir + "scoreNetwork/scoreN -s h -n %s -e %s -o %s -r %d -i %d -x %d -d %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file), 
		       "n1": src_dir + "scoreNetwork/scoreN -s 1 -n %s -e %s -o %s -r %d -i %d -x %d -d %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file), 
		       "nd": src_dir + "scoreNetwork/scoreN -s d -n %s -e %s -o %s &> %s" % (node_scores_file, edge_scores_as_node_scores_file, output_scores_file, score_log_file),
		       "nw": src_dir + "scoreNetwork/scoreN -s w -n %s -e %s -o %s -t %f &> %s" % (node_scores_file, edge_scores_file, output_scores_file, DEFAULT_NON_SEED_SCORE, score_log_file),
		       "nl": src_dir + "scoreNetwork/scoreN -s l -n %s -e %s -o %s -t %f &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_LINKER_THRESHOLD, score_log_file),
		       "nr": src_dir + "scoreNetwork/scoreN -s r -n %s -e %s -o %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file),
		       "nx": src_dir + "scoreNetwork/scoreN -s x -n %s -e %s -o %s &> %s" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file),
		       "nb": src_dir + "./netscore -c %s -i %s -o %s -t 0 -z 0 -nr 100 -r 1 -zp 0 -n %d -nd 2 -mx 1 -ms 3 -mn 0 -dn 2 -de 2 -mxe 0 -mne 0.00000001 -mnd 0.0000001 -mnde 0.0000001 -mnst 20 -mnste 20 -dxi 1 -dxn 0 -dxe 0 -e 0.0000001 &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, score_log_file),
		       "ff": src_dir + "./fFlow %s %s %s %d %f &> %s" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, DEFAULT_NON_SEED_SCORE, score_log_file),
		     }


    # Conduct experiment
    if MODE == "prepare":
	prepare(PPI, ASSOCIATION, biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, biana_network_file_filtered_by_reliability, network_file, network_file_filtered, input_log_file, node_file, seed_scores_file, network_file_identifier_type, node_description_file, association_scores_file, association_scores_file_identifier_type, input_dir, node_scores_file, edge_scores_file, interaction_relevance_file, interaction_relevance_file2, edge_scores_as_node_scores_file, sampled_file_prefix)
    elif MODE == "score":
	score(SCORING, score_commands, score_xval_commands, output_scores_file, log_file, job_file)
    elif MODE == "analyze":
	analyze(PPI, output_scores_file, log_file, node_scores_file, association_scores_file_identifier_type, node_description_file, network_file_identifier_type, association_scores_validation_file, r_script_file, predictions_file, labels_file, tex_script_file, output_log_file, output_dir, title)
    elif MODE == "all":
	prepare(PPI, ASSOCIATION, biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, biana_network_file_filtered_by_reliability, network_file, network_file_filtered, input_log_file, node_file, seed_scores_file, network_file_identifier_type, node_description_file, association_scores_file, association_scores_file_identifier_type, input_dir, node_scores_file, edge_scores_file, interaction_relevance_file, interaction_relevance_file2, edge_scores_as_node_scores_file, sampled_file_prefix)
	score(SCORING, score_commands, score_xval_commands, output_scores_file, log_file, job_file)
	analyze(PPI, output_scores_file, log_file, node_scores_file, association_scores_file_identifier_type, node_description_file, network_file_identifier_type, association_scores_validation_file, r_script_file, predictions_file, labels_file, tex_script_file, output_log_file, output_dir, title)
    else:
	raise ValueError("Unrecognized mode!")
    return


def prepare(PPI, ASSOCIATION, biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, biana_network_file_filtered_by_reliability, network_file, network_file_filtered, input_log_file, node_file, seed_scores_file, network_file_identifier_type, node_description_file, association_scores_file, association_scores_file_identifier_type, input_dir, node_scores_file, edge_scores_file, interaction_relevance_file, interaction_relevance_file2, edge_scores_as_node_scores_file, sampled_file_prefix):
    """
	Creates necessary files for scoring
    """
    # Create PPI network if necessary
    #if PPI.startswith("biana"):
    #	create_biana_network_files(biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, biana_network_file_filtered_by_reliability)

    # Filter network by degree and get largest connected component
    if not os.path.exists(network_file_filtered): 
	print "Filtering by degree", network_file, "->", network_file_filtered
	prepare_data.analyze_network(network_file)
	prepare_data.create_degree_filtered_network_file(network_file, network_file_filtered, ALLOWED_MAX_DEGREE)
    prepare_data.analyze_network(network_file_filtered, out_file = input_log_file)

    seed_to_score = None
    # Get node to association mapping
    if PPI.startswith("ori"):
	os.system("awk '{ if($2 == 1) print $1, $2}' %s > %s" % (node_file, seed_scores_file))
    elif PPI.startswith("piana_joan"): # or PPI.startswith("david"):
	#os.system("awk '{print $1, 1}' %s > %s" % (node_file, seed_scores_file))
    	all_nodes = set(prepare_data.get_nodes_in_network(network_file_filtered))
    	seed_nodes = prepare_data.get_nodes_from_nodes_file(node_file)
    	seed_to_score = dict([(node, 1) for node in seed_nodes])
    	prepare_data.create_node_scores_file(nodes = (all_nodes & seed_nodes), node_to_score = seed_to_score, node_scores_file = seed_scores_file, ignored_nodes = None, default_score = DEFAULT_NON_SEED_SCORE)

    if not os.path.exists(seed_scores_file): 
	seed_to_score = prepare_data.get_node_association_score_mapping(network_file = network_file_filtered, network_file_identifier_type = network_file_identifier_type, node_description_file = node_description_file, association_scores_file = association_scores_file, association_scores_file_identifier_type = association_scores_file_identifier_type, log_file = input_log_file, default_seed_score=DEFAULT_SEED_SCORE)

    # Create initial data analysis file
    if not os.path.exists(input_dir + "analyze_network.r") and seed_to_score is not None:
	prepare_data.create_R_analyze_network_script(network_file_filtered, seeds=seed_to_score.keys(), out_path=input_dir, title = PPI + "_" + ASSOCIATION)
	os.system("R CMD BATCH %s" % (input_dir + "analyze_network.r"))
	os.system("convert %sanalyze_network.eps %sanalyze_network.jpg" % (input_dir, input_dir))
	os.system("R CMD BATCH %s" % (input_dir + "analyze_network_log_scaled.r"))
	os.system("convert %sanalyze_network_log_scaled.eps %sanalyze_network_log_scaled.jpg" % (input_dir, input_dir))
	prepare_data.create_ARFF_network_metrics_file(network_file_filtered, seed_to_score, seed_to_score.keys(), input_dir + "analyze_network.arff")

    # Prepare scoring files
    prepare_scoring_files(seed_scores_file, network_file_filtered, seed_to_score, node_scores_file, edge_scores_file, interaction_relevance_file, interaction_relevance_file2, edge_scores_as_node_scores_file, sampled_file_prefix)
    return


def score(SCORING, score_commands, score_xval_commands, output_scores_file, log_file, job_file):
    """
	Runs or prints commands to run scoring method on the input files
    """
    score_original(SCORING, score_commands, output_scores_file, log_file, job_file)
    score_xval(SCORING, score_xval_commands, output_scores_file, log_file, job_file)
    return


def analyze(PPI, output_scores_file, log_file, node_scores_file, association_scores_file_identifier_type, node_description_file, network_file_identifier_type, association_scores_validation_file, r_script_file, predictions_file, labels_file, tex_script_file, output_log_file, output_dir, title):
    """
	Does cross validation and percentage analysis on the output files
    """
    analyzed = analyze_xval(r_script_file, output_scores_file, node_scores_file, predictions_file, labels_file, tex_script_file, output_log_file, output_dir, title, log_file) 
    if analyzed:
	analyze_xval_percentage(log_file, output_scores_file, node_scores_file, output_log_file)
    #analyze_original(PPI, output_scores_file, log_file, node_scores_file, association_scores_file_identifier_type, node_description_file, network_file_identifier_type, association_scores_validation_file)
    return


def generate_score_xval_command(SCORING, score_xval_commands, k):
    #return score_xval_command % (node_scores_file, k, edge_scores_file, output_scores_file, k, N_REPETITION, N_ITERATION, score_log_file, k)
    return score_xval_commands[SCORING].substitute(fold = '%d' % k)


def score_xval(SCORING, score_xval_commands, output_scores_file, log_file, job_file):
    f = open(log_file, "a")
    for k in range(1, N_X_VAL+1):
	if not os.path.exists(output_scores_file + ".%d" % k):
	    if only_print_command:
		print generate_score_xval_command(SCORING, score_xval_commands, k)
	    else:
		f.write("%s\n" % generate_score_xval_command(SCORING, score_xval_commands, k))
		if use_cluster:
		    #f2 = NamedTemporaryFile(delete=False)
		    #f2 = open("%s.%d" % (job_file, k), 'w') 
		    #f2.write("#!/usr/bin/env sh\n%s" % generate_score_xval_command(SCORING, score_xval_commands, k))
		    #f2.close()
		    #os.system( "qsub -o out.%d -e err.%d -l hostname=node52 -N %s -b y %s" % (k, k, SCORING, f2.name) )
		    #os.system( "qsub -cwd -o out.%d -e err.%d -l hostname=node52 -N %s -b y %s" % (k, k, SCORING, generate_score_xval_command(SCORING, score_xval_commands, k)) )
		    #if k % 2 != 0:
		    #	continue
		    os.system( "qsub -cwd -o out.%d -e err.%d -q bigmem -N %s -b y %s" % (k, k, SCORING, generate_score_xval_command(SCORING, score_xval_commands, k)) )
		    #os.unlink(f2.name)
		else:
		    os.system( generate_score_xval_command(SCORING, score_xval_commands, k) )
    f.close()
    return


def score_original(SCORING, score_commands, output_scores_file, log_file, job_file):
    if not os.path.exists(output_scores_file):
	f = open(log_file, "a")
	#print score_commands[SCORING]
	if only_print_command:
	    print score_commands[SCORING]
	else:
	    f.write("%s\n" % score_commands[SCORING])
	    if use_cluster:
		#f2 = open(job_file, 'w') 
		#f2.write("#!/usr/bin/env sh\n%s" % score_commands[SCORING])
		#f2.close()
		#os.system("qsub -o out -e err -l hostname=node34 -N %s -b y %s" % (SCORING, f2.name))
		#os.system("qsub -cwd -o out -e err -l hostname=node34 -N %s -b y %s" % (SCORING, score_commands[SCORING]))
		os.system("qsub -cwd -o out -e err -q bigmem -N %s -b y %s" % (SCORING, score_commands[SCORING]))
	    else:
		os.system(score_commands[SCORING])
	f.close()
    return


def analyze_xval(r_script_file, output_scores_file, node_scores_file, predictions_file, labels_file, tex_script_file, output_log_file, output_dir, title, log_file):
    if os.path.exists(r_script_file):
	return False
    list_node_scores_and_labels = []
    for k in range(1, N_X_VAL+1):
	node_validation_data = analyze_results.get_validation_node_scores_and_labels(file_result = output_scores_file+".%d"%k, file_seed_test_scores = node_scores_file+".%d.test"%k, file_node_scores = node_scores_file, n_random_negative_folds = N_RANDOM_NEGATIVE_FOLDS, default_score = DEFAULT_NON_SEED_SCORE, replicable = REPLICABLE)
	list_node_scores_and_labels.append(node_validation_data)
    analyze_results.create_ROCR_files(list_node_scores_and_labels, predictions_file, labels_file)
    analyze_results.create_R_script(r_script_file, output_dir, title) # os.path.basename(output_dir))
    os.system("R CMD BATCH %s" % (r_script_file))
    #os.system("convert %sperformance.eps %sperformance.jpg" % (output_dir, output_dir))
    #analyze_results.create_tex_script(tex_script_file, output_dir, title)
    analyze_results.record_performance_AUC_in_log_file(output_dir, output_log_file, title)

    #! now unnecessary since ROC curve analysis is done
    threshold_analysis = False
    if threshold_analysis:
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
    return True


def analyze_xval_percentage(log_file, output_scores_file, node_scores_file, output_log_file):
    f = open(log_file, "a")
    coverages = []
    for percentage in (10, 25, 50):
	#print "---- %s%%:" % percentage
	f.write("---- %s%%:\n" % percentage)
	n_seed_sum = 0.0
	n_seed_all_sum = 0.0
	values = []
	for k in range(1, N_X_VAL+1):
	    ##print output_scores_file+".%s.%d" % (SCORING, k), node_scores_file+".%d.test"%k 
	    n_seed, n_seed_all, n, i = analyze_results.calculate_seed_coverage_at_given_percentage(output_scores_file + ".%d" % k, node_scores_file+".%d.test" % k, percentage, DEFAULT_NON_SEED_SCORE)
	    ##print n_seed, n, i
	    if i > n+1:
		#print "Warning: Real coverage percentage is disputed due to equal scores!", i, n
		f.write("Warning: Real coverage percentage is disputed due to equal scores! %i %i\n" % (i, n))
	    n_seed_sum += n_seed
	    n_seed_all_sum += n_seed_all
	    values.append(float(n_seed)/n_seed_all)
	coverages.append("%f\t+/- %f" % calculate_mean_and_sigma.calc_mean_and_sigma(values))
	n_seed = n_seed_sum / N_X_VAL
	n_seed_all = n_seed_all_sum / N_X_VAL
	#print "Avg:", n_seed, "over", n
	f.write("Avg: %i over %i\n" % (n_seed, n))
	#coverages.append(float(n_seed)/n_seed_all)
    analyze_results.record_performance_coverage_in_log_file(output_log_file, coverages)
    f.close()
    return


def analyze_original(PPI, output_scores_file, log_file, node_scores_file, association_scores_file_identifier_type, node_description_file, network_file_identifier_type, association_scores_validation_file):
    if not os.path.exists(output_scores_file):
	raise Exception("Output score file does not exist!")
    f = open(log_file, "a")
    for percentage in (10, 25, 50):
	#print "---- %s:" % percentage
	f.write("---- %s:\n" % percentage)
	f.write("%s\n" % output_scores_file)
	coverage = analyze_results.calculate_seed_coverage_at_given_percentage(output_scores_file, node_scores_file, percentage, DEFAULT_NON_SEED_SCORE)
	f.write("%s\n" % str(coverage))
	if association_scores_file_identifier_type is not None and not os.path.exists(output_scores_file+"."+association_scores_file_identifier_type):
	    if PPI.startswith("biana"):
		prepare_data.convert_file_using_new_id_mapping(output_scores_file, node_description_file, network_file_identifier_type, "geneid", output_scores_file+".geneid")
		prepare_data.convert_file_using_new_id_mapping(output_scores_file+".geneid", gene_info_file, "geneid", association_scores_file_identifier_type, output_scores_file+"."+association_scores_file_identifier_type)
	    elif PPI.startswith("ori") or PPI.startswith("david") or PPI.startswith("piana_joan"):
		pass
	    else:
		prepare_data.convert_file_using_new_id_mapping(output_scores_file, node_description_file, network_file_identifier_type, association_scores_file_identifier_type, output_scores_file+"."+association_scores_file_identifier_type)
	if association_scores_validation_file is not None:
	    f.write("Validation seed coverage:")
	    if PPI.startswith("ori") or PPI.startswith("david") or PPI.startswith("piana_joan"):
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


def prepare_scoring_files(seed_scores_file, network_file_filtered, seed_to_score, node_scores_file, edge_scores_file, interaction_relevance_file, interaction_relevance_file2, edge_scores_as_node_scores_file, sampled_file_prefix):
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
    edges = None
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
    if not os.path.exists(edge_scores_as_node_scores_file):
	if edges is None:
	    edges = prepare_data.get_edges_in_network(network_file = network_file_filtered)
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

# obsolete
import score_network # obsolete

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
#netshort_results = output_dir + "netshort_onlyEdgeRelevance_noNodeInitialScoreAccumulation.txt"
#netrank_results = output_dir + "netrank_onlyEdgeRelevance_noNodeInitialScoreAccumulation.txt"
#netscore_results_1 = output_dir + "netscore_onlyEdgeRelevance_noNodeInitialScoreAccumulation_i1.txt"
#netscore_results_3 = output_dir + "netscore_onlyEdgeRelevance_noNodeInitialScoreAccumulation_i3.txt"
#netscore_results_no_1 = output_dir + "netscore_noEdgeRelevance_noNodeInitialScoreAccumulation_i1.txt"
#netscore_results_no_3 = output_dir + "netscore_noEdgeRelevance_noNodeInitialScoreAccumulation_i3.txt"
#netzcore_results_1 = output_dir + "netzcore_onlyEdgeRelevance_noNodeInitialScoreAccumulation_i1.txt"
#netzcore_results_3 = output_dir + "netzcore_onlyEdgeRelevance_noNodeInitialScoreAccumulation_i3.txt"


if __name__ == "__main__":
    main()

