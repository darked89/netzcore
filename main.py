#!/usr/bin/env python

#########################################################################
# Main workflow for scoring 
#
# eg 25/06/2009
#########################################################################

import prepare_data
import analyze_results
import calculate_mean_and_sigma
import os, os.path, time, subprocess
from string import Template
from sys import stdout as sys_stdout

from globals import *


def main(ppis, phenotypes, scoring_parameters, user_friendly_id=user_friendly_id):
    """
    Execute given experiments in the form of (ppi, phenotype, scoring_parameter) tuples
    """
    experiments = []
    for ppi in ppis:
	for phenotype in phenotypes:
	    for parameters in scoring_parameters:
		experiments.append((ppi, phenotype, parameters[0], parameters[1], parameters[2]))

    if MODE == "compare":
	compare_experiments(experiments, tex_format, functional_enrichment, user_friendly_id, summary_seed_cutoff, analysis_type = "common") # "user")
    elif MODE == "summary":
	sum_up_experiments(ppis, phenotypes, "auc", tex_format, user_friendly_id, summary_seed_cutoff)
	sum_up_experiments(ppis, phenotypes, "cov", tex_format, user_friendly_id, summary_seed_cutoff)
    else:
	#experiment_count = 0
	for experiment in experiments:
	    PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION = experiment
	    print "Running experiment:", experiment
	    if ignore_experiment_failures:
		try:
		    run_experiment(MODE, PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION, functional_enrichment, analyze_network, prepare_mutated)
		except:
		    print "!Problem!"
	    else:
		run_experiment(MODE, PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION, functional_enrichment, analyze_network, prepare_mutated)
	    #experiment_count += 1
	    if MODE in ("score", "all") and use_cluster and delay_experiment:
		delay = 10
		experiment_count = get_number_of_jobs_in_queues()
		while experiment_count > 20: #!
		    time.sleep(delay)
		    experiment_count = get_number_of_jobs_in_queues()
    
    #if MODE == "prepare" and use_cluster == False:
    #	os.system("scp -r ../data/input gaudi:netzcore/data/")
    return


def get_number_of_jobs_in_queues():
    p1 = subprocess.Popen(["qstat"], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
    experiment_count = int(p2.communicate()[0])
    return experiment_count


def decide_association_data(ASSOCIATION, PPI):
    """
	Decide disease association files (Association data to be used)
    """
    association_scores_validation_file = None
    candidates_file = None
    if ASSOCIATION == "custom": 
	association_dir = None 
	association_scores_file = None
	association_scores_file_identifier_type = None
    elif ASSOCIATION == "aneurysm":
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
    elif ASSOCIATION.startswith("chen_"):
	association_dir = data_dir + "chen_disease_data" + os.sep
	association_scores_file = association_dir + ASSOCIATION + ".txt"
	association_scores_file_identifier_type = "genesymbol"
    elif ASSOCIATION.startswith("navlakha_"):
	association_dir = data_dir + "navlakha" + os.sep 
	association_scores_file = association_dir + "associations" + os.sep + ASSOCIATION + ".txt"
	candidates_file = association_dir + "candidates" + os.sep + ASSOCIATION + ".txt"
	association_scores_file_identifier_type = "genesymbol"
    elif ASSOCIATION.startswith("hsdl_"):
	association_dir = data_dir + "tf_lineage" + os.sep + "hsdl_classification" + os.sep 
	association_scores_file = association_dir + ASSOCIATION + ".txt"
	association_scores_file_identifier_type = "genesymbol"
    elif ASSOCIATION.startswith("rob_"):
	association_dir = data_dir + "rob" + os.sep + "gene_sets" + os.sep 
	association_scores_file = association_dir + ASSOCIATION + ".txt"
	association_scores_file_identifier_type = "genesymbol"
    elif ASSOCIATION == "baldo_synthetic":
	association_dir = data_dir + "baldo_synthetic" + os.sep 
	association_scores_file = association_dir + "seeds.csv"
	association_scores_file_identifier_type = None
    elif ASSOCIATION.startswith("perturbed_omim_"):
	params = ASSOCIATION[len("perturbed_omim_"):]
	idx = params.rindex("_")
	i = params[idx+1:]
	idx2 = params.rindex("_", 0, idx)
	p = params[idx2+1:idx]
	p = p.lstrip("p")
	assoc = "omim_" + params[:idx2]
	#p, i = params[idx+1:].split("_")
	#print params, assoc
	#print p,i
	if PPI == "biana_no_tap_no_reliability":
	    association_dir = data_dir + "human_interactome_biana" + os.sep + "perturbed" + os.sep + assoc + os.sep + str(p) + os.sep
	    association_scores_file = association_dir + "sampled.txt." + str(i)
	    association_scores_file_identifier_type = None 
	elif PPI == "goh":
	    association_dir = data_dir + "goh_human_ppi" + os.sep + "perturbed" + os.sep + assoc + os.sep + str(p) + os.sep
	    association_scores_file = association_dir + "sampled.txt." + str(i)
	    association_scores_file_identifier_type = "genesymbol"
	else:
	    raise ValueError("Unrecognized ppi for perturbed association!")
    else:
	raise ValueError("Unrecognized association!")
    return (association_scores_file, association_scores_file_identifier_type, association_scores_validation_file, candidates_file, association_dir)


def decide_interaction_data(PPI, ASSOCIATION, association_scores_file):
    """ 
	Decide interaction data to be used
    """
    # Network specific
    global DEFAULT_NON_SEED_SCORE 
    specie = "Homo sapiens" # For GO functional enrichment analysis
    interaction_relevance_file = None
    interaction_relevance_file2 = None
    biana_network_file_filtered_by_method = None
    biana_network_file_filtered_by_reliability = None
    node_file = None
    node_description_file = None 
    network_dir = None
    # BIANA ppi
    if PPI.startswith("biana"):
	network_dir = biana_network_dir
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
	elif PPI == "biana_coexpression":
	    network_file = biana_network_file_prefix + "_coexpression_filtered.sif"
	elif PPI == "biana_coexpression_relevance":
	    network_file = biana_network_file_prefix + "_coexpression_filtered.sif"
	    interaction_relevance_file = biana_network_file_prefix + "_stringscore_coexpression.eda"
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
	    elif PPI.startswith("biana_no_tap_no_reliability_permuted_"):
		node_description_file = None 
		node_file = data_dir + "input" + os.sep + "biana_no_tap_no_reliability" + os.sep + ASSOCIATION + os.sep + "seed_scores.sif"
		params = PPI[len("biana_no_tap_no_reliability_permuted_"):]
		p, i = params.split("_")
		p = p.lstrip("p")
		network_dir = data_dir + "human_interactome_biana" + os.sep + "permuted" + os.sep + str(p) + os.sep
		network_file = network_dir + "sampled_graph.sif." + str(i) 
	    elif PPI.startswith("biana_no_tap_no_reliability_pruned_non_seed_interactions_"):
		node_description_file = None 
		node_file = data_dir + "input" + os.sep + "biana_no_tap_no_reliability" + os.sep + ASSOCIATION + os.sep + "seed_scores.sif"
		params = PPI[len("biana_no_tap_no_reliability_pruned_non_seed_interactions_"):]
		p, i = params.split("_")
		p = p.lstrip("p")
		network_dir = data_dir + "human_interactome_biana" + os.sep + "pruned_non_seed_interactions" + os.sep + ASSOCIATION + os.sep + str(p) + os.sep
		network_file = network_dir + "sampled_graph.sif." + str(i) 
	    elif PPI.startswith("biana_no_tap_no_reliability_pruned_"):
		node_description_file = None 
		node_file = data_dir + "input" + os.sep + "biana_no_tap_no_reliability" + os.sep + ASSOCIATION + os.sep + "seed_scores.sif"
		params = PPI[len("biana_no_tap_no_reliability_pruned_"):]
		p, i = params.split("_")
		p = p.lstrip("p")
		network_dir = data_dir + "human_interactome_biana" + os.sep + "pruned" + os.sep + ASSOCIATION + os.sep + str(p) + os.sep
		network_file = network_dir + "sampled_graph.sif." + str(i) 
	    else:
		raise ValueError("Unrecognized ppi!")
	else:
	    raise ValueError("Unrecognized ppi!")
	if ASSOCIATION.startswith("perturbed_"):
	    node_description_file = None 
	    node_file = association_scores_file
	    network_dir = data_dir + "input" + os.sep + PPI + os.sep
	    network_file = network_dir + "edge_scores.sif"
	if PPI.startswith("biana_no_tap_no_reliability_permuted_") or PPI.startswith("biana_no_tap_no_reliability_pruned_non_seed_interactions_") or PPI.startswith("biana_no_tap_no_reliability_pruned_"):
	    network_file_filtered = network_file
	else:
	    network_file_filtered = network_file[:-4] + "_degree_filtered.sif" # Using only the largest strongly connected component
    # Goh ppi
    elif PPI.startswith("goh"):
	network_file_identifier_type = "geneid"
	node_description_file = gene_info_file 
	if PPI == "goh":
	    network_dir = goh_network_dir
	    network_file = network_dir + "ppi.sif"
	    network_file_filtered = network_file[:-4] + "_degree_filtered.sif" 
	elif PPI == "goh_1e5":
	    network_file = goh_network_dir + "ppi.sif"
	    network_file_filtered = network_file[:-4] + "_degree_filtered.sif" 
	    DEFAULT_NON_SEED_SCORE = 0.00001 
	elif PPI.startswith("goh_permuted_"):
	    node_description_file = None 
	    node_file = data_dir + "input" + os.sep + "goh" + os.sep + ASSOCIATION + os.sep + "seed_scores.sif"
	    params = PPI[len("goh_permuted_"):]
	    p, i = params.split("_")
	    p = p.lstrip("p")
	    network_dir = goh_network_dir + "permuted" + os.sep + str(p) + os.sep
	    network_file = network_dir + "sampled_graph.sif." + str(i) 
	    network_file_filtered = network_file
	elif PPI.startswith("goh_pruned_"):
	    node_description_file = None 
	    node_file = data_dir + "input" + os.sep + "goh" + os.sep + ASSOCIATION + os.sep + "seed_scores.sif"
	    params = PPI[len("goh_pruned_"):]
	    p, i = params.split("_")
	    p = p.lstrip("p")
	    network_dir = goh_network_dir + "pruned" + os.sep + ASSOCIATION + os.sep + str(p) + os.sep
	    network_file = network_dir + "sampled_graph.sif." + str(i) 
	    network_file_filtered = network_file
	if ASSOCIATION.startswith("perturbed_"):
	    node_description_file = None 
	    node_file = association_scores_file
	    network_dir = data_dir + "input" + os.sep + PPI + os.sep
	    network_file = network_dir + "edge_scores.sif"
    # Entrez ppi
    elif PPI == "entrez":
	node_description_file = gene_info_file 
	network_file_identifier_type = "geneid"
	network_file = data_dir + "entrez_human_ppi" + os.sep + "ppi.sif"
	network_file_filtered = network_file[:-4] + "_degree_filtered.sif"
    # Navlakha HPRD ppi
    elif PPI == "hprd":
	node_description_file = None 
	network_file_identifier_type = "genesymbol"
	network_file = data_dir + "navlakha" + os.sep + "hprd.sif"
	network_file_filtered = network_file[:-4] + "_degree_filtered.sif"
	node_file = association_scores_file
    # Navlakha OPHID ppi
    elif PPI == "ophid":
	node_description_file = None 
	network_file_identifier_type = "genesymbol"
	network_file = data_dir + "navlakha" + os.sep + "ophid.sif"
	network_file_filtered = network_file[:-4] + "_degree_filtered.sif"
	node_file = association_scores_file
    # Rivasi ppi
    elif PPI == "rivasi":
	node_description_file = gene_info_file 
	network_file_identifier_type = "geneid"
	network_file = data_dir + "tf_lineage" + os.sep + "ppi.sif"
	network_file_filtered = network_file[:-4] + "_degree_filtered.sif"
    # Rhodes ppi
    elif PPI == "rhodes":
	network_dir = rhodes_network_dir
	node_description_file = gene_info_file 
	network_file_identifier_type = "geneid"
	network_file = network_dir + "ppi.sif"
	network_file_filtered = network_file[:-4] + "_degree_filtered.sif" 
	#interaction_relevance_file = network_file[:-4] + ".eda" # need to rescale / cluster scores because max/min >= 10000
    elif PPI.startswith("ori"):
	if PPI.startswith("ori_no_tap"):
	    network_base_dir = data_dir + "human_interactome_ori_no_tap" + os.sep
	else:
	    network_base_dir = data_dir + "human_interactome_ori" + os.sep
	node_description_file = None #biana_node_file_prefix + ".tsv"
	network_file_identifier_type = "user entity id"
	if PPI == "ori_network" or PPI == "ori_no_tap_network":
	    network_dir = network_base_dir + "network" + os.sep
	    node_file = network_dir + "aneurist.noa"
	    network_file = network_dir + "aneurist.sif"
	    network_file_filtered = network_file
	    interaction_relevance_file = network_dir + "aneurist.eda"
	    DEFAULT_NON_SEED_SCORE = 0.00001 
	elif PPI == "ori_coexpression_1e-2" or PPI == "ori_no_tap_coexpression_1e-2":
	    network_dir = network_base_dir + "coexpression" + os.sep
	    node_file = network_dir + "aneurist.noa"
	    network_file = network_dir + "aneurist.sif"
	    network_file_filtered = network_file
	    interaction_relevance_file = network_dir + "aneurist.eda"
	elif PPI == "ori_coexpression" or PPI == "ori_no_tap_coexpression":
	    network_dir = network_base_dir + "coexpression" + os.sep
	    node_file = network_dir + "aneurist.noa"
	    network_file = network_dir + "aneurist.sif"
	    network_file_filtered = network_file
	    interaction_relevance_file = network_dir + "aneurist.eda"
	    DEFAULT_NON_SEED_SCORE = 0.00001 
	elif PPI == "ori_coexpression_colocalization_1e-2" or PPI == "ori_no_tap_coexpression_colocalization_1e-2":
	    network_dir = network_base_dir + "coexpression_colocalization" + os.sep
	    node_file = network_dir + "aneurist.noa"
	    network_file = network_dir + "aneurist.sif"
	    network_file_filtered = network_file
	    interaction_relevance_file = network_dir + "aneurist.eda"
	elif PPI == "ori_coexpression_colocalization" or PPI == "ori_no_tap_coexpression_colocalization":
	    network_dir = network_base_dir + "coexpression_colocalization" + os.sep
	    node_file = network_dir + "aneurist.noa"
	    network_file = network_dir + "aneurist.sif"
	    network_file_filtered = network_file
	    interaction_relevance_file = network_dir + "aneurist.eda"
	    DEFAULT_NON_SEED_SCORE = 0.00001 
	elif PPI == "ori_colocalization" or PPI == "ori_no_tap_colocalization":
	    network_dir = network_base_dir + "colocalization" + os.sep
	    node_file = network_dir + "aneurist.noa"
	    network_file = network_dir + "aneurist.sif"
	    network_file_filtered = network_file
	    interaction_relevance_file = network_dir + "aneurist.eda"
	    DEFAULT_NON_SEED_SCORE = 0.00001 
    elif PPI == "javi":
	network_dir = data_dir + "human_interactome_javi" + os.sep
	node_description_file = None
	network_file_identifier_type = "userentityid"
	node_file = network_dir + "root_set.txt"
	network_file = network_dir + "human_network.tab"
	network_file_filtered = network_file
	#DEFAULT_NON_SEED_SCORE = 0.00001 
    elif PPI == "david": 
	network_dir = data_dir + "human_interactome_david" + os.sep
	network_file_identifier_type = "user entity id"
	network_file = network_dir + "human_network.sif"
	node_description_file = network_dir + "human_nodes.tsv"
	network_file_filtered = network_file[:-4] + "_degree_filtered.sif"
    elif PPI.startswith("david_"):
	network_dir = data_dir + "alzheimer_david10_Sep" + os.sep
	network_file_identifier_type = "user entity id"
	if PPI.startswith("david_homology"):
	    network_file = network_dir + "networks" + os.sep + "alzheimer_homology_network.sorted.sif"
	else:
	    network_file = network_dir + "networks" + os.sep + "alzheimer_human_network.sorted.sif"
	network_file_filtered = network_file[:-4] + "_degree_filtered.sif"
	type = PPI.split("_")[-1]
	global ONLY_LARGEST_COMPONENT 
	ONLY_LARGEST_COMPONENT = False
	#association_scores_file = data_dir + os.sep + "seeds" + "alzheimer_%s.list" % type
	node_file = network_dir + "seeds" + os.sep + "alzheimer.%s.ueID.list" % type
	#node_description_file = network_dir + "alzheimer_%s.tab" % type
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
    elif PPI.startswith("biogrid_yeast"):
	network_dir = data_dir + "rob" + os.sep
	node_description_file = None
	network_file_identifier_type = "orfname"
	node_file = association_scores_file
	network_file = network_dir + PPI + ".sif"
	network_file_filtered = network_file
    elif PPI == "yeastnet2":
	network_dir = data_dir + "rob" + os.sep
	node_description_file = None
	network_file_identifier_type = "orfname"
	node_file = association_scores_file
	network_file = network_dir + "yeastnet2.sif"
	interaction_relevance_file = network_dir + "yeastnet2.eda"
	network_file_filtered = network_file
    elif PPI == "baldo_synthetic":
	network_dir = data_dir + "baldo_synthetic" + os.sep
	node_description_file = None
	network_file_identifier_type = None
	node_file = association_scores_file
	network_file = network_dir + "network.sif"
	network_file_filtered = network_file
    else:
	raise ValueError("Unrecognized ppi!")

    return (interaction_relevance_file, interaction_relevance_file2, biana_network_file_filtered_by_method, \
	    biana_network_file_filtered_by_reliability, node_file, node_description_file, \
	    network_file_identifier_type, network_file, network_file_filtered, specie, network_dir)


def decide_title(PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION, N_LINKER_THRESHOLD):
    """
	Decide title for the run
    """
    title = "%s - %s - %s" % (PPI, ASSOCIATION, SCORING)

    if SCORING == "ns" or SCORING == "nh" or SCORING == "n1":
	title += " - r%d i%d" % (N_REPETITION, N_ITERATION)
    elif SCORING == "nz" or SCORING == "ff" or SCORING == "nb":
	title += " - i%d" % N_ITERATION
    elif SCORING == "nl":
	title += " - t%f" % N_LINKER_THRESHOLD

    return title


def decide_directory_hierarchy(PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION, N_LINKER_THRESHOLD):
    """
	Decide directory structure and create necessary folders
    """
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
    summary_dir = data_dir + "summary" + os.sep
    compare_dir = data_dir + "compare" + os.sep

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
    if not os.path.exists(summary_dir): 
	os.mkdir(summary_dir)
    if not os.path.exists(compare_dir): 
	os.mkdir(compare_dir)

    if SCORING == "ns" or SCORING == "nh" or SCORING == "n1":
	output_dir = output_dir + "r%di%d" % (N_REPETITION, N_ITERATION) + os.sep
	if not os.path.exists(output_dir): 
	    os.mkdir(output_dir)
    elif SCORING == "nz" or SCORING == "ff" or SCORING == "nb":
	output_dir = output_dir + "i%d" % N_ITERATION + os.sep
	if not os.path.exists(output_dir): 
	    os.mkdir(output_dir)
    elif SCORING == "nl":
	output_dir = output_dir + "t%f" % N_LINKER_THRESHOLD + os.sep
	if not os.path.exists(output_dir): 
	    os.mkdir(output_dir)

    return (input_dir, input_base_dir_network, sampling_dir, output_dir, output_base_dir_association, summary_dir, compare_dir)


def decide_scoring_and_analysis_files(input_dir, input_base_dir_network, sampling_dir, output_dir, output_base_dir_association):
    """
	Decide file names to be used in scoring and analysis
    """
    # Input/Output score file
    seed_scores_file = input_dir + "seed_scores.sif"
    node_scores_file = input_dir + "node_scores.sif"
    node_mapping_file = input_base_dir_network + "node_mapping.tsv"
    #edge_scores_file = input_dir + "edge_scores.sif"
    edge_scores_file = input_base_dir_network + "edge_scores.sif"
    edge_scores_as_node_scores_file = input_dir + "edge_scores_as_node_scores.sif"
    #reliability_filtered_edge_scores_file = input_base_dir_network + "reliability_filtered_edge_scores.sif"
    output_scores_file = output_dir + "node_scores.sif" # "_r%dn%d.sif" % (N_REPETITION, N_ITERATION)
    score_log_file = output_dir + "log.txt" # "_r%dn%d.txt" % (N_REPETITION, N_ITERATION)
    sampled_file_prefix = sampling_dir + "sampled_graph"

    # Log (README) files 
    log_file = output_dir + "README"
    input_log_file = input_dir + "README"
    output_log_file = output_base_dir_association + "README"
    job_file = output_dir + "job.sh"

    # Analysis files
    predictions_file = output_dir + "predictions.txt"
    labels_file = output_dir + "labels.txt"
    r_script_file = output_dir + "results.r"
    tex_script_file = output_dir + "results.tex"

    return (seed_scores_file, node_scores_file, node_mapping_file, edge_scores_file, edge_scores_as_node_scores_file, output_scores_file, \
	    score_log_file, sampled_file_prefix, log_file, input_log_file, job_file, output_log_file, predictions_file, \
	    labels_file, r_script_file, tex_script_file)


def decide_score_commands(node_scores_file, edge_scores_file, output_scores_file, edge_scores_as_node_scores_file, N_REPETITION, N_ITERATION, sampling_dir, score_log_file):
    """
	Decide commands to be used in scoring
    """
    score_xval_commands = { "ns": Template(src_dir + "scoreNetwork/scoreN -s s -n %s.$fold -e %s -o %s.$fold -r %d -i %d &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, score_log_file)),
			    "nz": Template(src_dir + "scoreNetwork/scoreN -s z -n %s.$fold -e %s -o %s.$fold -i %d -x %d -d %s &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file)),
			    "nh": Template(src_dir + "scoreNetwork/scoreN -s h -n %s.$fold -e %s -o %s.$fold -r %d -i %d -x %d -d %s &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file)),
			    "n1": Template(src_dir + "scoreNetwork/scoreN -s 1 -n %s.$fold -e %s -o %s.$fold -r %d -i %d -x %d -d %s &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_REPETITION, N_ITERATION, N_SAMPLE_GRAPH, sampling_dir, score_log_file)),
			    "nd": Template(src_dir + "scoreNetwork/scoreN -s d -n %s.$fold -e %s.$fold -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_as_node_scores_file, output_scores_file, score_log_file)),
			    "nw": Template(src_dir + "scoreNetwork/scoreN -s w -n %s.$fold -e %s -o %s.$fold -t %f &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, DEFAULT_NON_SEED_SCORE, score_log_file)),
			    "nl": Template(src_dir + "scoreNetwork/scoreN -s l -n %s.$fold -e %s -o %s.$fold -t %f &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_LINKER_THRESHOLD, score_log_file)),
			    "nr": Template(src_dir + "scoreNetwork/scoreN -s r -n %s.$fold -e %s -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file)), 
			    "nx": Template(src_dir + "scoreNetwork/scoreN -s x -n %s.$fold -e %s -o %s.$fold &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, score_log_file)),
			    "nb": Template(src_dir + "./netscore -c %s.$fold -i %s -o %s.$fold -t 0 -z 0 -nr 100 -r 1 -zp 0 -n %d -nd 2 -mx 1 -ms 3 -mn 0 -dn 2 -de 2 -mxe 0 -mne 0.00000001 -mnd 0.0000001 -mnde 0.0000001 -mnst 20 -mnste 20 -dxi 1 -dxn 0 -dxe 0 -e 0.0000001 &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, score_log_file)),
			    "ff": Template(src_dir + "./fFlow %s.$fold %s %s.$fold %d %f &> %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, N_ITERATION, DEFAULT_NON_SEED_SCORE, score_log_file)),
			    #"rw": Template("/sbi/users/emre/bin/R --slave --args %s.$fold %s %s.$fold < %sscoreNetwork/random_walk.r > %s.$fold" % (node_scores_file, edge_scores_file, output_scores_file, src_dir, score_log_file)),
			    "rw": Template("/sbi/users/emre/bin/R -f %sscoreNetwork/random_walk.r --args %s.$fold %s %s.$fold > %s.$fold" % (src_dir, node_scores_file, edge_scores_file, output_scores_file, score_log_file)),
			    "np": Template("/sbi/users/emre/bin/R -f %sscoreNetwork/random_walk.r --args %s.$fold %s %s.$fold 1 > %s.$fold" % (src_dir, node_scores_file, edge_scores_file, output_scores_file, score_log_file)),
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
		       "rw": "/sbi/users/emre/bin/R -f %sscoreNetwork/random_walk.r --args %s %s %s > %s" % (src_dir, node_scores_file, edge_scores_file, output_scores_file, score_log_file),
		       "np": "/sbi/users/emre/bin/R -f %sscoreNetwork/random_walk.r --args %s %s %s 1 > %s" % (src_dir, node_scores_file, edge_scores_file, output_scores_file, score_log_file),
		     }
    return score_xval_commands, score_commands


def run_experiment(MODE, PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION, functional_enrichment, analyze_network, prepare_mutated):

    # Create experiment parameters
    # Disease association files (Association data to be used)
    association_scores_file, association_scores_file_identifier_type, association_scores_validation_file, candidates_file, association_dir = decide_association_data(ASSOCIATION, PPI)

    # Interaction files (Interaction data to be used)
    (interaction_relevance_file, interaction_relevance_file2, biana_network_file_filtered_by_method, \
	    biana_network_file_filtered_by_reliability, node_file, node_description_file, \
	    network_file_identifier_type, network_file, network_file_filtered, specie, network_dir) = decide_interaction_data(PPI, ASSOCIATION, association_scores_file)

    # Human readable title for the run
    title = decide_title(PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION, N_LINKER_THRESHOLD)

    # Project directory structure
    (input_dir, input_base_dir_network, sampling_dir, output_dir, output_base_dir_association, summary_dir, compare_dir) = decide_directory_hierarchy(PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION, N_LINKER_THRESHOLD)

    # Input/Output score, logging and analysis files
    (seed_scores_file, node_scores_file, node_mapping_file, edge_scores_file, edge_scores_as_node_scores_file, output_scores_file, \
	score_log_file, sampled_file_prefix, log_file, input_log_file, job_file, output_log_file, predictions_file, \
	labels_file, r_script_file, tex_script_file) = decide_scoring_and_analysis_files(input_dir, input_base_dir_network, sampling_dir, output_dir, output_base_dir_association)

    # Scoring commands
    score_xval_commands, score_commands = decide_score_commands(node_scores_file, edge_scores_file, output_scores_file, edge_scores_as_node_scores_file, N_REPETITION, N_ITERATION, sampling_dir, score_log_file)

    global N_X_VAL, N_SEED
    # Conduct experiment
    if MODE == "prepare":
	prepare(PPI, ASSOCIATION, biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, biana_network_file_filtered_by_reliability, network_file, network_file_filtered, input_log_file, node_file, seed_scores_file, network_file_identifier_type, node_description_file, association_scores_file, association_scores_file_identifier_type, node_mapping_file, input_dir, node_scores_file, edge_scores_file, interaction_relevance_file, interaction_relevance_file2, edge_scores_as_node_scores_file, sampled_file_prefix, analyze_network, prepare_mutated, network_dir, association_dir)
    elif MODE == "score":
	N_SEED, n_linker, n_path = prepare_data.get_number_of_mapped_seeds(input_log_file)
	if leave_one_out_xval:
	    N_X_VAL = N_SEED
	score(SCORING, score_commands, score_xval_commands, output_scores_file, log_file, job_file)
    elif MODE == "analyze":
	N_SEED, n_linker, n_path = prepare_data.get_number_of_mapped_seeds(input_log_file)
	#if N_X_VAL is None:
	if leave_one_out_xval:
	    N_X_VAL = N_SEED
	analyze(PPI, SCORING, output_scores_file, log_file, node_scores_file, association_scores_file_identifier_type, node_mapping_file, node_description_file, network_file_identifier_type, association_scores_validation_file, r_script_file, predictions_file, labels_file, tex_script_file, output_log_file, output_dir, title, specie, network_file_filtered, candidates_file, functional_enrichment)
    elif MODE == "all":
	prepare(PPI, ASSOCIATION, biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, biana_network_file_filtered_by_reliability, network_file, network_file_filtered, input_log_file, node_file, seed_scores_file, network_file_identifier_type, node_description_file, association_scores_file, association_scores_file_identifier_type, node_mapping_file, input_dir, node_scores_file, edge_scores_file, interaction_relevance_file, interaction_relevance_file2, edge_scores_as_node_scores_file, sampled_file_prefix, analyze_network, prepare_mutated, network_dir, association_dir)
	score(SCORING, score_commands, score_xval_commands, output_scores_file, log_file, job_file)
	analyze(PPI, SCORING, output_scores_file, log_file, node_scores_file, association_scores_file_identifier_type, node_mapping_file, node_description_file, network_file_identifier_type, association_scores_validation_file, r_script_file, predictions_file, labels_file, tex_script_file, output_log_file, output_dir, title, specie, network_file_filtered, candidates_file, functional_enrichment)
    else:
	raise ValueError("Unrecognized mode!")
    return


def compare_experiments(experiments, tex_format=False, functional_enrichment=False, user_friendly_id="test", seed_cutoff=None, analysis_type = "common"):
    """
	Selects and checks functional annotation of common highest scoring nodes (mapping their genesymols) in different experiments
    """
    top_scoring_ids = None
    all_scoring_ids = set()
    species = set()
    prev_id_type = None
    scoring_to_values = {}
    for experiment in experiments:
	PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION = experiment
	# Disease association files (Association data to be used)
	association_scores_file, association_scores_file_identifier_type, association_scores_validation_file, candidates_file, association_dir = decide_association_data(ASSOCIATION, PPI)
	# Interaction files (Interaction data to be used)
	(interaction_relevance_file, interaction_relevance_file2, biana_network_file_filtered_by_method, \
	    biana_network_file_filtered_by_reliability, node_file, node_description_file, \
	    network_file_identifier_type, network_file, network_file_filtered, specie, network_dir) = decide_interaction_data(PPI, ASSOCIATION, association_scores_file)
	# Project directory structure
	(input_dir, input_base_dir_network, sampling_dir, output_dir, output_base_dir_association, summary_dir, compare_dir) = decide_directory_hierarchy(PPI, ASSOCIATION, SCORING, N_REPETITION, N_ITERATION, N_LINKER_THRESHOLD)
	# Input/Output score, logging and analysis files
	(seed_scores_file, node_scores_file, node_mapping_file, edge_scores_file, edge_scores_as_node_scores_file, output_scores_file, \
	score_log_file, sampled_file_prefix, log_file, input_log_file, job_file, output_log_file, predictions_file, \
	labels_file, r_script_file, tex_script_file) = decide_scoring_and_analysis_files(input_dir, input_base_dir_network, sampling_dir, output_dir, output_base_dir_association)
	if analysis_type == "common":
	    # Get high scoring node ids
	    selected_ids, all_ids = analyze_results.get_top_scoring_node_ids_at_given_percentage(output_scores_file, node_scores_file, node_mapping_file+"."+association_scores_file_identifier_type, DEFAULT_TOP_SCORING_PERCENTAGE, association_scores_file_identifier_type, default_non_seed_score = DEFAULT_NON_SEED_SCORE, exclude_seeds = True)
	    #print len(selected_ids), len(all_ids)
	    if top_scoring_ids is None:
		top_scoring_ids = set(selected_ids)
	    else:
		top_scoring_ids &= set(selected_ids)
	    all_scoring_ids |= set(all_ids)
	    species.add(specie)
	    if len(species) > 1:
		raise ValueError("Comparison can be made only on the data from same specie")
	    if prev_id_type is not None:
		if prev_id_type != association_scores_file_identifier_type:
		    raise ValueError("Comparison can be made only on the association data that provides the same nomenclature")
	    else:
		prev_id_type = association_scores_file_identifier_type
	elif analysis_type == "user": # avg ppv and sens at given thresholds
	    n_seed, n_linker, n_path = prepare_data.get_number_of_mapped_seeds(input_log_file)
	    if leave_one_out_xval and n_seed == 1: 
		#print "Skipping", PPI, ASSOCIATION, SCORING
		continue
	    elif seed_cutoff is not None and n_seed < seed_cutoff:
		continue
	    threshold_to_values = scoring_to_values.setdefault(SCORING, {})
	    f = open(r_script_file.rsplit('.')[0] + ".dat")
	    for line in f.readlines():
		if line.startswith(" ppv sens"):
		    continue
		tScore, ppv, sens = line.strip().split()#map(float, line.strip().split())
		if ppv == 'None' or sens == 'None':
		    continue
		tScore, ppv, sens = map(float, line.strip().split())
		#threshold_to_values[tScore] = (ppv_sum + ppv, sens_sum + sens)
		threshold_to_values.setdefault(tScore, []).append((ppv, sens))
	    f.close()
	else:
	    raise ValueErro("Unrecognized analysis_type")

    compare_dir += user_friendly_id + os.sep
    if not os.path.exists(compare_dir): 
	os.mkdir(compare_dir)

    if analysis_type == "user":
	f = open(compare_dir + "results.dat", "w")
	f.write(" ppv sens\n")
	for scoring, threshold_to_values in scoring_to_values.iteritems():
	    #print scoring
	    for tScore, values in threshold_to_values.iteritems():
		ppv_values, sens_values =(zip(*values))
		ppv_mean, ppv_sigma = calculate_mean_and_sigma.calc_mean_and_sigma(ppv_values)
		sens_mean, sens_sigma = calculate_mean_and_sigma.calc_mean_and_sigma(sens_values)
		f.write("%s_%s %f %f\n" % (scoring, tScore, ppv_mean, sens_mean))
	f.close()
	return

    f = open(compare_dir + "README", "a")
    f.write("%s\n" % (experiments))
    f.close()
    #out_file = sys_stdout
    out_file = open(compare_dir + "comparison.txt", "w")
    out_file_tex = open(compare_dir + "comparison.tex", "w")
    if "-" in top_scoring_ids:
	top_scoring_ids.remove("-")
    if "-" in all_scoring_ids:
	all_scoring_ids.remove("-")
    for id in top_scoring_ids: 
	out_file_tex.write("%s\\\\\n" % id)
	out_file.write("%s\n" % id)
    out_file_tex.write("\n")
    if functional_enrichment:
	if tex_format:
	    file = out_file_tex
	else:
	    file = out_file
	file.write("\nFUNCTIONAL ENRICHMENT OF COMMON HIGH SCORING NODES IN GIVEN EXPERIMENTS\n")
	file.write("%s common gene names/ids among %s\n\n" % (len(top_scoring_ids), len(all_scoring_ids)))
	analyze_results.check_functional_enrichment(list(top_scoring_ids), list(all_scoring_ids), prev_id_type, file.write, specie = species.pop(), mode = "unordered", tex_format=tex_format)
    out_file.close()
    out_file_tex.close()
    return


def sum_up_experiments(ppis, phenotypes, analysis_type="auc", tex_format=False, user_friendly_id="test", seed_cutoff=None):
    """
	Gives an averaged performance summary of ppi and association data over different scoring methods
	analysis_type: "auc" or "cov"
    """
    phenotypes_to_skip = set()
    phenotype_to_seed_values = {}
    if seed_cutoff is not None:
	for ppi in ppis:
	    for phenotype in phenotypes:
		if phenotype in phenotypes_to_skip:
		    continue
		(input_dir, input_base_dir_network, sampling_dir, output_dir, output_base_dir_association, summary_dir, compare_dir) = decide_directory_hierarchy(ppi, phenotype, "nx", 1, 1, N_LINKER_THRESHOLD)
		(seed_scores_file, node_scores_file, node_mapping_file, edge_scores_file, edge_scores_as_node_scores_file, output_scores_file, \
		score_log_file, sampled_file_prefix, log_file, input_log_file, job_file, output_log_file, predictions_file, \
		labels_file, r_script_file, tex_script_file) = decide_scoring_and_analysis_files(input_dir, input_base_dir_network, sampling_dir, output_dir, output_base_dir_association)
		n_seed, n_linker, n_path = [None] * 3
		n_seed, n_linker, n_path = prepare_data.get_number_of_mapped_seeds(input_log_file)
		if n_seed < seed_cutoff: #if n_seed >= seed_cutoff: 
		    phenotypes_to_skip.add(phenotype)
		    print "Skipping", phenotype, n_seed
		    continue
		else:
		    phenotype_to_seed_values.setdefault(phenotype, []).append((n_seed, n_linker, n_path))
	phenotypes = list(set(phenotypes)-phenotypes_to_skip)

    ppi_phenotype_auc_container = dict([ (ppi, dict([ (phenotype, {}) for phenotype in phenotypes ])) for ppi in ppis ])
    for ppi in ppis:
	for phenotype in phenotypes:
	    # Project directory structure
	    (input_dir, input_base_dir_network, sampling_dir, output_dir, output_base_dir_association, summary_dir, compare_dir) = decide_directory_hierarchy(ppi, phenotype, "nx", 1, 1, N_LINKER_THRESHOLD)
	    # Input/Output score, logging and analysis files
	    (seed_scores_file, node_scores_file, node_mapping_file, edge_scores_file, edge_scores_as_node_scores_file, output_scores_file, \
	    score_log_file, sampled_file_prefix, log_file, input_log_file, job_file, output_log_file, predictions_file, \
	    labels_file, r_script_file, tex_script_file) = decide_scoring_and_analysis_files(input_dir, input_base_dir_network, sampling_dir, output_dir, output_base_dir_association)
	    #print output_log_file
	    f = open(output_log_file)
	    for line in f.readlines():
		words = line.strip().split("\t")
		sub_words = words[0].split('-')
		method = sub_words[2].strip()
		if len(sub_words) > 3:
		    parameter = sub_words[3].strip()
		auc = float(words[1])
		auc_dev = float(words[2].strip("+/-"))
		cov = float(words[3])
		cov_dev = float(words[4].strip("+/-"))
		if analysis_type == "auc":
		    #if auc<=0.5 and n_seed <= 20 and (method == "ns" or method == "nd"):
		    #	print ppi, phenotype, n_seed, method, auc
		    ppi_phenotype_auc_container[ppi][phenotype][method] = auc
		else:
		    ppi_phenotype_auc_container[ppi][phenotype][method] = cov
	    f.close()

    scoring_methods_real = set()
    for scoring in scoring_methods:
    	for ppi, phenotype_container in ppi_phenotype_auc_container.iteritems():
    	    for phenotype, method_to_auc in phenotype_container.iteritems():
    		scoring_methods_real.add(scoring)
    #print scoring_methods_real 

    #common_methods = [] + scoring_methods
    common_methods = list(scoring_methods_real)
    for scoring in scoring_methods:
	for ppi, phenotype_container in ppi_phenotype_auc_container.iteritems():
	    for phenotype, method_to_auc in phenotype_container.iteritems():
		if not method_to_auc.has_key(scoring):
		    #print scoring, ppi, phenotype
		    if scoring in common_methods:
			common_methods.remove(scoring)
    #print common_methods

    #i = 0
    for method in reversed(["ns", "nz", "nd", "ff", "nr", "rw", "np"]):
	if method in common_methods:
	    common_methods.remove(method)
	    common_methods.insert(0, method) #i
	    #i += 1


    summary_dir += user_friendly_id + os.sep
    if not os.path.exists(summary_dir): 
	os.mkdir(summary_dir)

    log_file = open(summary_dir + "README", "a")
    log_file.write("%s\n%s\n" % (phenotypes, ppis))
    log_file.write("Skipped %s\n" % (phenotypes_to_skip))
    log_file.close()

    #out_file = sys_stdout
    file_name = analysis_type + "_phenotypes" 
    #if file_name > 60:
    #	file_name_new =  "cropped_" + file_name[:60]
    #	print "Cropping\n", file_name, "\nto\n", file_name_new
    #	file_name = file_name_new 
    out_file = open(summary_dir + file_name + ".txt", "w")

    if analysis_type == "auc":
	out_file.write("\nAVERAGE AUC OVER DIFFERENT PHENOTYPES\n")
    else:
	#out_file.write("\nAVERAGE SEED COVERAGE AT %d%% OVER DIFFERENT PHENOTYPES\n" % DEFAULT_TOP_SCORING_PERCENTAGE)
	out_file.write("\nAVERAGE SEED COVERAGE OVER DIFFERENT PHENOTYPES\n")

    if tex_format:
	out_file.write("\n%s\n" % " & ".join(common_methods))
	data_file = open(summary_dir + file_name + ".dat", 'w')
	data_file.write("\t%s\n" % "\t".join(common_methods))
	for ppi, phenotype_container in ppi_phenotype_auc_container.iteritems():
	    ppi_methods = [ (0, 0) ]*len(common_methods)
	    for i, scoring in enumerate(common_methods):
		auc_list = []
		for phenotype, method_to_auc in phenotype_container.iteritems():
		    auc_list.append(method_to_auc[scoring])
		mean, sigma = calculate_mean_and_sigma.calc_mean_and_sigma(auc_list)
		ppi_methods[i] = (mean, sigma)
	    out_file.write("%s & %s\\\\\n" % (ppi, " & ".join([ "%.2f ($\\pm$%.2f)" % (100*m, 100*s) for m,s in ppi_methods ])))
	    data_file.write("%s\t%s\n" % (ppi, "\t".join(map(lambda x: str(100*x[0]), ppi_methods))))
	data_file.close()
    else:
	for scoring in common_methods:
	    out_file.write("\n%s\n" % scoring)
	    for ppi, phenotype_container in ppi_phenotype_auc_container.iteritems():
		auc_list = []
		for phenotype, method_to_auc in phenotype_container.iteritems():
		    auc_list.append(method_to_auc[scoring])
		mean, sigma = calculate_mean_and_sigma.calc_mean_and_sigma(auc_list)
		out_file.write("%s:\t%f\t+/- %f\n" % (ppi, mean, sigma))
    out_file.close()

    # for perturbed seeds
    file_name = analysis_type + "_perturbed" 
    data_file2 = open(summary_dir + file_name + ".dat", "w")
    data_file2.write("\t%s\n" % "\t".join(common_methods))
    for ppi, phenotype_container in ppi_phenotype_auc_container.iteritems():
	iter_to_score_to_pheno = {}
	for scoring in common_methods:
	    for phenotype, method_to_auc in phenotype_container.iteritems():
		if "perturbed" not in phenotype:
		    break
		pheno, iter = phenotype.rsplit("_", 1)
		iter_to_score_to_pheno.setdefault(iter, {}).setdefault(scoring, []).append(method_to_auc[scoring])
	for iter, score_to_pheno in iter_to_score_to_pheno.iteritems():
	    aucs = []
	    for scoring in common_methods:
		mean, sigma = calculate_mean_and_sigma.calc_mean_and_sigma(score_to_pheno[scoring])
		aucs.append((mean, sigma))
	    data_file2.write("%s_%s\t%s\n" % (ppi, iter, "\t".join(map(lambda x: str(100*x[0]), aucs))))
    data_file2.close()

    file_name = analysis_type + "_ppis" #+ "-".join(ppis)
    out_file = open(summary_dir + file_name + ".txt", "w")
    if analysis_type == "auc":
	out_file.write("\nAVERAGE AUC OVER DIFFERENT PPIS\n")
    else:
	#out_file.write("\nAVERAGE SEED COVERAGE AT %d%% OVER DIFFERENT PPIS\n" % DEFAULT_TOP_SCORING_PERCENTAGE)
	out_file.write("\nAVERAGE SEED COVERAGE OVER DIFFERENT PPIS\n") 
    if tex_format:
	out_file.write("\n%s\n" % " & ".join(common_methods))
	data_file = open(summary_dir + file_name + ".dat", 'w')
	data_file.write("\t%s\n" % "\t".join(common_methods))
	for phenotype in phenotypes:
	    phenotype_methods = [ (0, 0) ]*len(common_methods)
	    for i, scoring in enumerate(common_methods):
		auc_list = []
		for ppi, phenotype_container in ppi_phenotype_auc_container.iteritems():
		    auc_list.append(phenotype_container[phenotype][scoring])
		mean, sigma = calculate_mean_and_sigma.calc_mean_and_sigma(auc_list)
		phenotype_methods[i] = (mean, sigma)
	    out_file.write("%s & %s\\\\\n" % (phenotype, " & ".join([ "%.2f ($\\pm$%.2f)" % (100*m, 100*s) for m,s in phenotype_methods ])))
	    data_file.write("%s\t%s\n" % (phenotype, "\t".join(map(lambda x: str(100*x[0]), phenotype_methods))))
	data_file.close()
    else:
	for scoring in common_methods:
	    out_file.write("\n%s\n" % scoring)
	    for phenotype in phenotypes:
		auc_list = []
		for ppi, phenotype_container in ppi_phenotype_auc_container.iteritems():
		    auc_list.append(phenotype_container[phenotype][scoring])
		mean, sigma = calculate_mean_and_sigma.calc_mean_and_sigma(auc_list)
		out_file.write("%s:\t%f\t+/- %f\n" % (phenotype, mean, sigma))
    out_file.close()

    if seed_cutoff is not None:
	data_file = open(summary_dir + "seeds.dat", 'w')
	data_file.write("\tn_seed\tn_linker\tn_path\n")
	for phenotype, values in phenotype_to_seed_values.iteritems():
	    n_seeds, n_linkers, n_paths = zip(*values)
	    mean_seeds, sigma = calculate_mean_and_sigma.calc_mean_and_sigma(n_seeds)
	    #print phenotype, values
	    if n_linkers[0] is None: # If n_linker has not been calculated - though only checks the first 
		mean_linkers, sigma = None, None
	    else:
		mean_linkers, sigma = calculate_mean_and_sigma.calc_mean_and_sigma(n_linkers)
	    if n_paths[0] is None:
		mean_paths, sigma = None, None
	    else:
		mean_paths, sigma = calculate_mean_and_sigma.calc_mean_and_sigma(n_paths)
	    data_file.write("%s\t%s\t%s\t%s\n" % (phenotype, mean_seeds, mean_linkers, mean_paths))
	data_file.close()

    return



def prepare(PPI, ASSOCIATION, biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, biana_network_file_filtered_by_reliability, network_file, network_file_filtered, input_log_file, node_file, seed_scores_file, network_file_identifier_type, node_description_file, association_scores_file, association_scores_file_identifier_type, node_mapping_file, input_dir, node_scores_file, edge_scores_file, interaction_relevance_file, interaction_relevance_file2, edge_scores_as_node_scores_file, sampled_file_prefix, analyze_network, prepare_mutated, network_dir, association_dir):
    """
	Creates necessary files for scoring
    """
    # Create PPI network if necessary
    #if PPI.startswith("biana"):
    #	create_biana_network_files(biana_node_file_prefix, biana_network_file_prefix, biana_network_file_filtered_by_method, biana_network_file_filtered_by_reliability)

    # Filter network by degree and get largest connected component
    if not os.path.exists(network_file_filtered): 
	print "Filtering by degree", network_file, "->", network_file_filtered
	if analyze_network:
	    prepare_data.analyze_network(network_file)
	prepare_data.create_degree_filtered_network_file(network_file, network_file_filtered, ALLOWED_MAX_DEGREE, ONLY_LARGEST_COMPONENT)
    #prepare_data.analyze_network(network_file_filtered, out_file = input_log_file)

    seed_to_score = None
    # Get node to association mapping
    if node_description_file is None: 
	#if PPI.startswith("ori"):
	#    os.system("awk '{ if($2 == 1) print $1, $2}' %s > %s" % (node_file, seed_scores_file))
	#elif PPI.startswith("piana_joan") or PPI == "javi": # or PPI.startswith("david"):
	#else: 
	    #os.system("awk '{print $1, 1}' %s > %s" % (node_file, seed_scores_file))
	all_nodes = set(prepare_data.get_nodes_in_network(network_file_filtered))
	seed_nodes = prepare_data.get_nodes_from_nodes_file(node_file)
	seed_nodes = all_nodes & seed_nodes
	seed_to_score = dict([(node, 1) for node in seed_nodes])
	prepare_data.create_node_scores_file(nodes = (all_nodes & seed_nodes), node_to_score = seed_to_score, node_scores_file = seed_scores_file, ignored_nodes = None, default_score = DEFAULT_NON_SEED_SCORE)
	if input_log_file is not None:
	    f = open(input_log_file, 'a')
	    f.write("Covered gene products (seed nodes): %s among %s\n" % (len(all_nodes & seed_nodes), len(seed_nodes)))
	    f.close()
	if analyze_network:
	    prepare_data.analyze_network(network_file_filtered, out_file = input_log_file, seeds = seed_to_score.keys())


    if node_description_file is not None and not os.path.exists(seed_scores_file): 
	seed_to_score = prepare_data.get_node_association_score_mapping(network_file = network_file_filtered, network_file_identifier_type = network_file_identifier_type, node_description_file = node_description_file, association_scores_file = association_scores_file, association_scores_file_identifier_type = association_scores_file_identifier_type, log_file = input_log_file, default_seed_score=DEFAULT_SEED_SCORE)
	if analyze_network:
	    prepare_data.analyze_network(network_file_filtered, out_file = input_log_file, seeds = seed_to_score.keys())

    global N_X_VAL, N_SEED
    N_SEED, n_linker, n_path = prepare_data.get_number_of_mapped_seeds(input_log_file)
    if leave_one_out_xval:
	N_X_VAL = N_SEED

    # Create initial data analysis file
    if analyze_network and not os.path.exists(input_dir + "analyze_network.r") and seed_to_score is not None:
	prepare_data.create_R_analyze_network_script(network_file_filtered, seeds=seed_to_score.keys(), out_path=input_dir, title = PPI + "_" + ASSOCIATION)
	#os.system("R CMD BATCH %s" % (input_dir + "analyze_network.r"))
	#os.system("convert %sanalyze_network.eps %sanalyze_network.jpg" % (input_dir, input_dir))
	#os.system("R CMD BATCH %s" % (input_dir + "analyze_network_log_scaled.r"))
	#os.system("convert %sanalyze_network_log_scaled.eps %sanalyze_network_log_scaled.jpg" % (input_dir, input_dir))
	#prepare_data.create_ARFF_network_metrics_file(network_file_filtered, seed_to_score, seed_to_score.keys(), input_dir + "analyze_network.arff")

    global N_X_VAL, N_SEED
    N_SEED, n_linker, n_path = prepare_data.get_number_of_mapped_seeds(input_log_file)
    if N_X_VAL is None:
	N_X_VAL = N_SEED

    # Prepare scoring files
    prepare_scoring_files(PPI, seed_scores_file, network_file_filtered, seed_to_score, node_scores_file, association_scores_file_identifier_type, node_mapping_file, node_description_file, network_file_identifier_type, edge_scores_file, interaction_relevance_file, interaction_relevance_file2, edge_scores_as_node_scores_file, sampled_file_prefix)

    # Creating mutated (permuted/pruned) networks for further significance assessment runs
    if prepare_mutated is not None:
	g = prepare_data.get_network_as_graph(edge_scores_file, use_edge_data = True)
	if seed_to_score is None:
	    seed_to_score = prepare_data.get_node_to_score_from_node_scores_file(seed_scores_file)
	if prepare_mutated == "perturbed":
	    print "Creating perturbed seeds for network"
	    output_dir = network_dir + "perturbed" + os.sep 
	    if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	    output_dir += ASSOCIATION + os.sep
	    if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	    for percentage in xrange(10,100,10):
		output_dir_inter = output_dir + str(percentage) + os.sep 
		if not os.path.exists(output_dir_inter):
		    os.mkdir(output_dir_inter)
		output_prefix = output_dir_inter + "sampled.txt."
		if os.path.exists(output_prefix + "1"):
		    break
		prepare_data.sample_perturbed_seeds_at_percentage(seed_to_score, g.nodes(), N_SAMPLE_GRAPH, percentage, output_prefix)
	elif prepare_mutated == "permuted":
	    print "Creating permuted networks"
	    output_dir = network_dir + "permuted" + os.sep
	    if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	    for percentage in xrange(10,110,10):
		output_dir_inter = output_dir + str(percentage) + os.sep
		if not os.path.exists(output_dir_inter):
		    os.mkdir(output_dir_inter)
		output_prefix = output_dir_inter + "sampled_graph.sif."
		if os.path.exists(output_prefix + "1"):
		    break
		prepare_data.sample_permuted_network_at_percentage(g, N_SAMPLE_GRAPH, percentage, output_prefix)
	elif prepare_mutated == "pruned_non_seed_interactions":
	    print "Creating pruned_non_seed_interactions networks"
	    output_dir = network_dir + "pruned_non_seed_interactions" + os.sep
	    if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	    output_dir += ASSOCIATION + os.sep
	    if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	    for percentage in xrange(10,100,10):
		output_dir_inter = output_dir + str(percentage) + os.sep 
		if not os.path.exists(output_dir_inter):
		    os.mkdir(output_dir_inter)
		output_prefix = output_dir_inter + "sampled_graph.sif."
		if os.path.exists(output_prefix + "1"):
		    break
		prepare_data.sample_pruned_network_at_percentage(g, N_SAMPLE_GRAPH, percentage, output_prefix, reserved_nodes=set(seed_to_score.keys()))
	elif prepare_mutated == "pruned":
	    print "Creating pruned networks"
	    output_dir = network_dir + "pruned" + os.sep
	    if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	    output_dir += ASSOCIATION + os.sep
	    if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	    for percentage in xrange(10,100,10):
		output_dir_inter = output_dir + str(percentage) + os.sep 
		if not os.path.exists(output_dir_inter):
		    os.mkdir(output_dir_inter)
		output_prefix = output_dir_inter + "sampled_graph.sif."
		if os.path.exists(output_prefix + "1"):
		    break
		prepare_data.sample_pruned_network_at_percentage(g, N_SAMPLE_GRAPH, percentage, output_prefix, reserved_nodes=None)
    return


def score(SCORING, score_commands, score_xval_commands, output_scores_file, log_file, job_file):
    """
	Runs or prints commands to run scoring method on the input files
    """
    if N_SEED == 1:
	return
    if score_with_all_seeds:
    	score_original(SCORING, score_commands, output_scores_file, log_file, job_file)
    else:
	score_xval(SCORING, score_xval_commands, output_scores_file, log_file, job_file) 
    return


def analyze(PPI, SCORING, output_scores_file, log_file, node_scores_file, association_scores_file_identifier_type, node_mapping_file, node_description_file, network_file_identifier_type, association_scores_validation_file, r_script_file, predictions_file, labels_file, tex_script_file, output_log_file, output_dir, title, specie, network_file, candidates_file, functional_enrichment):
    """
	Does cross validation and percentage analysis on the output files
    """
    if N_SEED == 1:
	return
    analyzed = analyze_xval(SCORING, r_script_file, output_scores_file, node_scores_file, predictions_file, labels_file, tex_script_file, output_log_file, output_dir, title, log_file, candidates_file, analysis_type = "auc") #"user") 
    if analyzed:
	analyze_xval_percentage(log_file, output_scores_file, node_scores_file, output_log_file)
	# No need to analyze original by default, do it explicitly if need be
	#analyze_original(PPI, output_scores_file, log_file, node_scores_file, association_scores_file_identifier_type, node_mapping_file, node_description_file, network_file_identifier_type, association_scores_validation_file, specie, network_file, functional_enrichment)
    return


def generate_score_xval_command(SCORING, score_xval_commands, k):
    #return score_xval_command % (node_scores_file, k, edge_scores_file, output_scores_file, k, N_REPETITION, N_ITERATION, score_log_file, k)
    return score_xval_commands[SCORING].substitute(fold = '%d' % k)


def score_xval(SCORING, score_xval_commands, output_scores_file, log_file, job_file):
    qname = None
    if use_cluster:
	if SCORING in ("nz", "ns"):
	    qname = "bigmem"
	elif SCORING in ("nr", "ff", "nd"):
	    qname = "bigmem" #"sbi-short"
	else: #elif SCORING in ("nd", "nw"):
	    qname = "sbi"
	qname = "sbi" 
    f = open(log_file, "a")
    for k in range(1, N_X_VAL+1):
	if not os.path.exists(output_scores_file + ".%d" % k): 
	    if only_print_command:
		print generate_score_xval_command(SCORING, score_xval_commands, k)
	    else:
		f.write("%s\n" % generate_score_xval_command(SCORING, score_xval_commands, k))
		if use_cluster:
		    #os.system( "qsub -cwd -o out.%d -e err.%d -l hostname=node52 -N %s -b y %s" % (k, k, SCORING, generate_score_xval_command(SCORING, score_xval_commands, k)) )
		    #if k % 2 != 0:
		    #	continue
		    os.system( "qsub -cwd -o out.%d -e err.%d -q %s -N %s -b y %s" % (k, k, qname, SCORING, generate_score_xval_command(SCORING, score_xval_commands, k)) )
		else:
		    os.system( generate_score_xval_command(SCORING, score_xval_commands, k) )
    f.close()
    return


def score_original(SCORING, score_commands, output_scores_file, log_file, job_file):
    qname = None
    if use_cluster:
	if SCORING in ("nz", "ns"):
	    qname = "bigmem"
	elif SCORING in ("nr", "ff", "nd"):
	    qname = "bigmem" #"sbi-short"
	else:
	    qname = "sbi"
	qname = "sbi" 
    if not os.path.exists(output_scores_file):
	f = open(log_file, "a")
	#print score_commands[SCORING]
	if only_print_command:
	    print score_commands[SCORING]
	else:
	    f.write("%s\n" % score_commands[SCORING])
	    if use_cluster:
		#os.system("qsub -cwd -o out -e err -l hostname=node34 -N %s -b y %s" % (SCORING, score_commands[SCORING]))
		os.system("qsub -cwd -o out -e err -q %s -N %s -b y %s" % (qname, SCORING, score_commands[SCORING]))
	    else:
		os.system(score_commands[SCORING])
	f.close()
    return


def analyze_xval(SCORING, r_script_file, output_scores_file, node_scores_file, predictions_file, labels_file, tex_script_file, output_log_file, output_dir, title, log_file, candidates_file, analysis_type="auc"):
    if analysis_type == "auc":
	#if os.path.exists(r_script_file):
	#    return False
	if not os.path.exists(r_script_file): 
	    list_node_scores_and_labels = []
	    previous_negative_sample_size = None
	    for k in range(1, N_X_VAL+1):
		node_validation_data, previous_negative_sample_size = analyze_results.get_validation_node_scores_and_labels(file_result = output_scores_file+".%d"%k, file_seed_test_scores = node_scores_file+".%d.test"%k, file_node_scores = node_scores_file, n_random_negative_folds = N_RANDOM_NEGATIVE_FOLDS, default_score = DEFAULT_NON_SEED_SCORE, replicable = REPLICABLE, candidates_file = candidates_file, previous_negative_sample_size=previous_negative_sample_size)
		list_node_scores_and_labels.append(node_validation_data)
	    analyze_results.create_ROCR_files(list_node_scores_and_labels, predictions_file, labels_file)
	    analyze_results.create_R_script(r_script_file, output_dir, title, only_auc=only_auc) # If only_auc is True only area under ROC curve is checked, other graphs are not drawn
	if not os.path.exists(os.path.dirname(r_script_file)+os.sep+"auc.txt"):
	    os.system("R CMD BATCH %s" % (r_script_file))
	    #os.system("convert %sperformance.eps %sperformance.jpg" % (output_dir, output_dir))
	    #analyze_results.create_tex_script(tex_script_file, output_dir, title)
	    analyze_results.record_performance_AUC_in_log_file(output_dir, output_log_file, title)
	    return True
	return False
    elif analysis_type == "user":
	#if not os.path.exists(candidates_file):
	#    print candidates_file
	#return
	# For user defined threshold performance analysis
	r_data_file = r_script_file.rsplit(".")[0] + ".dat"
	if os.path.exists(r_data_file):
	    f = open(r_data_file, "a")
	    #return False
	else:
	    f = open(r_data_file, "w")
	    f.write(" ppv sens\n")
	#thresholds = [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975 ]
	thresholds = THRESHOLDS[SCORING]
	counts_at_thresholds = [ [0.0] * 4 for tScore in thresholds] #nTP_sum, nFP_sum, nFN_sum, nTN_sum = 0.0, 0.0, 0.0, 0.0
	for k in range(1, N_X_VAL+1):
	    dictNodeResult, setNodeTest, non_seeds = analyze_results.get_values_from_files_for_performance_metric_counts(file_result = output_scores_file+".%d" % k, file_seed_test_scores = node_scores_file+".%d.test"%k, file_node_scores = node_scores_file, default_score = DEFAULT_NON_SEED_SCORE, candidates_file = candidates_file)
	    for i, tScore in enumerate(thresholds):
		nTP, nFP, nFN, nTN = analyze_results.calculate_performance_metric_counts(dictNodeResult, setNodeTest, non_seeds, score_threshold = tScore, n_random_negative_folds = N_RANDOM_NEGATIVE_FOLDS, replicable = REPLICABLE)
		for j, v in enumerate([nTP, nFP, nFN, nTN]):
		    counts_at_thresholds[i][j] += v
	for i, tScore in enumerate(thresholds):
	    nTP, nFP, nFN, nTN = counts_at_thresholds[i]
	    (acc, sens, spec, ppv) = analyze_results.calculatePerformance(nTP, nFP, nFN, nTN)
	    f.write("%f %s %s\n" % (tScore, ppv, sens))
	f.close()
	return False
    ## OLD - slightly slower
    elif analysis_type == "user2":
	# For user defined threshold performance analysis
	r_data_file = r_script_file.rsplit(".")[0] + ".dat2"
	if os.path.exists(r_data_file):
	    return False
	#thresholds = [ 0.02, 0.9 ] #[0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975 ]
	thresholds = THRESHOLDS[SCORING]
	f = open(r_data_file, "w")
	f.write(" ppv sens\n")
	#for tScore in [ i*0.01 for i in xrange(0, 100, 5)]:
	for tScore in thresholds:
	    #f.write("---- t: %s\n" % tScore)
	    nTP_sum, nFP_sum, nFN_sum, nTN_sum = 0.0, 0.0, 0.0, 0.0
	    for k in range(1, N_X_VAL+1):
		##print output_scores_file+".ns.%d"%k, node_scores_file+".%d.test"%k 
		nTP, nFP, nFN, nTN = analyze_results.calculate_performance_metric_counts_using_files(file_result = output_scores_file+".%d" % k, file_seed_test_scores = node_scores_file+".%d.test"%k, file_node_scores = node_scores_file, score_threshold = tScore, n_random_negative_folds = N_RANDOM_NEGATIVE_FOLDS, default_score = DEFAULT_NON_SEED_SCORE, replicable = REPLICABLE, candidates_file = candidates_file)
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
	    #f.write("Avg. A: %s S: %s P: %s\n" % (acc, sens, ppv))
	    f.write("%f %s %s\n" % (tScore, ppv, sens))
	    # Calculate 1-Specificity (1-TN/N[FP+TN] = FP/N) (aka FPR) 
	    fpr = 1-spec
	f.close()
	return False


def analyze_xval_percentage(log_file, output_scores_file, node_scores_file, output_log_file):
    f = open(log_file, "a")
    coverages = []
    f.write("\nXVAL SEED COVERAGE ANALYSIS\n")
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


def analyze_original(PPI, output_scores_file, log_file, node_scores_file, association_scores_file_identifier_type, node_mapping_file, node_description_file, network_file_identifier_type, association_scores_validation_file, specie, network_file, functional_enrichment = False):
    if not os.path.exists(output_scores_file):
	raise Exception("Output score file does not exist!")

    f = open(log_file, "a")

    #if association_scores_file_identifier_type is not None and not os.path.exists(output_scores_file+"."+association_scores_file_identifier_type):
    #	if PPI.startswith("biana"):
    #	    prepare_data.convert_file_using_new_id_mapping(output_scores_file, node_description_file, network_file_identifier_type, "geneid", output_scores_file+".geneid")
    #	    prepare_data.convert_file_using_new_id_mapping(output_scores_file+".geneid", gene_info_file, "geneid", association_scores_file_identifier_type, output_scores_file+"."+association_scores_file_identifier_type)
    #	elif PPI.startswith("ori") or PPI.startswith("david") or PPI.startswith("piana_joan"):
    #	    pass
    #	else:
    #	    prepare_data.convert_file_using_new_id_mapping(output_scores_file, node_description_file, network_file_identifier_type, association_scores_file_identifier_type, output_scores_file+"."+association_scores_file_identifier_type)

    f.write("\nSEED COVERAGE ANALYSIS\n")
    for percentage in (10, 25, 50):
	f.write("---- %s:\n" % percentage)
	#f.write("%s\n" % output_scores_file)
	coverage = analyze_results.calculate_seed_coverage_at_given_percentage(output_scores_file, node_scores_file, percentage, DEFAULT_NON_SEED_SCORE)
	f.write("%s\n" % str(coverage))
	#if association_scores_validation_file is not None:
	#    f.write("\nValidation seed coverage:\n")
	#    if PPI.startswith("ori") or PPI.startswith("david") or PPI.startswith("piana_joan"):
	#   pass
	#    else:
	#	f.write("%s\n" % str(analyze_results.calculate_seed_coverage_at_given_percentage(output_scores_file+"."+association_scores_file_identifier_type, association_scores_validation_file, percentage, DEFAULT_NON_SEED_SCORE)))

    if functional_enrichment and association_scores_file_identifier_type is not None and os.path.exists(node_mapping_file+"."+association_scores_file_identifier_type):
    	f.write("\nFUNCTIONAL ENRICHMENT ANALYSIS (OVER TOP %d%% SCORING NODES)\n" % DEFAULT_TOP_SCORING_PERCENTAGE)
    	analyze_results.check_functional_enrichment_at_given_percentage(output_scores_file, node_scores_file, node_mapping_file+"."+association_scores_file_identifier_type, DEFAULT_TOP_SCORING_PERCENTAGE, association_scores_file_identifier_type, f.write, DEFAULT_NON_SEED_SCORE, exclude_seeds = True, specie = specie, mode = "ordered")

	#for percentage in (10, 25, 50):
	#    f.write("---- %s:\n" % percentage)
	#    f.close()
	#    analyze_results.check_functional_enrichment_at_given_percentage(output_scores_file, node_scores_file, output_scores_file+"."+association_scores_file_identifier_type, percentage, association_scores_file_identifier_type, log_file, DEFAULT_NON_SEED_SCORE, exclude_seeds)
	#    f = open(log_file, "a")

    	f.write("\nFUNCTIONAL ENRICHMENT ANALYSIS OF MODULES (OVER TOP %d%% SCORING NODES)\n" % DEFAULT_TOP_SCORING_PERCENTAGE)
	analyze_results.check_functional_enrichment_of_high_scoring_modules(network_file, "greedy", output_scores_file, node_scores_file, node_mapping_file+"."+association_scores_file_identifier_type, DEFAULT_TOP_SCORING_PERCENTAGE, association_scores_file_identifier_type, f.write, DEFAULT_NON_SEED_SCORE, exclude_seeds = False, specie = specie, mode = "unordered")

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


def prepare_scoring_files(PPI, seed_scores_file, network_file_filtered, seed_to_score, node_scores_file, association_scores_file_identifier_type, node_mapping_file, node_description_file, network_file_identifier_type, edge_scores_file, interaction_relevance_file, interaction_relevance_file2, edge_scores_as_node_scores_file, sampled_file_prefix):
    #print "s:%d x:%d" % (N_SEED, N_X_VAL)
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
	prepare_data.generate_cross_validation_node_score_files(nodes = nodes, seed_to_score = seed_to_score, node_scores_file = node_scores_file, xval = N_X_VAL, default_score = DEFAULT_NON_SEED_SCORE, replicable = REPLICABLE)

    # Create node id to genesymbol mapping
    if association_scores_file_identifier_type is not None and not os.path.exists(node_mapping_file+"."+association_scores_file_identifier_type): 
	if node_description_file is not None:
	    print "Creating node mapping file", node_mapping_file+"."+association_scores_file_identifier_type
	    if PPI.startswith("biana"):
		if association_scores_file_identifier_type == "genesymbol":
		    #prepare_data.convert_file_using_new_id_mapping(node_scores_file, node_description_file, network_file_identifier_type, "geneid", node_mapping_file+".geneid", id_to_id_mapping=True)
		    #prepare_data.convert_file_using_new_id_mapping(node_mapping_file+".geneid", gene_info_file, "geneid", association_scores_file_identifier_type, node_mapping_file+"."+association_scores_file_identifier_type, id_to_id_mapping=True)
		    prepare_data.convert_file_using_new_id_mapping(node_scores_file, node_description_file, network_file_identifier_type, association_scores_file_identifier_type, node_mapping_file+"."+association_scores_file_identifier_type, id_to_id_mapping=True, intermediate_mapping_file=gene_info_file, intermediate_mapping_id_type="geneid")
		else:
		    prepare_data.convert_file_using_new_id_mapping(node_scores_file, node_description_file, network_file_identifier_type, association_scores_file_identifier_type, node_mapping_file+"."+association_scores_file_identifier_type, id_to_id_mapping=True)
	    else:
		prepare_data.convert_file_using_new_id_mapping(node_scores_file, node_description_file, network_file_identifier_type, association_scores_file_identifier_type, node_mapping_file+"."+association_scores_file_identifier_type, id_to_id_mapping=True)

    # Create edge scores (original + xval) files as well as node scores as edge scores files
    edges = None
    edge_to_score = None
    if not os.path.exists(edge_scores_file): 
	print "Creating edge score files", edge_scores_file
	edges = prepare_data.get_edges_in_network(network_file = network_file_filtered)
	if all([interaction_relevance_file, interaction_relevance_file2]):
	    edge_to_score = prepare_data.get_edge_to_score_from_sif_attribute_file(interaction_relevance_file)
	    edge_to_score = dict([(e, (sum([float(i) for i in v])/len(v))/1000 + 1) for e,v in edge_to_score.iteritems()])
	    edge_to_score2 = prepare_data.get_edge_to_score_from_sif_attribute_file(interaction_relevance_file2)
	    edge_to_score2 = dict([(e, (sum([float(i) for i in v])/len(v))/1000) for e,v in edge_to_score2.iteritems()])
	    edge_to_score = dict([(e, v+edge_to_score2[e]) for e,v in edge_to_score.iteritems()])
	elif interaction_relevance_file is not None:
	    edge_to_score = dict([(e, 1) for e in edges])
	    edge_to_score_string = prepare_data.get_edge_to_score_from_sif_attribute_file(interaction_relevance_file)
	    # Normalizing string score (0.001 is added to avoid 0 score edges)
	    #edge_to_score = dict([(e, (sum([float(i) for i in v])/len(v))/1000 + 1) for e,v in edge_to_score_string.iteritems()])
	    for e,v in edge_to_score_string.iteritems():
		score = 1 + (sum([float(i) for i in v]) / len(v)) # / 1000 #! this is for STRING only
		a,b = e
		if e in edge_to_score:
		    edge_to_score[e] = score
		elif (b,a) in edge_to_score:
		    edge_to_score[(b,a)] = score
	else:
	    edge_to_score = dict([(e, 1) for e in edges])
	prepare_data.create_edge_scores_file(network_file = network_file_filtered, edge_scores_file = edge_scores_file, edge_to_score = edge_to_score, default_score=DEFAULT_NON_SEED_SCORE)
    if not os.path.exists(edge_scores_as_node_scores_file): 
	if edges is None:
	    edges = prepare_data.get_edges_in_network(network_file = network_file_filtered, data=True)
	    if edge_to_score is None:
		edge_to_score = dict([ ((u,v), score) for u,v,score in edges ])
	    edges = [ (u,v) for u,v,score in edges ]
	seed_to_score = prepare_data.get_node_to_score_from_node_scores_file(seed_scores_file)
	prepare_data.create_edge_scores_as_node_scores_file(edges = edges, node_to_score = seed_to_score, edge_to_score = edge_to_score, edge_scores_file = edge_scores_as_node_scores_file, ignored_nodes = None, default_score = DEFAULT_NON_SEED_SCORE)
	prepare_data.generate_cross_validation_edge_score_as_node_score_files(edges = edges, seed_to_score = seed_to_score, edge_to_score = edge_to_score, edge_scores_file = edge_scores_as_node_scores_file, xval = N_X_VAL, default_score = DEFAULT_NON_SEED_SCORE, replicable = REPLICABLE)
    # Create random network files
    if not os.path.exists(sampled_file_prefix + ".sif.1"): 
	print "Creating sampled networks"
	prepare_data.sample_network_preserving_topology(edge_scores_file, N_SAMPLE_GRAPH, sampled_file_prefix + ".sif.")
    return



if __name__ == "__main__":
    multiple = False 
    if multiple == False:
	main(ppis, phenotypes, scoring_parameters)
    else:
	for i_parameter in range(10,110,10): 
	    ufi = user_friendly_id + "_permuted_p%d-omim" % i_parameter
	    #ufi = user_friendly_id + "_pruned_p%d-omim" % i_parameter
	    #ufi = user_friendly_id + "-omim_perturbed_p%d" % i_parameter
	    #ppis = ["biana_no_tap_no_reliability_pruned_non_seed_interactions_p%s_%s" % (p, i) for p in xrange(i_parameter,i_parameter+10,10) for i in xrange(1,101) ] 
	    ppis = ["biana_no_tap_no_reliability_permuted_p%s_%s" % (p, i) for p in xrange(i_parameter,i_parameter+10,10) for i in xrange(1,101) ] 
	    #ppis = ["goh_permuted_p%s_%s" % (p, i) for p in xrange(i_parameter,i_parameter+10,10) for i in xrange(1,101) ] 
	    #ppis = ["goh_pruned_p%s_%s" % (p, i) for p in xrange(i_parameter,i_parameter+10,10) for i in xrange(1,101) ] 
	    #phenotypes = [ "perturbed_%s_p%i_%i" % (d, p, i) for d in omim_phenotypes for p in xrange(i_parameter,i_parameter+10,10) for i in xrange(1,101) ]  
	    #phenotypes = [ "perturbed_%s_p%i_%i" % (d, p, i) for d in ["omim_systemic_lupus_erythematosus"] for p in xrange(i_parameter,i_parameter+10,10) for i in xrange(11,101) ]  
	    main(ppis, phenotypes, scoring_parameters, ufi)

