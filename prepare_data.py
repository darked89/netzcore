#########################################################################
# Methods for preparing input files for scoring methods
#
# eg 23/06/2009
#########################################################################

from biana.utilities import biana_output_converter
from biana.utilities import graph_utilities as network_utilities

def generate_cross_validation_node_files(g, node_to_score, seeds, node_file_netzcore_prefix, edge_file_netzcore_relevance_prefix, arff_file_prefix, edge_file_netzcore):
    # seeds = [ v for v, s in node_to_score.iteritems() if s != 0 ] # was assuming that non-seeds should have zero score
    for k, training, test in k_fold_cross_validation(seeds, 10, randomize = True, useSeed = True):
	assign_node_scores(g = g, node_to_score = node_to_score, out_file = node_file_netzcore_prefix+"_%i.txt"%k, as_edge_relevance_score = False, seeds_ignored=test)
	assign_node_scores(g = g, node_to_score = node_to_score, out_file = edge_file_netzcore_relevance_prefix+"_%i.txt"%k, as_edge_relevance_score = True, seeds_ignored=test)
	create_edge_file_from_weight_and_score_files(edge_file_netzcore_relevance_prefix+"_%i.txt"%k, edge_file_netzcore)
	network_utilities.create_arff_file_with_network_metrics(g, node_to_score, training, arff_file_prefix+"_%i.arff"%k)
    return

def create_edge_file_from_weight_and_score_files(edge_file_weights, edge_file_scores, out_file):
    setNode, setEdge, dictNode, dictEdgeWeight = network_utilities.get_nodes_and_edges_from_sif_file(file_name = edge_file_weights[:-3]+"sif", store_edge_type = True)
    setNode, setEdge, dictNode, dictEdge = network_utilities.get_nodes_and_edges_from_sif_file(file_name = edge_file_scores[:-3]+"sif", store_edge_type = True)
    f = open(out_file, "w")
    for e, s in dictEdge.iteritems():
	u,v=e
	s = float(s)
	w = dictEdgeWeight[e]
	s = s*100 + 1
	#if s == 0:
	#    s = 0.1
	w /= s
	#print w
	f.write("%s %s %s\n" % (u, w, v))
    return 

def k_fold_cross_validation(X, K, randomize = False, useSeed = False):
	"""
	By John Reid (code.activestate.com)
	Generates K (training, validation) pairs from the items in X.

	Each pair is a partition of X, where validation is an iterable
	of length len(X)/K. So each training iterable is of length (K-1)*len(X)/K.

	If randomise is true, a copy of X is shuffled before partitioning,
	otherwise its order is preserved in training and validation.
	"""
	#if randomize: from random import shuffle; X=list(X); shuffle(X)
	if randomize: 
	    from random import shuffle, seed
	    X=list(X)
	    if useSeed: 
		random.seed(123)
	    shuffle(X)
	for k in xrange(K):
	    training = [x for i, x in enumerate(X) if i % K != k]
	    validation = [x for i, x in enumerate(X) if i % K == k]
	    yield k, training, validation


def assign_edge_reliability_scores(g, network_file_method_attribute, network_file_source_attribute, network_file_pubmed_attribute, out_file):
    # linear combination of 
    # pubmed/5
    # db/3
    # method/2
    # jaccard/0.4
    # cc1*cc2/0.4
    edge_to_methods = network_utilities.get_edge_values_from_sif_attribute_file(file_name = network_file_method_attribute, store_edge_type = False)
    edge_to_sources = network_utilities.get_edge_values_from_sif_attribute_file(file_name = network_file_source_attribute, store_edge_type = False)
    edge_to_pubmeds = network_utilities.get_edge_values_from_sif_attribute_file(file_name = network_file_pubmed_attribute, store_edge_type = False)
    edge_to_jaccard = network_utilities.get_jaccard_index_map(g)
    node_to_ccoef = network_utilities.get_clustering_coefficient_map(g)
    f = open(out_file, 'w')
    for u,v in g.edges_iter():
	score = 0.0
	if edge_to_methods.has_key((u,v)):
	    score += len(edge_to_methods[(u,v)])/2.0
	elif edge_to_methods.has_key((v,u)):
	    score += len(edge_to_methods[(v,u)])/2.0
	if edge_to_sources.has_key((u,v)):
	    score += len(edge_to_sources[(u,v)])/3.0
	elif edge_to_sources.has_key((v,u)):
	    score += len(edge_to_sources[(v,u)])/3.0
	if edge_to_pubmeds.has_key((u,v)):
	    score += len(edge_to_pubmeds[(u,v)])/5.0
	elif edge_to_pubmeds.has_key((v,u)):
	    score += len(edge_to_pubmeds[(v,u)])/5.0
	score += node_to_ccoef[u]*node_to_ccoef[v]/0.4
	score += edge_to_jaccard[(u,v)]/0.4
	f.write("%s\t%s\t%f\n" % (u, v, score))
    f.close()
    return

def assign_node_scores(g, node_to_score, out_file, as_edge_relevance_score=False, seeds_ignored = None):
    f = open(out_file, 'w')
    if as_edge_relevance_score:
	for u,v in g.edges_iter():
	    if seeds_ignored is not None and u in seeds_ignored:
		score_u = 0
	    else:
		score_u = node_to_score[u]
	    if seeds_ignored is not None and v in seeds_ignored:
		score_v = 0
	    else:
		score_v = node_to_score[v]
	    f.write("%s\t%s\t%f\n" % (u, v, (score_u + score_v) / 2))
    else:
	for v in g.nodes_iter():
	    if seeds_ignored is not None and v in seeds_ignored:
		score_v = 0
	    else:
		score_v = node_to_score[v]
	    f.write("%s\t%f\n" % (v, score_v))
	f.close()
    return

def check_seed_coverage_and_get_seed_nodes(g, node_file, node_file_score):
    setNode, setDummy, dictNode, dictEdge = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_file_score, store_edge_type = False)
    node_to_genes, gene_to_nodes = biana_output_converter.get_user_entity_gene_mapping(node_file)
    covered_genes = set()
    seeds = set()
    for v in g.nodes_iter():
	gene_with_score = setNode & node_to_genes[v]
	covered_genes |= gene_with_score
	if len(gene_with_score) > 0:
	    seeds.add(v)
	if len(gene_with_score) > 1:
	    print "More than one gene??", gene_with_score
    print "Covered genes:", len(covered_genes)
    print "Covered gene products (seed nodes):", len(seeds)
    return seeds

def map_scores_to_biana_nodes(node_file, node_file_score):
    setNode, setDummy, dictNode, dictEdge = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_file_score)
    node_to_genes, gene_to_nodes = biana_output_converter.get_user_entity_gene_mapping(node_file)
    node_to_score = {}
    covered_genes = set()
    seeds = set()
    for id, genes in node_to_genes.iteritems():
	gene_with_score = setNode & node_to_genes[id]
	covered_genes |= gene_with_score
	if len(gene_with_score) == 0:
	    node_to_score[id] = 0
	else:
	    #if len(gene_with_score) > 1:
	    #	print "More than one gene??", gene_with_score
	    seeds.add(id)
	    i = 0
	    score = 0.0
	    for gene in covered_genes:
		i+=1
		score += float(dictNode[gene])
	    node_to_score[id] = score/i
	    if score/i <= 0:
		print "non-positive seed score", gene, score/i
    return node_to_score, seeds 

def create_method_filtered_network_files(network_file_method_attribute):
    biana_output_converter.filter_network_by_interaction_type(network_attribute_file_name = network_file_method_attribute, network_out_file_name="/home/emre/arastirma/data/human_interactome_biana/human_y2h.sif", interaction_type="y2h")
    biana_output_converter.filter_network_by_interaction_type(network_attribute_file_name = network_file_method_attribute, network_out_file_name="/home/emre/arastirma/data/human_interactome_biana/human_tap.sif", interaction_type="tap")
    biana_output_converter.filter_network_by_interaction_type(network_attribute_file_name = network_file_method_attribute, network_out_file_name="/home/emre/arastirma/data/human_interactome_biana/human_no_tap.sif", interaction_type="tap", reverse_selection=True)
    return

def extract_test_data_from_arff_files(arff_file_all, arff_file_filtered, arff_file_training, arff_file_test):
    f_all = open(arff_file_all)
    f_training = open(arff_file_training)
    f_filtered = open(arff_file_filtered)
    f_out = open(arff_file_test, 'w')
    line_all = f_all.readline()
    line_training = f_training.readline()
    line_filtered = f_filtered.readline()
    while line_all.startswith("@"):
	f_out.write(line_filtered)
	line_all = f_all.readline()
	line_training = f_training.readline()
	line_filtered = f_filtered.readline()
    while line_all:
	words_all = line_all.split(",")
	words_training = line_training.split(",")
	if(words_all[1] != words_training[1]):
	    f_out.write(line_filtered)
	line_all = f_all.readline()
	line_training = f_training.readline()
	line_filtered = f_filtered.readline()
    f_out.close()
    return

def generate_cross_validation_test_test_arff_files(k, arff_file_all, arff_file_filtered, arff_file_training_prefix):
    for i in xrange(k):
	extract_test_data_from_arff_files(arff_file_all, arff_file_filtered, arff_file_training_prefix+"_%i.arff"%i, arff_file_training_prefix+"_test_%i.arff"%i)
    return

def sample_network_preserving_topology(network_sif_file, n_sample, output_prefix):
    g = network_utilities.create_network_from_sif_file(network_file_in_sif = network_sif_file, weighted = True, delim = " ")
    for i in xrange(1,n_sample+1):
	g_sampled = network_utilities.randomize_graph(graph=g, randomization_type="preserve_topology_and_node_degree")
	#network_utilities.output_network_in_sif(g_sampled, output_prefix+"%s"%i)
	network_utilities.output_network_in_txt(g_sampled, output_prefix+"%s"%i)
    return

if __name__ == "__main__":
    pass

