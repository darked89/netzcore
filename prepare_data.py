#########################################################################
# Methods for preparing input files for scoring methods
#
# eg 23/06/2009
#########################################################################

from biana.utilities import biana_output_converter
from biana.utilities import TsvReader
from biana.utilities import graph_utilities as network_utilities


def get_nodes_in_network(network_file):
    g = network_utilities.create_network_from_sif_file(network_file)
    return g.nodes()


def get_edges_in_network(network_file):
    g = network_utilities.create_network_from_sif_file(network_file)
    return g.edges()

#def get_degrees(network_file):
#    g = network_utilities.create_network_from_sif_file(network_file)
#    return g.degree()

def get_node_to_score_from_node_scores_file(node_scores_file):
    nodes, set_dummy, node_to_score, dict_dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)
    return node_to_score

def get_nodes_from_nodes_file(node_scores_file):
    nodes, set_dummy, dict_dummy, dict_dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)
    return nodes

def get_edge_to_score_from_sif_attribute_file(interaction_relevance_file):
    return network_utilities.get_edge_values_from_sif_attribute_file(file_name = interaction_relevance_file, store_edge_type = False)

def generate_cross_validation_edge_score_as_node_score_files(edges, seed_to_score, edge_scores_file, xval = 5, default_score = 0):
    seeds = seed_to_score.keys()
    for k, training, test in k_fold_cross_validation(seeds, xval, randomize = True, replicable = True):
	create_edge_scores_as_node_scores_file(edges = edges, node_to_score = seed_to_score, edge_scores_file = edge_scores_file+".%i"%k, ignored_nodes = test, default_score = default_score)
    return


def generate_cross_validation_node_score_files(nodes, seed_to_score, node_scores_file, xval = 5, default_score = 0):
    seeds = seed_to_score.keys()
    for k, training, test in k_fold_cross_validation(seeds, xval, randomize = True, replicable = True):
	create_node_scores_file(nodes = nodes, node_to_score = seed_to_score, node_scores_file = node_scores_file+".%i"%k, ignored_nodes = test, default_score = default_score )
	create_node_scores_file(nodes = test, node_to_score = seed_to_score, node_scores_file = node_scores_file+".%i.test"%k, ignored_nodes = None, default_score = default_score )
    return


def k_fold_cross_validation(X, K, randomize = False, replicable = False):
    """
    By John Reid (code.activestate.com)
    Generates K (training, validation) pairs from the items in X.

    Each pair is a partition of X, where validation is an iterable
    of length len(X)/K. So each training iterable is of length (K-1)*len(X)/K.

    If randomise is true, a copy of X is shuffled before partitioning,
    otherwise its order is preserved in training and validation.

    If replicable is true, an internal number (123) is used to create the same random splits at each call
    """
    #if randomize: from random import shuffle; X=list(X); shuffle(X)
    if randomize: 
	from random import shuffle, seed
	X=list(X)
	if replicable: 
	    seed(123)
	shuffle(X)
    for k in xrange(K):
	training = [x for i, x in enumerate(X) if i % K != k]
	validation = [x for i, x in enumerate(X) if i % K == k]
	yield k+1, training, validation
    return


def sample_network_preserving_topology(network_sif_file, n_sample, output_prefix):
    g = network_utilities.create_network_from_sif_file(network_file_in_sif = network_sif_file, use_edge_data = True)#, delim = " ")
    for i in xrange(1,n_sample+1):
	g_sampled = network_utilities.randomize_graph(graph=g, randomization_type="preserve_topology_and_node_degree")
	network_utilities.output_network_in_sif(g_sampled, output_prefix+"%s"%i)
    return


def create_edge_scores_file(network_file, edge_scores_file, edge_to_score = None, default_score=1):
    g = network_utilities.create_network_from_sif_file(network_file)
    f = open(edge_scores_file, 'w')
    for u,v in g.edges_iter():
	score = default_score
	if edge_to_score.has_key((u,v)):
	    score = edge_to_score[(u,v)]
	elif edge_to_score.has_key((v,u)):
	    score = edge_to_score[(v,u)]
	#else:
	#    print (u,v), 
	f.write("%s %f %s\n" % (u, score, v))
    f.close()
    return


def old_create_edge_scores_as_node_scores_file(edges, node_to_score, edge_scores_file, ignored_nodes = None, default_score = 0):
    """
	Creates edge score file from node association scores, intended comparing netshort with other algorithms without using other edge reliability/relevance score
    """
    g = network_utilities.create_network_from_sif_file(network_file)
    setNode, setDummy, dictNode, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)
    f = open(edge_scores_file, 'w')
    for u,v in g.edges_iter():
	if ignored_nodes is not None and u in ignored_nodes:
	    score_u = default_score
	else:
	    score_u = dictNode[u]
	if ignored_nodes is not None and v in ignored_nodes:
	    score_v = default_score
	else:
	    score_v = dictNode[v]
	f.write("%s %f %s\n" % (u, (score_u + score_v) / 2, v))
    f.close()
    return


def create_edge_scores_as_node_scores_file(edges, node_to_score, edge_scores_file, ignored_nodes = None, default_score = 0):
    """
	Creates edge score file from node association scores, intended comparing netshort with other algorithms without using other edge reliability/relevance score
    """
    f = open(edge_scores_file, 'w')
    for u,v in edges:
	if ignored_nodes is not None and u in ignored_nodes:
	    score_u = default_score
	else:
	    if node_to_score.has_key(u):
		score_u = node_to_score[u]
	    else:
		score_u = default_score
	if ignored_nodes is not None and v in ignored_nodes:
	    score_v = default_score
	else:
	    if node_to_score.has_key(v):
		score_v = node_to_score[v]
	    else:
		score_v = default_score
	f.write("%s %f %s\n" % (u, (score_u + score_v) / 2, v))
    f.close()
    return

def create_node_scores_file(nodes, node_to_score, node_scores_file, ignored_nodes = None, default_score = 0):
    """
	Creates node association/relevance score file for the nodes in the network 
	ignored_nodes: nodes in this set are ignored (assigned default_score)
    """
    f = open(node_scores_file, 'w')
    for v in nodes:
	if ignored_nodes is not None and v in ignored_nodes:
	    score_v = default_score
	else:
	    if node_to_score.has_key(v):
		score_v = node_to_score[v]
	    else:
		score_v = default_score
	f.write("%s %f\n" % (v, score_v))
    f.close()
    return 

def convert_file_using_new_id_mapping(file_to_be_converted, node_description_file, from_id_type, to_id_type, new_file):
    """
	Maps nodes given as from_id_type to their correspondants in to_id_type using node_description_file
	Can convert node / network file in sif format (node file with data, network file with data in the middle)
    """
    nodes, edges, node_to_data, edge_to_data = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_to_be_converted, store_edge_type = True)

    #node_id_to_new_ids, dummy = biana_output_converter.get_attribute_to_attribute_mapping(node_description_file, from_id_type, to_id_type, keys_to_include=nodes, include_inverse_mapping = False)
    reader = TsvReader.TsvReader(node_description_file, inner_delim = ",")
    columns, node_id_to_new_ids = reader.read(fields_to_include = [from_id_type, to_id_type], keys_to_include=nodes, merge_inner_values = True)

    f = open(new_file, 'w')
    if edges is None:
	for v in nodes:
	    if node_to_data.has_key(v):
		if node_id_to_new_ids.has_key(v):
		    vals = reduce(lambda x,y: x+y, node_id_to_new_ids[v])
		    for id in vals:
			f.write("%s %s\n" % (id, node_to_data[v]))
	    else:
		if node_id_to_new_ids.has_key(v):
		    vals = reduce(lambda x,y: x+y, node_id_to_new_ids[v])
		    for id in vals:
			f.write("%s\n", id)
    else:
	for e in edges:
	    u,v = e
	    if edge_to_data.has_key(e):
		if node_id_to_new_ids.has_key(u):
		    vals = reduce(lambda x,y: x+y, node_id_to_new_ids[u])
		    for id in vals:
			if node_id_to_new_ids.has_key(v):
			    vals2 = reduce(lambda x,y: x+y, node_id_to_new_ids[v])
			    for id2 in vals2:
				f.write("%s %s %s\n" % (id, edge_to_data[e], id2))
    f.close()
    return

def get_node_association_score_mapping(network_file, network_file_identifier_type, node_description_file, association_scores_file, association_scores_file_identifier_type, log_file = None, default_seed_score=1.0):
    """
	Maps genes and their scores to nodes in the network using given association_scores_file, correspondance identifiers
    """
    g = network_utilities.create_network_from_sif_file(network_file)
    nodes = g.nodes()
    setNode, setDummy, dictNode, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = association_scores_file, store_edge_type = False)
    if dictNode is None:
	dictNode = dict([ (v, default_seed_score) for v in setNode ])
    node_to_genes, gene_to_nodes = biana_output_converter.get_attribute_to_attribute_mapping(node_description_file, network_file_identifier_type, association_scores_file_identifier_type, keys_to_include=set(nodes))
    covered_genes = set()
    seeds = set()
    node_to_score = {}
    if log_file is not None:
	log_fd = open(log_file, "a")
    else:
	log_fd = None
    for v in nodes: 
	gene_with_score = setNode & node_to_genes[v]
	covered_genes |= gene_with_score

	if len(gene_with_score) > 0:
	    seeds.add(v)
	    if len(gene_with_score) > 1:
		#print "More than one gene:", gene_with_score, "for", v
		if log_fd is not None:
		    log_fd.write("More than one gene: %s for %s\n" % (gene_with_score, v))
	    i=0
	    score = 0
	    for gene in covered_genes:
		i+=1
		score += float(dictNode[gene])
	    score /= i
	    if score <= 0:
		#print "non-positive seed score", v, score, "genes:", node_to_genes[v]
		if log_fd is not None:
		    log_fd.write("non-positive seed score %s %s genes: %s\n" % (v, score, node_to_genes[v]))
	    node_to_score[v] = score
	#else:
	#    score = default_score
	#node_to_score[v] = score

    #print "Covered genes (seed genes):", len(covered_genes), "among", len(setNode)
    #print "Covered gene products (seed nodes):", len(seeds), "among", g.number_of_nodes()
    if log_fd is not None:
	log_fd.write("Covered genes (seed genes): %s among %s\n" % (len(covered_genes), len(setNode)))
	log_fd.write("Covered gene products (seed nodes): %s among %s\n" % (len(seeds), g.number_of_nodes()))
	log_fd.close()
    return node_to_score


def analyze_network(network_file, out_file = None):
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data = False)
    network_utilities.analyze_network(g, out_file = out_file)
    return


def create_R_analyze_network_script(network_file, seeds=None, out_path="./", title = ""):
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data = False)
    network_utilities.create_R_analyze_network_script(g, seeds, out_path, title)
    network_utilities.create_R_analyze_network_script(g, seeds, out_path, title, scale_by_log=True)
    return


def create_ARFF_network_metrics_file(network_file, node_to_score, seeds, arff_file_name):
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data = False)
    network_utilities.create_ARFF_network_metrics_file(g, node_to_score, seeds, arff_file_name)
    return


def create_degree_filtered_network_file(network_file, network_file_filtered, degree):
    """
	Creates a network file removing nodes that has connections more than given degree
    """
    # Load network
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data = True)
    #network_utilities.analyze_network(g)
    g = network_utilities.filter_network(g = g, degree_threshold = degree, largest_connected_component=True) 
    # Get degrees of highly connected nodes
    #network_utilities.analyze_network(g)
    network_utilities.output_network_in_sif(g, network_file_filtered)
    return


def create_method_filtered_network_files(network_file_prefix):
    """
	Creates 3 new sif files from the given network file where interactions coming from non-y2h, non-tap and tap are filtered.
    """
    network_file_method_attribute = network_file_prefix + "_method_id.eda"
    biana_output_converter.filter_network_by_interaction_type(network_attribute_file_name = network_file_method_attribute, network_out_file_name = network_file_prefix + "_y2h.sif", interaction_type="y2h")
    biana_output_converter.filter_network_by_interaction_type(network_attribute_file_name = network_file_method_attribute, network_out_file_name = network_file_prefix + "_tap.sif", interaction_type="tap")
    biana_output_converter.filter_network_by_interaction_type(network_attribute_file_name = network_file_method_attribute, network_out_file_name = network_file_prefix + "_no_tap.sif", interaction_type="tap", reverse_selection=True)
    return


def create_edge_reliability_filtered_network_file(network_file, network_file_prefix, out_file):
    # linear combination of 
    # pubmed/5
    # db/3
    # method/2
    # jaccard/0.4
    # cc1*cc2/0.4
    edge_to_methods = network_utilities.get_edge_values_from_sif_attribute_file(file_name = network_file_prefix + "_method_id.eda", store_edge_type = False)
    edge_to_sources = network_utilities.get_edge_values_from_sif_attribute_file(file_name = network_file_prefix + "_source.eda", store_edge_type = False)
    edge_to_pubmeds = network_utilities.get_edge_values_from_sif_attribute_file(file_name = network_file_prefix + "_pubmed.eda", store_edge_type = False)
    #edge_to_jaccard = network_utilities.get_jaccard_index_map(g)
    #node_to_ccoef = network_utilities.get_clustering_coefficient_map(g)
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data = True)
    f = open(out_file, 'w')
    for u,v in g.edges_iter():
	score = 0.0
	if edge_to_methods.has_key((u,v)):
	    score += len(edge_to_methods[(u,v)])/2.0
	    score_method = len(edge_to_methods[(u,v)])
	elif edge_to_methods.has_key((v,u)):
	    score += len(edge_to_methods[(v,u)])/2.0
	    score_method = len(edge_to_methods[(v,u)])
	if edge_to_sources.has_key((u,v)):
	    score += len(edge_to_sources[(u,v)])/3.0
	    score_source = len(edge_to_sources[(u,v)])
	elif edge_to_sources.has_key((v,u)):
	    score += len(edge_to_sources[(v,u)])/3.0
	    score_source = len(edge_to_sources[(v,u)])
	if edge_to_pubmeds.has_key((u,v)):
	    score += len(edge_to_pubmeds[(u,v)])/5.0
	    score_pubmed = len(edge_to_pubmeds[(u,v)])
	elif edge_to_pubmeds.has_key((v,u)):
	    score += len(edge_to_pubmeds[(v,u)])/5.0
	    score_pubmed = len(edge_to_pubmeds[(v,u)])
	#score += node_to_ccoef[u]*node_to_ccoef[v]/0.4
	#score += edge_to_jaccard[(u,v)]/0.4
	#f.write("%s\t%s\t%f\n" % (u, v, score))
	if score_source + score_pubmed > 2: # and score_pubmed > 2:
	    f.write("%s\t%s\t%s\n" % (u, g.get_edge(u,v), v))
    f.close()
    return


def create_human_interactome_using_biana(node_file_prefix="human_nodes", network_files_prefix="human_network", network_type="experimental", load_from_saved_session = False): 
    """
	Creates human ppi network files using BIANA
	Node file: node_file_prefix.tsv (as tab seperated values)
	Edge files: network_files_prefix.sif (network in sif format), network_file_prefix_attribute.eda (Cytoscape edge attribute files for each attribute)
	Network type: "experimental" (type interaction, all ppi db but string) | "functional" (type functional_association, only string) | "all" (experimental + functional)
    """
    from biana.biana_commands import available_sessions, create_new_session, save_session, load_session
    from biana.BianaObjects import BianaSessionManager
    
    DB_NAME = "test_biana" 
    DB_HOST = "127.0.0.1" 
    DB_USER = "biana_user" 
    DB_PASS = "biana_pass" 
    UNIFICATION_PROTOCOL = "(p2)(noself)(noprevious)uniprot_seqtax_geneid_scoppdb" 

    identifier_description_list = [("taxid", "9606")]
    level = 0
    relation_type_list = []
    relation_attribute_restriction_list = []
    node_attributes=["uniprotaccession", "geneid", "ensembl", "genesymbol", "go", "hgnc", "disease", "proteinsequence"]
    relation_attributes=["pubmed", "method_id"]
    if network_type == "functional":
	relation_type_list = ["functional_association"] 
	relation_attributes.extend(["stringscore", "stringscore_neighborhood", "stringscore_fusion", "stringscore_cooccurence", "stringscore_coexpression", "stringscore_experimental", "stringscore_db", "stringscore_textmining"])
	# Uncommet below for only relations from STRING with score > 700
	#relation_attribute_restriction_list=[("stringscore",">700")]
    elif network_type == "experimental":
	relation_type_list = ["interaction"] #,"reaction"] # ["interaction","reaction","complex","pathway"]
    elif network_type == "all":
	relation_type_list ["interaction", "functional_association"]
    else:
	raise ValueError("Unrecognized network_type!")

    if load_from_saved_session:
    	load_session(network_files_prefix + "_session.dat")
	objSession = available_sessions["biana_session"]
    else:
	# Create a new BIANA session
	create_new_session(sessionID="biana_session",dbname=DB_NAME, dbhost=DB_HOST, dbuser=DB_USER, dbpassword=DB_PASS, unification_protocol=UNIFICATION_PROTOCOL)
	objSession = available_sessions["biana_session"]

	# Create a new User Entity Set (group of biomolecules) in this session
	uESet1 = objSession.create_new_user_entity_set( identifier_description_list=identifier_description_list, attribute_restriction_list=[], id_type="embedded", new_user_entity_set_id="User_Entity_Set_1")
	print "# of nodes: ", objSession.get_user_entity_set("User_Entity_Set_1").getSize()

	# Fetch relations of biomolecules in the set 
	objSession.create_network( user_entity_set_id = "User_Entity_Set_1" , level = level, relation_type_list=relation_type_list, relation_attribute_restriction_list = relation_attribute_restriction_list, include_relations_last_level=True, use_self_relations=True)
	print "# of edges: ", objSession.get_user_entity_set("User_Entity_Set_1").getNumberEdges()

	# Save the session
	save_session("biana_session", network_files_prefix + "_session.dat")

    # Output set details - all nodes
    objSession.output_user_entity_set_details(user_entity_set_id = "User_Entity_Set_1", out_method = open(node_file_prefix + ".tsv", "w").write, attributes=node_attributes, include_level_info=False, include_degree_info=False, level=None, only_selected=False, output_format="tabulated", include_tags_info = False, include_tags_linkage_degree_info=[], substitute_node_attribute_if_not_exists=False, output_1_value_per_attribute=False, include_command_in_rows=False, output_only_native_values=False)

    objSession.output_user_entity_set_details(user_entity_set_id = "User_Entity_Set_1", out_method = open(node_file_prefix + "_only_native.tsv", "w").write, attributes=["uniprotaccession"], include_level_info=False, include_degree_info=False, level=None, only_selected=False, output_format="tabulated", include_tags_info = False, include_tags_linkage_degree_info=[], substitute_node_attribute_if_not_exists=False, output_1_value_per_attribute=False, include_command_in_rows=False, output_only_native_values=True)

    # Export all the information of this set and its network to a file in a tab separated format
    objSession.output_user_entity_set_network_in_sif_format(user_entity_set_id = "User_Entity_Set_1", output_path = "./", output_prefix = network_files_prefix, node_attributes = [], participant_attributes = [], relation_attributes = relation_attributes, output_1_value_per_attribute = False, include_tags = False) 

    ## Export sequences of entries in this set as a FASTA file
    #objSession.output_user_entity_set_sequences_in_fasta(user_entity_set_id="User_Entity_Set_1", only_selected = False, out_method=open("human_fasta.txt",'w').write, type="proteinsequence") #, attributes=["uniprotaccession"])

    ## Export all the information of this set and its network to a file in a tab separated format
    #available_sessions["biana_session"].output_user_entity_set_network(user_entity_set_id = "User_Entity_Set_1", out_method = open(network_files_prefix + ".tsv", "w").write, node_attributes = [], participant_attributes = [], relation_attributes = ["method_id"], allowed_relation_types = "all", include_relation_ids = False, include_participant_ids = True, include_relation_type = True, include_relation_sources = True, output_1_value_per_attribute = False, output_format = "tabulated", only_selected = False, substitute_node_attribute_if_not_exists = False, include_participant_tags = False, include_relation_tags = False, include_unconnected_nodes=False)
    return


def old_create_human_interactome_using_biana(node_file_prefix="human_nodes", network_files_prefix="human_network", include_string=False): #, load_from_saved_session=False):
    """
	Creates human ppi network files using BIANA
	Node file: node_file_prefix.tsv (as tab seperated values)
	Edge files: network_files_prefix.sif (network in sif format), network_file_prefix_attribute.eda (Cytoscape edge attribute files for each attribute)
    """
    from biana.biana_commands import available_sessions, create_new_session, save_session, load_session
    from biana.BianaObjects import BianaSessionManager
    
    DB_NAME = "test_biana" 
    DB_HOST = "127.0.0.1" 
    DB_USER = "biana_user" 
    DB_PASS = "biana_pass" 
    UNIFICATION_PROTOCOL = "(p2)(noself)(noprevious)uniprot_seqtax_geneid_scoppdb" 

    identifier_description_list = [("taxid", "9606")]
    relation_type_list = ["interaction"] #,"reaction"] # ["interaction","reaction","complex","pathway"]
    level = 0

    #if load_from_saved_session:
    #	load_session(network_files_prefix + "_session.dat")
    #	return

    # Create a new BIANA session
    create_new_session(sessionID="biana_session",dbname=DB_NAME, dbhost=DB_HOST, dbuser=DB_USER, dbpassword=DB_PASS, unification_protocol=UNIFICATION_PROTOCOL)
    objSession = available_sessions["biana_session"]

    # Create a new User Entity Set (group of biomolecules) in this session
    uESet1 = objSession.create_new_user_entity_set( identifier_description_list=identifier_description_list, attribute_restriction_list=[], id_type="embedded", new_user_entity_set_id="User_Entity_Set_1")
    uESet1 = objSession.get_user_entity_set("User_Entity_Set_1")

    print "# of nodes: ", uESet1.getSize()

    if include_string:
	# Duplicate the User Entity Set (group of biomolecules) in this session
	uESet2 = objSession.duplicate_user_entity_set(user_entity_set_id="User_Entity_Set_1", new_user_entity_set_id="User_Entity_Set_2")

	# Fetch relations (except STRING) of biomolecules in 1st set
	objSession.create_network( user_entity_set_id = "User_Entity_Set_1" , level = level, relation_type_list=relation_type_list, relation_attribute_restriction_list=[], include_relations_last_level=True, use_self_relations=True)

	print "# of experimental edges: ", uESet1.getNumberEdges()

	# Fetch relations (only STRING with scores > 700) of biomolecules in 2nd set
	relation_type_list = ["functional_association"]
	#objSession.create_network( user_entity_set_id = "User_Entity_Set_2" , level = level, relation_type_list=relation_type_list, relation_attribute_restriction_list=[("stringscore",">700")], include_relations_last_level=True, use_self_relations=True)
	objSession.create_network( user_entity_set_id = "User_Entity_Set_2" , level = level, relation_type_list=relation_type_list, relation_attribute_restriction_list=[], include_relations_last_level=True, use_self_relations=True)
	print "# of string edges: ", uESet2.getNumberEdges()

	uESet3 = objSession.get_union_of_user_entity_set_list(user_entity_set_list=["User_Entity_Set_1", "User_Entity_Set_2"], include_relations=True, new_user_entity_set_id="User_Entity_Set_3") 
	objSession.remove_user_entity_set(user_entity_set_id="User_Entity_Set_1")
	objSession.remove_user_entity_set(user_entity_set_id="User_Entity_Set_2")

	# Output set details - all nodes
	objSession.output_user_entity_set_details(user_entity_set_id = "User_Entity_Set_3", out_method = open(node_file_prefix+".tsv", "w").write, attributes=["uniprotaccession", "geneid", "ensembl", "genesymbol", "go", "hgnc", "disease", "proteinsequence"], include_level_info=False, include_degree_info=False, level=None, only_selected=False, output_format="tabulated", include_tags_info = False, include_tags_linkage_degree_info=[], substitute_node_attribute_if_not_exists=False, output_1_value_per_attribute=False, include_command_in_rows=False)
	
	# Export all the information of this set and its network to a file in a tab separated format
	#objSession.output_user_entity_set_network(user_entity_set_id = "User_Entity_Set_3", out_method=open("network.out", "w").write, node_attributes = ["uniprotaccession"], participant_attributes = [], relation_attributes=["stringscore"], include_relation_ids=False, include_participant_ids=False, include_relation_type=True, include_relation_sources=True, output_1_value_per_attribute=False, output_format="tabulated", include_participant_tags=False, include_relation_tags=False) 
	objSession.output_user_entity_set_network_in_sif_format(user_entity_set_id = "User_Entity_Set_3", output_path = "./", output_prefix = network_files_prefix, node_attributes = [], participant_attributes = [], relation_attributes=["pubmed", "method_id","stringscore"], output_1_value_per_attribute=False, include_tags=False) 
    else:
	
	# Fetch relations (except STRING) of biomolecules in 1st set
	objSession.create_network( user_entity_set_id = "User_Entity_Set_1" , level = level, relation_type_list=relation_type_list, relation_attribute_restriction_list=[], include_relations_last_level=True, use_self_relations=True)
	print "# of experimental edges: ", uESet1.getNumberEdges()
	## Get tags for each node - seed, connected, whatever...
	#user_entities_to_print = set(available_sessions["biana_session"].get_user_entity_set(user_entity_set_id = "User_Entity_Set_1").get_user_entity_ids(level=0)) | set(available_sessions["biana_session"].get_user_entity_set(user_entity_set_id = "User_Entity_Set_1").get_user_entity_ids_by_linker_degree_cutoff(2))
	#available_sessions["biana_session"].select_user_entities_from_user_entity_set("User_Entity_Set_1", user_entities_to_print, clear_previous_selection = True)
	#objSession.tag_selected_user_entities(user_entity_set_id="User_Entity_Set_1", tag="seed+connected")

	# Output set details - all nodes
	objSession.output_user_entity_set_details(user_entity_set_id = "User_Entity_Set_1", out_method = open(node_file_prefix + ".tsv", "w").write, attributes=["uniprotaccession", "geneid", "ensembl", "genesymbol", "go", "hgnc", "disease", "proteinsequence"], include_level_info=False, include_degree_info=True, level=None, only_selected=False, output_format="tabulated", include_tags_info = True, include_tags_linkage_degree_info=[], substitute_node_attribute_if_not_exists=False, output_1_value_per_attribute=False, include_command_in_rows=False)

	## Export sequences of entries in this set as a FASTA file
	#objSession.output_user_entity_set_sequences_in_fasta(user_entity_set_id="User_Entity_Set_1", only_selected = False, out_method=open("human_fasta.txt",'w').write, type="proteinsequence") #, attributes=["uniprotaccession"])
	
	# Export all the information of this set and its network to a file in a tab separated format
	objSession.output_user_entity_set_network_in_sif_format(user_entity_set_id = "User_Entity_Set_1", output_path = "./", output_prefix = network_files_prefix, node_attributes = [], participant_attributes = [], relation_attributes=["pubmed", "method_id"], output_1_value_per_attribute=False, include_tags=False) 
	available_sessions["biana_session"].output_user_entity_set_network(user_entity_set_id = "User_Entity_Set_1", out_method = open(network_files_prefix + ".tsv", "w").write, node_attributes = [], participant_attributes = [], relation_attributes = ["method_id"], allowed_relation_types = "all", include_relation_ids = False, include_participant_ids = True, include_relation_type = True, include_relation_sources = True, output_1_value_per_attribute = False, output_format = "tabulated", only_selected = False, substitute_node_attribute_if_not_exists = False, include_participant_tags = False, include_relation_tags = False, include_unconnected_nodes=False)
    # Save the session
    save_session("biana_session", network_files_prefix + "_session.dat")
    return


def old_generate_cross_validation_node_score_files(g, node_to_score, seeds, node_file_netzcore_prefix, edge_file_netzcore_relevance_prefix, arff_file_prefix, edge_file_netzcore):
    # seeds = [ v for v, s in node_to_score.iteritems() if s != 0 ] # was assuming that non-seeds should have zero score
    for k, training, test in k_fold_cross_validation(seeds, 10, randomize = True, useSeed = True):
	assign_node_scores(g = g, node_to_score = node_to_score, out_file = node_file_netzcore_prefix+"_%i.txt"%k, as_edge_relevance_score = False, seeds_ignored=test)
	assign_node_scores(g = g, node_to_score = node_to_score, out_file = edge_file_netzcore_relevance_prefix+"_%i.txt"%k, as_edge_relevance_score = True, seeds_ignored=test)
	create_edge_file_from_weight_and_score_files(edge_file_netzcore_relevance_prefix+"_%i.txt"%k, edge_file_netzcore)
	network_utilities.create_arff_file_with_network_metrics(g, node_to_score, training, arff_file_prefix+"_%i.arff"%k)
    return


def old_create_edge_file_from_weight_and_score_files(edge_file_weights, edge_file_scores, out_file):
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


def old_extract_test_data_from_arff_files(arff_file_all, arff_file_filtered, arff_file_training, arff_file_test):
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


def old_generate_cross_validation_test_test_arff_files(k, arff_file_all, arff_file_filtered, arff_file_training_prefix):
    for i in xrange(k):
	extract_test_data_from_arff_files(arff_file_all, arff_file_filtered, arff_file_training_prefix+"_%i.arff"%i, arff_file_training_prefix+"_test_%i.arff"%i)
    return


def old_assign_edge_reliability_scores(g, network_file_method_attribute, network_file_source_attribute, network_file_pubmed_attribute, out_file):
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


if __name__ == "__main__":
    pass

