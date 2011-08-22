#########################################################################
# Methods for analyzing output files from scoring methods
#
# eg 13/08/2009
#########################################################################

from biana.utilities import graph_utilities as network_utilities
from funcassociate import client
from biana.utilities import TsvReader


def score_mcl(node_scores_file, network_file, output_scores_file, module_file, default_non_seed_score):
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)
    #modules = get_modules_of_graph(g, "mcl", inflation=2) # if edge weight based clustering is desired
    seeds, nodes = get_seeds_from_node_scores_file(node_scores_file, default_non_seed_score)
    modules = get_modules_from_file(module_file)
    f = open(output_scores_file, 'w')
    node_to_score = {}
    #selected = set()
    for module in modules:
	module = set(module)
	#common = module&seeds
	#if 100*float(len(common))/len(module) > threshold:
	    #selected |= module
	#score = float(len(common))/len(module)
	#n = len(module)-len(common)
	#if n == 0:
	#    continue
	#score = 1.0/n
	for node in module:
	    #node_to_score[node] = score
	    neighbors = set(g.neighbors(node))
	    common = neighbors & module
	    if node in common:
		common.remove(node)
	    #if len(common) == 0:
	    #	continue
	    score = float(len(common&seeds)) / len(module)
	    node_to_score[node] = score
    for node in nodes:
	if node in node_to_score:
	    f.write("%s\t%f\n" % (node, node_to_score[node]))
	else:
	    f.write("%s\t0.0\n" % node)
    f.close()
    return

def get_seed_and_all_nodes_and_ids(node_scores_file, node_mapping_file, default_non_seed_score):
    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)
    dummy, dummy, initial_node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)
    seeds = set()
    seed_ids = set()
    all_ids = set()
    for node in initial_node_to_score:
	if initial_node_to_score[node] > default_non_seed_score:
	    seeds.add(node)
	    if node in id_to_mapped_ids:
		seed_ids.add(id_to_mapped_ids[node][0])
	if node in id_to_mapped_ids:
	    all_ids.add(id_to_mapped_ids[node][0])
    return seeds, seed_ids, all_ids

def get_seeds_from_node_scores_file(node_scores_file, default_non_seed_score):
    nodes, dummy, initial_node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)
    seeds = set()
    for node in initial_node_to_score:
	if initial_node_to_score[node] > default_non_seed_score:
	    seeds.add(node)
    return seeds, nodes

def get_modules_from_file(output_file):
    f = open(output_file)
    modules = []
    for line in f:
	words = line.strip().split("\t")
	modules.append(words)
    f.close()
    return modules

def get_modules_of_graph(sub_graph, module_detection_type, output_file=None, inflation=1.7):
    if module_detection_type == "connected":
	modules = network_utilities.get_connected_components(sub_graph, return_as_graph_list=True)
    elif module_detection_type == "mcl":
	from os import system
	if output_file is None:
	    temp_file_name = ".temp_module_file.txt89734234"
	else:
	    temp_file_name = output_file
	f = open(temp_file_name, 'w')
	nodes = set()
	for node1, node2, data in sub_graph.edges(data=True):
	    nodes.add(node1)
	    nodes.add(node2)
	    f.write("%s\t%s\t%f\n" % (node1, node2, data))
	for node in sub_graph.nodes():
	    if node not in nodes:
		f.write("%s\n" % node)
	f.close()
	# Optimum inflation parameter was 1.7-1.8 in a recent comparison paper
	system("mcl %s --abc -I %f -o %s 2>> %s" % (temp_file_name, inflation, temp_file_name + ".mcl", temp_file_name + ".err"))
	f = open(temp_file_name + ".mcl")
	modules = []
	for line in f:
	    words = line.strip().split("\t")
	    modules.append(words)
	f.close()
    else:
	raise ValueError("Unrecognized module detection type")
    #print len(modules), map(len, modules)
    return modules

def get_seed_function_enrichment_in_modules(modules, sub_graph, seeds, seed_ids, all_ids, id_to_mapped_ids, output_method, go_ids, id_type = "genesymbol", specie = "Homo sapiens", p_value_cutoff = 0.05): 
    """
	Get functional enrichment of seed within modules 
    """
    n = len(seed_ids)
    n_seed, n_all = 0, 0
    n_go_term, n_go_term_all = 0, 0
    n_module = 0 
    # FuncAssoc row: n_gene n_all LOD p_val adj_p_val go_id go_name
    #response = check_functional_enrichment(list(seed_ids), list(all_ids), id_type, None, specie = specie, mode = "unordered")
    #for row in response:
    #	print row 
    #go_ids = set([ row[6] for row in response if float(row[5]) <= p_value_cutoff ])
    #print len(go_ids), len(response)
    i = 0
    all_go_terms_in_modules = set()
    seed_go_terms_in_modules = set()
    for module in modules:
	if len(module) == 1:
	    continue
	ids_in_module = set([ id_to_mapped_ids[node][0] for node in module if node in id_to_mapped_ids])
	k = len(seed_ids & ids_in_module)
	if k < 2:
	    continue
	N = len(ids_in_module)
	i += 1
	n_seed += k
	n_all += N
	output_method("%s\n" % ", ".join(ids_in_module))
	response_module = check_functional_enrichment(list(ids_in_module), list(all_ids), id_type, output_method, specie = specie, mode = "unordered")
	for row in response_module:
	    if float(row[5]) <= p_value_cutoff:
		if row[6] in go_ids:
		    n_go_term += 1 
		    seed_go_terms_in_modules.add(row[6])
		n_go_term_all += 1
		all_go_terms_in_modules.add(row[6])
    n_module = i
    if n_go_term_all == 0:
	ratio = 0
    else:
	ratio = n_go_term/float(n_go_term_all)
    #return (ratio, n_seed, n_all, n_module, n)
    #return (ratio, n_go_term, n_go_term_all, n_module, len(go_ids))
    return (ratio, len(seed_go_terms_in_modules), len(all_go_terms_in_modules), n_module, len(go_ids))


def get_connected_seed_enrichment_in_modules(modules, sub_graph, seeds, seed_ids, all_ids, id_to_mapped_ids, output_method = None, enrichment_type = "seed"): #"non_seed_in_between_connections"):
    """
	Calculates enrichment of modules in terms of seeds according to one of the following metrics:
	seed, non_seed_connections, non_seed_in_between_connections, seed_excluding_connections

    """
    from scipy.stats import hypergeom
    M = len(all_ids) 
    n = len(seed_ids)
    i = 0
    n_seed = 0
    n_all = 0
    n_module = 0 
    e_non_seed = 0
    e_all = 0
    n_non_seed_in_between = 0
    n_in_between = 0
    for module in modules:
	if len(module) == 1:
	    continue
	ids_in_module = set([ id_to_mapped_ids[node][0] for node in module if node in id_to_mapped_ids])
	k = len(seed_ids & ids_in_module)
	if k < 2:
	    continue
	N = len(ids_in_module)
	if output_method is not None:
	    p_val = sum(hypergeom.pmf(x, M, n, N) for x in range(k,n+1))
	    output_method("Module_%d: %d %d %d %e\n" % (i, k, N, n, p_val))
	i += 1
	#f = open(".%d.mcl" % i, 'w') 
	#[ f.write("%s\n" % node) for node in module ]
	#f.close()
	n_seed += k
	n_all += N
	if enrichment_type == "non_seed_in_between_connections":
	    sub_sub_graph = network_utilities.get_subgraph(sub_graph, module)
	    seeds_in_module = seeds & set(module)
	    if True: # based on shortest paths between all nodes
		for a, node1 in enumerate(module):
		    for b, node2 in enumerate(module):
			if a < b:
			    path = network_utilities.get_shortest_path_between(sub_sub_graph, node1, node2)
			    if path:
				n_non_seed_in_between += len(set(path)-seeds_in_module)
				n_in_between += len(path) - 2 
	    else: # based on shortest paths between all seeds
		for a, seed1 in enumerate(seeds_in_module):
		    for b, seed2 in enumerate(seeds_in_module):
			if a < b:
			    path = network_utilities.get_shortest_path_between(sub_sub_graph, seed1, seed2)
			    if path:
				n_non_seed_in_between += len(set(path)-seeds_in_module)
				n_in_between += len(path) - 2 
			    #print path
			    #print n_non_seed_in_between, n_in_between
	elif enrichment_type == "non_seed_connections":
	    sub_sub_graph = network_utilities.get_subgraph(sub_graph, module)
	    connections = sub_sub_graph.edges(seeds)
	    non_self_connections = [ connection for connection in connections if connection[0] != connection[1]]
	    non_seed_connections = [ connection for connection in non_self_connections if connection[0] not in seeds or connection[1] not in seeds ]
	    e_non_seed += len(non_seed_connections)
	    e_all += len(non_self_connections)
	elif enrichment_type == "seed_excluding_connections":
	    connections = sub_graph.edges(module)
	    non_self_connections = [ connection for connection in connections if connection[0] != connection[1]]
	    non_seed_connections = [ connection for connection in non_self_connections if connection[0] not in seeds and connection[1] not in seeds ]
	    e_non_seed += len(non_seed_connections)
	    e_all += len(non_self_connections)
    n_module = i

    if enrichment_type == "non_seed_in_between_connections":
	if n_in_between == 0:
	    ratio = 0
	else:
	    ratio = n_non_seed_in_between/float(n_in_between)
    elif enrichment_type in ("non_seed_connections", "seed_excluding_connections"):
	if e_all == 0:
	    ratio = 0
	else:
	    ratio = e_non_seed/float(e_all)
    elif enrichment_type == "seed":
	ratio = n_seed/float(n)

    return (ratio, n_seed, n_all, n_module, n)

def check_connected_seed_enrichment_in_network(network_file, node_scores_file, node_mapping_file, id_type, output_method, default_non_seed_score, module_detection_type = "connected", go_ids = None, specie = "Homo sapiens", p_value_cutoff = 0.05):
    """
	Get number of seeds included in one level neighborhood of seeds
	module_detection_type: connected | mcl
    """
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)
    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)
    dummy, dummy, initial_node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)

    seeds = set()
    seed_ids = set()
    for node in g.nodes():
	if initial_node_to_score[node] > default_non_seed_score:
	    seeds.add(node)
	    if node in id_to_mapped_ids:
		seed_ids.add(id_to_mapped_ids[node][0])

    all_ids = set()
    non_mapped_nodes = set()
    for node in initial_node_to_score:
	if node in id_to_mapped_ids:
	    all_ids.add(id_to_mapped_ids[node][0])
	else:
	    non_mapped_nodes.add(node)
    #print len(non_mapped_nodes), "nodes do not exist in the mapping file"

    sub_graph = network_utilities.get_subgraph(g, seeds)
    modules = get_modules_of_graph(sub_graph, module_detection_type)
    if go_ids is None:
	(ratio, n_seed, n_all, n_module, n) = get_connected_seed_enrichment_in_modules(modules, sub_graph, seeds, seed_ids, all_ids, id_to_mapped_ids, output_method)
    else:
	(ratio, n_seed, n_all, n_module, n) = get_seed_function_enrichment_in_modules(modules, sub_graph, seeds, seed_ids, all_ids, id_to_mapped_ids, output_method, go_ids, id_type = id_type, specie = specie, p_value_cutoff = p_value_cutoff)
    return (ratio, n_seed, n_all, n_module, n)


def check_connected_seed_enrichment_of_neighbors_in_network(network_file, node_scores_file, node_mapping_file, id_type, output_method, default_non_seed_score, module_detection_type = "connected", go_ids = None, specie = "Homo sapiens", p_value_cutoff = 0.05):
    """
	Get number of seeds included in one level neighborhood of seeds
	module_detection_type: connected | mcl
    """
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)
    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)
    dummy, dummy, initial_node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)

    seeds = set()
    seed_ids = set()
    for node in g.nodes():
	if initial_node_to_score[node] > default_non_seed_score:
	    seeds.add(node)
	    if node in id_to_mapped_ids:
		seed_ids.add(id_to_mapped_ids[node][0])

    neighbors = set()
    for seed in seeds:
	neighbors |= set(g.neighbors(seed))

    selected_nodes = seeds|neighbors
    selected_ids = set(id_to_mapped_ids[node][0] for node in selected_nodes if node in id_to_mapped_ids)
    all_ids = set()
    non_mapped_nodes = set()
    for node in initial_node_to_score:
	if node in id_to_mapped_ids:
	    all_ids.add(id_to_mapped_ids[node][0])
	else:
	    non_mapped_nodes.add(node)
    #print len(non_mapped_nodes), "nodes do not exist in the mapping file"

    sub_graph = network_utilities.get_subgraph(g, selected_nodes)
    modules = get_modules_of_graph(sub_graph, module_detection_type)
    if go_ids is None:
	(ratio, n_seed, n_all, n_module, n) = get_connected_seed_enrichment_in_modules(modules, sub_graph, seeds, seed_ids, all_ids, id_to_mapped_ids, output_method)
    else:
	(ratio, n_seed, n_all, n_module, n) = get_seed_function_enrichment_in_modules(modules, sub_graph, seeds, seed_ids, all_ids, id_to_mapped_ids, output_method, go_ids, id_type = id_type, specie = specie, p_value_cutoff = p_value_cutoff)
    return (ratio, n_seed, n_all, n_module, n) 

def check_connected_seed_enrichment_of_modules_of_given_nodes(nodes, network_file, node_scores_file, node_mapping_file, id_type, output_method, default_non_seed_score, module_detection_type = "connected", go_ids = None, specie = "Homo sapiens", p_value_cutoff = 0.05):
    """
	Get number of seeds included in modules (connected nodes) identified by given nodes
    """
    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)
    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)
    dummy, dummy, initial_node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)

    seeds = set()
    seed_ids = set()
    for node in g.nodes():
	if initial_node_to_score[node] > default_non_seed_score:
	    seeds.add(node)
	    if node in id_to_mapped_ids:
		seed_ids.add(id_to_mapped_ids[node][0])

    selected_nodes = nodes
    selected_ids = set(id_to_mapped_ids[node][0] for node in selected_nodes if node in id_to_mapped_ids)
    all_ids = set()
    non_mapped_nodes = set()
    for node in initial_node_to_score:
	if node in id_to_mapped_ids:
	    all_ids.add(id_to_mapped_ids[node][0])
	else:
	    non_mapped_nodes.add(node)
    #print len(non_mapped_nodes), "nodes do not exist in the mapping file"

    sub_graph = network_utilities.get_subgraph(g, selected_nodes)
    modules = get_modules_of_graph(sub_graph, module_detection_type)
    if go_ids is None:
	(ratio, n_seed, n_all, n_module, n) = get_connected_seed_enrichment_in_modules(modules, sub_graph, seeds, seed_ids, all_ids, id_to_mapped_ids, output_method)
    else:
	(ratio, n_seed, n_all, n_module, n) = get_seed_function_enrichment_in_modules(modules, sub_graph, seeds, seed_ids, all_ids, id_to_mapped_ids, output_method, go_ids, id_type = id_type, specie = specie, p_value_cutoff = p_value_cutoff)
    return (ratio, n_seed, n_all, n_module, n) 

def check_connected_seed_enrichment_of_high_scoring_modules(network_file, output_scores_file, node_scores_file, node_mapping_file, cutoff, id_type, output_method, default_non_seed_score, module_detection_type = "connected", go_ids = None, specie = "Homo sapiens", p_value_cutoff = 0.05):
    """
	Get number of seeds included in modules (connected nodes) identified by the scoring method
    """
    selected_ids, all_ids, seed_ids, selected_nodes = get_top_scoring_node_ids_at_given_cutoff(output_scores_file, node_scores_file, node_mapping_file, cutoff, id_type, default_non_seed_score, exclude_seeds = False, one_gene_per_node = True)

    g = network_utilities.create_network_from_sif_file(network_file, use_edge_data=True)
    dummy, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = output_scores_file, store_edge_type = False)
    sub_graph = network_utilities.get_subgraph(g, selected_nodes)
    modules = get_modules_of_graph(sub_graph, module_detection_type)

    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)

    dummy, dummy, initial_node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)
    seeds = set()
    for node in g.nodes():
	if initial_node_to_score[node] > default_non_seed_score:
	    seeds.add(node)

    seed_ids = set(seed_ids)
    all_ids = set(all_ids)
    if go_ids is None:
	(ratio, n_seed, n_all, n_module, n) = get_connected_seed_enrichment_in_modules(modules, sub_graph, seeds, seed_ids, all_ids, id_to_mapped_ids, output_method)
    else:
	(ratio, n_seed, n_all, n_module, n) = get_seed_function_enrichment_in_modules(modules, sub_graph, seeds, seed_ids, all_ids, id_to_mapped_ids, output_method, go_ids, id_type = id_type, specie = specie, p_value_cutoff = p_value_cutoff)
    return (ratio, n_seed, n_all, n_module, n) #min(vals)

def check_functional_enrichment_of_high_scoring_modules(network_file, module_detection_type, output_scores_file, node_scores_file, node_mapping_file, cutoff, id_type, output_method, default_non_seed_score = 0, exclude_seeds = False, specie = "Homo sapiens", mode = "unordered"):
    """
	Check functional enrichment of highest scoring modules at given cutoff
    """
    g = network_utilities.create_network_from_sif_file(network_file)

    dummy, dummy, initial_node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)

    dummy, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = output_scores_file, store_edge_type = False)
    ids, n, i = get_top_scoring_nodes_at_given_cutoff(node_to_score, cutoff)

    if module_detection_type != "greedy":
	raise ValueError("Only greedy all highest neighbor node inclusion is supported!")

    sub_graph = network_utilities.get_subgraph(g, ids)
    modules = network_utilities.get_connected_components(sub_graph, return_as_graph_list=True)
    #print "NetworkX way:"
    print len(modules), map(len, modules)

    #modules = get_high_scoring_modules(g, node_to_score, ids)
    #print "Handcrafted way:"
    #print len(modules), map(len, modules)

    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)

    all_ids = reduce(lambda x,y: x+y, id_to_mapped_ids.values())
    for module in modules:
	if len(module) == 1:
	    continue
	if exclude_seeds:
	    ids_in_module = [ node for node in module if initial_node_to_score[node] <= default_non_seed_score ]
	else:
	    ids_in_module = module
	selected_ids = reduce(lambda x,y: x+y, [ id_to_mapped_ids[node] for node in ids_in_module if node in id_to_mapped_ids])
	output_method("Module: %d genes among %d\n" % (len(selected_ids), len(set(all_ids))))
	output_method("%s\n" % ", ".join(selected_ids))
	check_functional_enrichment(selected_ids, all_ids, id_type, output_method, specie = specie, mode = mode)
    return


def get_high_scoring_modules(g, node_to_score, ids):
    """
	Can use ids instead of node_to_score & min_score, helper function can be removed and made inline 
    """
    # Using MCL for modules
    return get_modules_of_graph(g, "mcl")
    
    # Below not tested, likely to yield in the connected component of the highest scoring node
    def get_high_scoring_neighbors(v, g, node_to_score, min_score):
	neighbors = set()
	for u in g.neighbors(v):
	    if node_to_score[u] >= min_score:
		neighbors.add(u)
	return neighbors

    min_score = node_to_score[ids[25]]
    modules = []
    current_module = set()
    ids_included_in_modules = set()
    #for u,v in g.edges_iter():
    for v in ids:
	if v in ids_included_in_modules:
	    continue
	else:
	    modules.append(current_module)
	    current_module = set()
	current_module.add(v)
	ids_included_in_modules.add(v) 
	current_module_copy = set() | current_module
	for u in current_module_copy:
	    neighbors = get_high_scoring_neighbors(u, g, node_to_score, min_score)
	    for n in neighbors:
		if n in ids_included_in_modules:
		    continue
		ids_included_in_modules.add(n) 
		current_module.add(n)
    return modules


def get_id_to_mapped_id_mapping(node_mapping_file):
    reader = TsvReader.TsvReader(node_mapping_file, inner_delim = ",")
    columns, id_to_mapped_ids = reader.read(fields_to_include = None, merge_inner_values = True)
   
    id_to_mapped_ids_formatted = {}
    for node, vals in id_to_mapped_ids.iteritems():
	vals = reduce(lambda x,y: x+y, vals)
	if "-" in vals:
	    vals.remove("-")
	if len(vals) < 1:
	    continue
	id_to_mapped_ids_formatted[node] = vals

    return id_to_mapped_ids_formatted


def get_top_scoring_node_ids_at_given_cutoff(output_scores_file, node_scores_file, node_mapping_file, cutoff, id_type, default_non_seed_score=0, exclude_seeds=False, one_gene_per_node=True):
    """
	Get ids of highest scoring nodes at given cutoff (if ends with % taken as percentage, otherwise taken as score cutoff)
    """
    dummy, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = output_scores_file, store_edge_type = False)

    ids, n, i = get_top_scoring_nodes_at_given_cutoff(node_to_score, cutoff)

    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)
   
    selected_ids = []
    all_ids = []
    seed_ids = []
    selected_node_ids = []

    dummy, dummy, initial_node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)

    for node in ids:
	if exclude_seeds:
	    if initial_node_to_score[node] > default_non_seed_score:
		#print id_to_mapped_ids[node][0]
		continue
	selected_node_ids.append(node)
	if node not in id_to_mapped_ids:
	    continue
	if one_gene_per_node:
	    vals = [ id_to_mapped_ids[node][0] ]
	else:
	    vals = id_to_mapped_ids[node]
	selected_ids.extend(vals)

    non_mapped_nodes = set()
    for node in initial_node_to_score:
    #for val_list in id_to_mapped_ids.values():
	if node in id_to_mapped_ids:
	    val_list = id_to_mapped_ids[node]
	    if one_gene_per_node:
		all_ids.extend([val_list[0]])
	    else:
		all_ids.extend(val_list)
	else:
	    non_mapped_nodes.add(node)
    #print len(non_mapped_nodes), "nodes do not exist in the mapping file"

    for node, score in initial_node_to_score.iteritems():
	if score > default_non_seed_score:
	    if node not in id_to_mapped_ids:
		continue
	    if one_gene_per_node:
		vals = [ id_to_mapped_ids[node][0] ]
	    else:
		vals = id_to_mapped_ids[node]
	    seed_ids.extend(vals)

    return selected_ids, all_ids, seed_ids, selected_node_ids


def check_functional_enrichment_at_given_cutoff(output_scores_file, node_scores_file, node_mapping_file, cutoff, id_type, output_method, default_non_seed_score=0, exclude_seeds=False, specie = "Home sapiens", mode="unordered"):
    """
	Check functional enrichment of highest scoring nodes at given cutoff
    """
    selected_ids, all_ids, seed_ids, selected_node_ids = get_top_scoring_node_ids_at_given_cutoff(output_scores_file, node_scores_file, node_mapping_file, cutoff, id_type, default_non_seed_score, exclude_seeds)
    #print len(selected_ids), len(all_ids)
    output_method("%d gene names/ids among %d\n" % (len(selected_ids), len(all_ids)))
    check_functional_enrichment(selected_ids, all_ids, id_type, output_method, specie = specie, mode = mode)
    return


def check_functional_enrichment(subset_gene_ids, gene_ids, id_type, output_method, specie = "Homo sapiens", mode = "unordered", request_info=False, tex_format=False):
    """
	Check GO functional enrichment using funcassociate web service
	gene_ids is a list of gene symbols (without whitespace) or gene ids
	id_type
    """
    if id_type == "geneid":
	id_type = "entrezgene"
    elif id_type == "genesymbol":
	id_type = "hgnc_symbol"
    elif id_type == "uniprotaccession":
	id_type = "uniprot_accession"
    elif id_type == "uniprotentry":
	id_type = "uniprot_id"
    elif id_type == "sgd":
	id_type = "sgd_systematic"
    else:
	raise ValueError("Unrecognized id_type: %s" % id_type)

    reps = 1500
    client_funcassociate = client.FuncassociateClient()
    response = client_funcassociate.functionate(query = subset_gene_ids,
                             species = specie,
                             namespace = id_type,
                             genespace = gene_ids,
                             mode = mode,
                             reps = reps)

    if output_method is None:
	return response["over"]

    #headers = ["N", "M", "X", "LOD", "P", "P_adj", "attrib ID", "attrib name"]
    headers = [ "# of high scoring genes", "# of total genes in high scoring subquery genes", "# of total genes", "Log of odds ratio", "P-value", "Adjusted p-value", "GO term ID", "Go term name" ]

    #if mode == "unordered":
    #	headers.pop(1)
    headers.pop(1) # Now that column is always present independent of the mode
    if tex_format:
	output_method("%s\\\\\n" % " & ".join(headers))
    else:
	output_method("%s\n" % "\t".join(headers))

    zero = "< %f" % (1.0/float(reps))

    for row in response["over"]:
	if mode == "unordered":
	    row = row[:1] + row[2:] #row.pop(1)
        if row[4] is 0:
            row[4] = zero
	if mode == "unordered":
	    interval = range(2,5)
	else:
	    interval = range(3,6)
	#print row
	for i in interval:
	    if isinstance(row[i], str) and row[i].startswith("<"):
		#print row[i]
		val = float(row[i].lstrip("<"))
		if tex_format:
		    row[i] = "$<$%.5f" % val
		else:
		    row[i] = "<%.5f" % val
	    else:
		row[i] = "%.5f" % row[i]
	if tex_format:
	    output_method("%s\\\\\n" % " & ".join(map(str, row)))
	else:
	    output_method("%s\n" % "\t".join(map(str, row)))

    if request_info:
	output_method("\nREQUEST INFO\n")
	info = response["request_info"]
	for k in info.keys():
	    output_method("%s: %s\n" % (k, info[k]))
    return response["over"]


def record_performance_AUC_in_log_file(absolute_dir, log_file, title):
    f = open(absolute_dir + "auc.txt")
    line = f.readline()
    words = line.strip().split()
    mean, sigma = words[1].strip("\""), words[2].strip("\"")
    f.close()
    f = open(log_file, "a")
    f.write("%s\t%s\t+/- %s" % (title, mean, sigma)) 
    f.close()
    return


def record_performance_coverage_in_log_file(log_file, coverages):
    f = open(log_file, "a")
    f.write("\t%s\n" % "\t".join(str(i) for i in coverages)) 
    f.close()
    return


def create_tex_script(fileName, absolute_dir, title):
    f = open(fileName, "w")
    f.write("\\frame {\n") 
    f.write("\\frametitle{%s}\n" % title) 
    #\vspace{-0.35cm}
    f.write("\\begin{figure}\n") 
    f.write("\t\\includegraphics[scale=0.9]{%sperformance.eps}\n" % absolute_dir)
    f.write("\\end{figure}\n") 
    f.write("}\n") 
    f.close()
    return

def create_R_script(fileName, absolute_dir, title, only_auc=False):
    f = open(fileName, "w")
    f.write("library(ROCR)\n") 
    if only_auc:
	f.write("v<-read.table(\"%spredictions.txt\")\n" % absolute_dir) 
	f.write("l<-read.table(\"%slabels.txt\")\n" % absolute_dir)
	f.write("pred<-prediction(v, l)\n")
	f.write("perfAUC<-performance(pred, \"auc\")\n")
	f.write("e=c(); n=c(); x=0; for ( i in perfAUC@y.values ) { x<-x+1;  e[x] <- i; n[x]<-x }\n")
	f.write("sink(\"%sauc.txt\", append=TRUE, split=TRUE)\n" % absolute_dir)
	f.write("paste(format(mean(e), digits=3), format(sd(e), digits=3), sep=\" \")\n") 
	f.write("sink()\n")
	f.close()
	return
    f.write("f<-function(perf) { d<-rep(0, times=length(perf@x.values[[1]]))\n")
    f.write("for(vals in perf@x.values) { d<-d+vals; }\n")
    f.write("for(i in 1:length(d)) { if(d[i] == Inf) { d[i]<-0; }; }\n")
    f.write("d<-d/length(perf@x.values); d<-max(d); return(d); }\n")
    f.write("extractvals<-function(perf) { x=c(); y=c(); for ( i in 1:length(perf@y.values) ) { x[i]<-perf@x.values[i]; y[i]<-perf@y.values[i] }; list(x,y); }\n")
    f.write("v<-read.table(\"%spredictions.txt\")\n" % absolute_dir) 
    f.write("l<-read.table(\"%slabels.txt\")\n" % absolute_dir)
    f.write("pred<-prediction(v, l)\n")
    #if image_type = "eps":
    f.write("#postscript(\"%sperformance.eps\", width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = \"special\", title = \"%s\")\n" % (absolute_dir, title))
    #else:
    #	f.write("bitmap(\"%sperformance.jpg\", res = 1200, height = 6, width = 6, type = \"jpeg\", horizontal = FALSE, onefile = FALSE, paper = \"special\", title = \"%s\")" % (absolute_dir, title)) 
    f.write("par(mfrow=c(2,2))\n")
    f.write("perfROC<-performance(pred, \"tpr\", \"fpr\")\n")
    f.write("plot(perfROC, lwd=2, col=2, xlab=\"False Positive Rate\", ylab=\"True Positive Rate\", main=\"ROC curve\", plotCI.col=2, avg=\"vertical\", spread.estimate=\"stddev\", show.spread.at=seq(0,1,by=0.20))\n")
    #f.write("title(\"ROC\")\n")
    f.write("legend(\"bottomright\", c(\"(Avg. over xval folds)\"), lty=c(1), col=c(2))\n") 
    #f.write("legend(\"bottomright\", c(\"(Avg. over xval folds)\", paste(\"AUC:\", format(mean(y), digits=2), sep=\" \")), lty=c(1,0), col=c(2,1))\n") 
    f.write("perfPPV<-performance(pred, \"ppv\")\n")
    f.write("perfSens<-performance(pred, \"sens\")\n")
    f.write("d<-f(perfPPV)\n")
    #f.write("plot(perfPPV, lwd=2, col=2, ylab=\"PPV & Sens\", main=\"PPV vs Sensitivity\", plotCI.col=2, avg=\"vertical\", spread.estimate=\"stddev\", show.spread.at=seq(0,1,by=0.20))\n")        
    f.write("plot(perfPPV, lwd=2, col=2, ylab=\"PPV & Sens\", main=\"PPV vs Sensitivity\", plotCI.col=2, avg=\"vertical\", spread.estimate=\"stddev\", show.spread.at=seq(0,d,by=d/6))\n")        
    f.write("d<-f(perfSens)\n")
    f.write("plot(perfSens, lwd=2, col=3, plotCI.col=3, avg=\"vertical\", spread.estimate=\"stddev\", show.spread.at=seq(0,d,by=d/6), add=TRUE)\n") 
    f.write("vals<-extractvals(performance(pred, \"prbe\")); x<-unlist(vals[1]); y<-unlist(vals[2])\n")
    f.write("legend(\"bottomright\", c(\"PPV\", \"Sens\", paste(\"(\", format(mean(x), digits=2), format(mean(y), digits=2), \")\", sep=\" \")), lty=c(1,1,0), col=c(2,3,1))\n") 
    #f.write("vals<-extractvals(performance(pred, \"mat\")); x<-unlist(vals[1]); y<-unlist(vals[2])\n")
    f.write("perfAUC<-performance(pred, \"auc\")\n")
    f.write("e=c(); n=c(); x=0; for ( i in perfAUC@y.values ) { x<-x+1;  e[x] <- i; n[x]<-x }; barplot(e, names=n, ylim=c(0,1),ylab= \"AUC\",xlab=\"Fold\", main=\"Area under ROC curve (AUC)\")\n")
    f.write("legend(\"topright\", c(paste(\"(Avg: \", format(mean(e), digits=3), \")\",sep=\"\")), lty=c(), col=c())\n") 
    f.write("sink(\"%sauc.txt\", append=TRUE, split=TRUE)\n" % absolute_dir)
    f.write("paste(format(mean(e), digits=3), format(sd(e), digits=3), sep=\" \")\n") 
    f.write("sink()\n")
    f.write("perfRMSE<-performance(pred, \"rmse\")\n")
    f.write("e=c(); n=c(); x=0; for ( i in perfRMSE@y.values ) { x<-x+1;  e[x] <- i; n[x]<-x }; barplot(e, names=n, ylim=c(0,1),ylab= \"RMSE\",xlab=\"Fold\"); title(\"Root mean square error (RMSE)\")\n")
    f.write("legend(\"topright\", c(paste(\"(Avg: \", format(mean(e), digits=3), \")\", sep=\"\")), lty=c(), col=c())\n") 
    f.write("mtext(\'%s\', outer=TRUE, line=-1)\n" % title) 
    f.write("dev.off()\n")
    f.write("perfFPR<-performance(pred, \"fpr\")\n")
    f.write("plot(perfFPR, lwd=2, col=3, plotCI.col=3, avg=\"vertical\")\n")
    f.write("level<-0.05\nepsilon<-0.0001\n")
    f.write("m<-c(); for(i in 1:5) { m[i]<-mean(perfFPR@x.values[[i]][perfFPR@y.values[[i]]<level+epsilon & perfFPR@y.values[[i]]>level-epsilon])}\n")
    f.write("sink(\"%scutoff.txt\", append=TRUE, split=TRUE)\n" % absolute_dir)
    f.write("paste(format(mean(m), digits=3), format(sd(m), digits=3), sep=\" \")\n") 
    f.write("sink()\n")
    f.close()
    #os.system("R CMD BATCH %s" % "*.R") 
    return

def create_ROCR_files(list_node_scores_and_labels, file_predictions, file_labels):
    """
	list_node_scores_and_labels: list of node (score, label) tuple (corresponding to each validation node) list (corresponding to xval fold)
    """
    fileOut = open(file_predictions, 'w')
    fileOut2 = open(file_labels, 'w')
    firstTime = True
    for i, node_scores_and_labels in enumerate(zip(*list_node_scores_and_labels)):
	if i == 0:
            for j in xrange(len(node_scores_and_labels)):
                fileOut.write("\tFold_" + str(j+1))  
                fileOut2.write("\tFold_" + str(j+1))  
            fileOut.write("\n")  
            fileOut2.write("\n")  
        fileOut.write("%d"%i)
        fileOut2.write("%d"%i)
        for (score, label) in node_scores_and_labels:
            fileOut.write("\t" + str(score))
            fileOut2.write("\t" + str(label))
        fileOut.write("\n")
        fileOut2.write("\n")
    fileOut.close()
    fileOut2.close()
    return


def get_validation_node_scores_and_labels(file_result, file_seed_test_scores, file_node_scores, n_random_negative_folds = None, default_score = 0, replicable = 123, candidates_file = None):
    """
	Returns a list of scores and labels [ ([0-1], [01]) ] for validation
	file_result: File to parse output scores 
	file_seed_test_scores: File to parse test seeds
	file_node_scores: File to parse all non seeds
	n_random_negative_folds: Number of non-seed scores to be averaged to be assigned as negative instance
				 If None calculated to cover as much as non-seed scores as possible
				 If 0 all negative data is used
	default_score: All nodes that have a higher score than this score in file_node_scores will be considered as seeds
    """
    dictNodeResult, setNodeTest, non_seeds = get_values_from_files_for_performance_metric_counts(file_result, file_seed_test_scores, file_node_scores, default_score, candidates_file=None)
    node_validation_data = [ (dictNodeResult[id], 1) for id in setNodeTest ]

    if candidates_file is not None:
	setCandidates, setDummy, dictCandidates, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = candidates_file, store_edge_type = False)
	dictNodeResult = dictNodeResult.fromkeys(setCandidates)

    if n_random_negative_folds == 0:
	#node_validation_data.extend([(dictNodeResult[id], 0) for id in non_seeds ])
	node_validation_data.extend([(dictNodeResult[id], 0) for id in set(dictNodeResult.keys()) & set(non_seeds) ])
    else:
	n_actual_folds = 0
	negative_sample_size = len(setNodeTest)
	negative_scores = [ 0 ] * negative_sample_size 
	for sample in generate_samples_from_list_without_replacement(non_seeds, negative_sample_size, n_random_negative_folds, replicable = replicable):
	    for i, id in enumerate(sample):
		negative_scores[i] += dictNodeResult[id] 
	    n_actual_folds += 1
	node_validation_data.extend(map(lambda x: (x/n_actual_folds, 0), negative_scores))
    return node_validation_data


def generate_samples_from_list_without_replacement(elements, sample_size, n_folds = None, replicable = None):
    """
	Iteratively returns (yields) n_folds sublists of elements with a size of sample_size 
	n_folds: If None calculated to cover as much elements as possible
	replicable: If not None uses this replicable as the seed for random
    """
    from random import shuffle, seed
    if replicable is not None:
	seed(replicable)
    shuffle(elements)
    if n_folds is None:
	from math import ceil
	#n_folds = len(elements) / sample_size
	n_folds = int(ceil(float(len(elements)) / sample_size))
    for i in range(n_folds):
	if (i+1)*sample_size < len(elements):
	    yield elements[i*sample_size:(i+1)*sample_size]
	else:
	    yield elements[i*sample_size:]
    return


def get_non_seeds_from_node_scores_file(file_node_scores, default_score = 0):
    """
	file_node_scores: initial node/seed scores
    """
    setNode, setDummy, dictNode, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_node_scores, store_edge_type = False)
    non_seeds = [ id for id, score in dictNode.iteritems() if score <= default_score ]
    return non_seeds

def calculate_performance_metric_counts_using_navlakha(dictNodeResult, setNodeTest, score_threshold, non_seeds = None):
    (nP, nN) = (0.0, 0.0)
    for id in setNodeTest: 
	if non_seeds is not None:
	    if id not in non_seeds:
		continue
	score = dictNodeResult[id]
	if score >= score_threshold:
	    nP += 1
	else:
	    nN += 1
    return (nP, nN)

def calculate_performance_metric_counts_using_candidates(dictNodeResult, setNodeTest, non_seeds, score_threshold):
    (nTP, nFP, nFN, nTN) = (0.0, 0.0, 0.0, 0.0)
    for id, score in dictNodeResult.iteritems(): # if candidates based - for each candidate
        if id in setNodeTest: # in the initial association file
            if score >= score_threshold:
                nTP += 1
            else:
                nFN += 1
	else: 
	    if id in non_seeds: # not in the initial association file but in the candidates
		if score >= score_threshold:
		    nFP += 1 
		else:
		    nTN += 1 
    return (nTP, nFP, nFN, nTN)

def calculate_performance_metric_counts_using_random_negatives(dictNodeResult, setNodeTest, non_seeds, score_threshold, n_random_negative_folds = None, replicable=123):
    (nTP, nFP, nFN, nTN) = (0.0, 0.0, 0.0, 0.0)
    for id, score in dictNodeResult.iteritems(): # if candidates based - for each candidate
        if id in setNodeTest: # in the initial association file
            if score >= score_threshold:
                nTP += 1
            else:
                nFN += 1

    if n_random_negative_folds == 0:
	for id, score in dictNodeResult.iteritems():
	    if id in non_seeds:
		if score >= score_threshold:
		    nFP += 1
		else:
		    nTN += 1
    else:
	n_actual_folds = 0
	for sample in generate_samples_from_list_without_replacement(non_seeds, len(setNodeTest), n_random_negative_folds, replicable = replicable):
	    setNegative = set(sample)
	    n_actual_folds += 1
	    for id, score in dictNodeResult.iteritems():
		if id in setNegative:
		    if score >= score_threshold:
			nFP += 1
		    else:
			nTN += 1
	nFP /= n_actual_folds
	nTN /= n_actual_folds
    return (nTP, nFP, nFN, nTN)


def calculate_performance_metric_counts_using_files(file_result, file_seed_test_scores, file_node_scores, score_threshold, n_random_negative_folds = None, default_score = 0, replicable=123, candidates_file = None):
    """
	Calculate and return TP, FP, TN, FN (FP and TN based on random selected non-seed nodes)
	file_result: output scores file
	file_seed_test_scores: seed nodes separated for test
	file_node_scores: initial node/seed scores
	score_threshold: threshold for considering data as P, N
	n_random_negative_folds: Number of negative data selection folds (each fold contain same number of test nodes) to be averaged for FP and TN calculation
				 FP and TN are not necesserily integers when n_random_negative_folds > 1
				 If None calculated to cover as much as non-seed scores as possible
				 If 0 all negative data is used
    """
    dictNodeResult, setNodeTest, non_seeds = get_values_from_files_for_performance_metric_counts(file_result, file_seed_test_scores, file_node_scores, default_score, candidates_file)
    return calculate_performance_metric_counts_using_random_negatives(dictNodeResult, setNodeTest, non_seeds, score_threshold, n_random_negative_folds, replicable)


def get_values_from_files_for_performance_metric_counts(file_result, file_seed_test_scores, file_node_scores, default_score, candidates_file = None):
    setNodeResult, setDummy, dictNodeResult, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_result, store_edge_type = False)
    setNodeTest, setDummy, dictNodeTest, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_seed_test_scores, store_edge_type = False)
    non_seeds = get_non_seeds_from_node_scores_file(file_node_scores, default_score = default_score)
    if candidates_file is not None:
	setNode, setDummy, dictNode, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = candidates_file, store_edge_type = False)
	#dictNodeResult = dictNodeResult.fromkeys(setNode)
	nodes = dictNodeResult.keys()
	for node in nodes:
	    if node not in setNode:
		del dictNodeResult[node]
    return dictNodeResult, setNodeTest, non_seeds


def calculatePerformance(nTP, nFP, nFN, nTN):
    try:
	acc = (nTP + nTN) / (nTP + nFP + nTN + nFN)
    except ZeroDivisionError:
	acc = None
    try:
        sens = nTP / (nTP + nFN)
    except:
        sens = None
    try:
        spec = nTN / (nTN + nFP)
    except:
        spec = None
    try:
        ppv = nTP / (nTP + nFP)
    except:
        ppv = None

    #if spec is not None:
    #    return (sens, (1-spec))
    #else:
    #    return (sens, None)

    return (acc, sens, spec, ppv)


def get_top_scoring_nodes_at_given_cutoff(node_to_score, cutoff):
    """
	Returns highest scoring nodes at given cutoff (if ends with % taken as percentage, otherwise taken as score cutoff)
    """

    result_scores = node_to_score.items()
    result_scores.sort(lambda x,y: cmp(y[1], x[1]))

    ids = [] # the order is important

    if str(cutoff).endswith("%"):
	percentage = float(cutoff.rstrip("%"))
	i = 0
	n = len(result_scores)*percentage/100
	last_score = result_scores[0][1]
	for id, score in result_scores:
	    i+=1
	    if i>n:
		if last_score != score:
		    break
	    ids.append(id)
	    last_score = score
    else:
	cutoff = float(cutoff)
	i = 0
	for id, score in result_scores:
	    i += 1
	    if score >= cutoff:
		ids.append(id)
	n = i

    return ids, n, i


def calculate_seed_coverage_at_given_cutoff(file_result, file_seed_scores, cutoff, default_score=0):
    """
	Calculates number of seed nodes included in the given cutoff of high ranking nodes in result file 
	Returns a tuple containing number of seed nodes and all nodes at that cutoff (if ends with % taken as percentage, otherwise taken as score cutoff)
    """

    n_seed = 0

    dummy, dummy, node_to_score_initial, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_seed_scores, store_edge_type = False)
    # getting seeds (using initial non-seed score assumption)
    seeds = set([ id for id, score in node_to_score_initial.iteritems() if score > default_score ])

    dummy, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_result, store_edge_type = False)
    ids, n, i = get_top_scoring_nodes_at_given_cutoff(node_to_score, cutoff)
    for id in ids:
	if id in seeds:
    	    n_seed += 1
    return n_seed, len(seeds), n, i
    
    #setDummy, setDummy, dictNodeResult, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_result, store_edge_type = False)
    #setDummy, setDummy, dictNode, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_seed_scores, store_edge_type = False)
    
    #result_scores = dictNodeResult.items()
    #result_scores.sort(lambda x,y: cmp(y[1], x[1]))

    #i = 0
    #n = len(result_scores)*percentage/100
    #n_seed = 0

    #last_score = result_scores[0][1]
    #for id, score in result_scores:
    #	i+=1
    #	if i>n:
    #	    if last_score != score:
    #		#print "i,n:", i, n
    #		#print "id,score,last:", id, score, last_score
    #		break
    #	# checking whether node was a seed (using initial non-seed score assumption)
    #	if dictNode.has_key(id) and dictNode[id] > default_score:
    #	    n_seed += 1
    #	last_score = score
    
    #return n_seed, len(dictNode), n, i


