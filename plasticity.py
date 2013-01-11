import os

DATA_DIR = "../data/"

omim_phenotypes = ["alzheimer", "breast cancer", "diabetes", "anemia", "myopathy", "obesity", "parkinson disease", "prostate cancer", "hypertension", "leukemia", "lung cancer", "ataxia", "epilepsy", "schizophrenia", "cardiomyopathy", "cataract", "lymphoma", "mental retardation", "systemic lupus erythematosus"] # insulin, neuropathy, asthma, spastic paraplegia
omim_phenotypes = [ "omim_" + "_".join(p.split()) for p in omim_phenotypes ]

def main():
    #output_ns_related_modules()
    case_study_pruned_networks()
    #case_study_neighborhood()
    #get_go_function_counts()
    #get_similarity_of_go_terms_for_diseases()
    ##get_similarity_of_go_terms_for_omim_diseases()
    #get_number_of_seed_connecting_edges()
    #get_number_of_seed_connecting_paths()
    #check_modularization()
    return

def get_similarity_of_go_terms_for_diseases():
    
    base_dir = DATA_DIR + "/module/biana_no_tap-omim/"

    phenotype_to_functions = get_go_function_counts()
    phenotypes = phenotype_to_functions.keys()
    #for disease, values in disease_to_go_terms.iteritems():
    #	seed_terms, module_terms, all_terms = values

    # Functional similarity matrix
    f = open(base_dir + "functional_similarity.dat", 'w')
    f2 = open(base_dir + "functional_similarity_matrix.dat", 'w')
    f3 = open(base_dir + "functional_similarity_seed_matrix.dat", 'w')
    f2.write("%s\n" % " ".join(phenotypes))
    f3.write("%s\n" % " ".join(phenotypes))
    f.write("phenotype1 phenotype2 n_go_1 n_go_2 n_go_intersection n_go_union jaccard\n")
    #print phenotype_to_functions
    for i, pheno1 in enumerate(phenotypes):
	f2.write("%s " % pheno1)
	f3.write("%s " % pheno1)
	values = []
	values2 = []
	for j, pheno2 in enumerate(phenotypes):
	    s1, m1, g1 = phenotype_to_functions[pheno1]
	    s2, m2, g2 = phenotype_to_functions[pheno2]
	    if len(g1|g2) != 0:
		jaccard = float(len(g1&g2))/len(g1|g2)
	    else:
		jaccard = 0
	    if len(s1|s2) != 0:
		jaccard2 = float(len(s1&s2))/len(s1|s2)
	    else:
		jaccard2 = 0
	    if i < j:
		f.write("%s %s %d %d %d %d %f\n" % (pheno1, pheno2, len(g1), len(g2), len(g1&g2), len(g1|g2), jaccard))
	    values.append(jaccard)
	    values2.append(jaccard2)
	f2.write("%s\n" % " ".join(map(str, values)))
	f3.write("%s\n" % " ".join(map(str, values2)))
    f.close()
    f2.close()
    f3.close()

    from toolbox import OboParser
    g=OboParser.getOboGraph("/home/emre/arastirma/celldiff/data/GO/gene_ontology.1_2.obo")

    # Phenotype vs function matrix
    f = open(base_dir + "phenotype_vs_functions.dat", 'w')
    f.write("%s\n" % " ".join(phenotypes))
    all_go_ids = reduce(lambda x,y: x|y, zip(*phenotype_to_functions.values())[2])
    for go_id in all_go_ids:
	f.write("%s " % "_".join(g.node[go_id]['n'].split(" "))) 
	values = []
	for pheno in phenotypes:
	    if go_id in phenotype_to_functions[pheno][2]:
		values.append(1)
	    else:
		values.append(0)
	f.write("%s\n" % " ".join(map(str, values)))
    f.close()

    f = open(base_dir + "phenotype_vs_functions_seed.dat", 'w')
    f.write("%s\n" % " ".join(phenotypes))
    all_seed_go_ids = reduce(lambda x,y: x|y, zip(*phenotype_to_functions.values())[0])
    for go_id in all_seed_go_ids:
	f.write("%s " % "_".join(g.node[go_id]['n'].split(" ")))
	values = []
	for pheno in phenotypes:
	    if go_id in phenotype_to_functions[pheno][0]:
		values.append(1)
	    else:
		values.append(0)
	f.write("%s\n" % " ".join(map(str, values)))
    f.close()
    return

def check_modularization():
    from toolbox import mcl_utilities
    base_dir = DATA_DIR + "input/biana_no_tap_no_reliability/"
    network_file = base_dir + "edge_scores.sif"
    module_file = base_dir + "modules.txt"
    g = mcl_utilities.create_network_from_sif_file(network_file, use_edge_data=True)
    if not os.path.exists(module_file):
	modules = mcl_utilities.get_modules_of_graph(g, "mcl", output_file=module_file, inflation=1.7)
    else:
	modules = mcl_utilities.get_modules_from_file(module_file)                  
    print map(lambda x: len(x), modules)
    return

def get_number_of_seed_connecting_paths():
    from toolbox import network_utilities as gu
    #base_dir = DATA_DIR + "input/biana_no_tap_no_reliability/"
    base_dir = DATA_DIR + "input/goh/"
    network_file = base_dir + "edge_scores.sif"
    output_file = base_dir + "path_counts.txt"
    g = gu.create_network_from_sif_file(network_file, use_edge_data=True)
    f = open(output_file, 'w') 
    for phenotype in omim_phenotypes:
	print phenotype
	node_file = base_dir + phenotype + "/seed_scores.sif"
	seeds = [ line.strip().split()[0] for line in open(node_file) ]
	n = float(len(seeds))
	count = 0
	path_length = 0
	for i, seed1 in enumerate(seeds):
	    for j, seed2 in enumerate(seeds):
		if i<j:
		    #count += len(find_all_paths(g, seed1, seed2, []))
		    for path in all_shortest_paths(g, seed1, seed2):
			count += 1
			path_length += len(path) - 1
	#print count, count/n, float(path_length)/count
	f.write("%s %d %f %f\n" % (phenotype, count, count/n, float(path_length)/count))
    f.close()
    return

# From networkx 1.7
def all_shortest_paths(G, source, target, weight=None):
    """Compute all shortest paths in the graph.

Parameters
----------
G : NetworkX graph

source : node
Starting node for path.

target : node
Ending node for path.

weight : None or string, optional (default = None)
If None, every edge has weight/distance/cost 1.
If a string, use this edge attribute as the edge weight.
Any edge attribute not present defaults to 1.

Returns
-------
paths: generator of lists
A generator of all paths between source and target.

Examples
--------
>>> G=nx.Graph()
>>> G.add_path([0,1,2])
>>> G.add_path([0,10,2])
>>> print([p for p in nx.all_shortest_paths(G,source=0,target=2)])
[[0, 1, 2], [0, 10, 2]]

Notes
-----
There may be many shortest paths between the source and target.

See Also
--------
shortest_path()
single_source_shortest_path()
all_pairs_shortest_path()
"""
    import networkx as nx
    if weight is not None:
        pred,dist = nx.dijkstra_predecessor_and_distance(G,source,weight=weight)
    else:
        pred = nx.predecessor(G,source)
    if target not in pred:
        raise nx.NetworkXException()
    stack = [[target,0]]
    top = 0
    while top >= 0:
        node,i = stack[top]
        if node == source:
          yield [p for p,n in reversed(stack[:top+1])]
        if len(pred[node]) > i:
            top += 1
            if top == len(stack):
                stack.append([pred[node][i],0])
            else:
                stack[top] = [pred[node][i],0]
        else:
            stack[top-1][1] += 1
            top -= 1
    return


def get_number_of_seed_connecting_edges():
    from toolbox import network_utilities as gu
    #base_dir = DATA_DIR + "input/biana_no_tap_no_reliability/"
    base_dir = DATA_DIR + "input/goh/"
    network_file = base_dir + "edge_scores.sif"
    output_file = base_dir + "edge_counts.txt"
    g = gu.create_network_from_sif_file(network_file, use_edge_data=True)
    f = open(output_file, 'w')
    n = float(g.number_of_edges())
    for phenotype in omim_phenotypes:
	node_file = base_dir + phenotype + "/seed_scores.sif"
	seeds = [ line.strip().split()[0] for line in open(node_file) ]
	n_seed = 0
	n_seed_nonseed = 0
	n_nonseed = 0
	for u, v in g.edges():
	    if u in seeds and v in seeds:
		n_seed += 1
	    elif u in seeds or v in seeds:
		n_seed_nonseed += 1
	    else:
		n_nonseed += 1
	#n_seed /= n
	#n_seed_nonseed /= n
	#n_nonseed /= n
	f.write("%s %d %d %d\n" % (phenotype, n_seed, n_seed_nonseed, n_nonseed))
    f.close()
    return


def get_similarity_of_go_terms_for_omim_diseases():
    base_dir = DATA_DIR + "module/biana_no_tap-omim/"
    f = open(base_dir + "modules.txt")
    phenotype_to_functions = {}
    phenotype_to_genes = {}
    flag = False
    phenotypes = []
    go_id_to_go_term = {}
    for line in f:
	if line.endswith(" ns\n"):
	    flag = True
	    pheno = line.strip().split(" ")[-2]
	    phenotype_to_functions[pheno] = set()
	    phenotype_to_genes[pheno] = set()
	    phenotypes.append(pheno)
	    continue
	if flag:
	    if line.startswith("biana_no_tap"):
		flag = False
		continue
	    words = line.strip().split("\t")
	    if len(words) > 5 and words[5].startswith("GO:"):
		# Get go terms
		go_id_to_go_term[words[5]] = words[6]
		phenotype_to_functions[pheno].add(words[5])
	    elif line[0] != "#":
		# Get genes
		words = line.strip().split(", ")
		for word in words:
		    phenotype_to_genes[pheno].add(word)
    f.close()
    #print len(phenotype_to_functions["omim_alzheimer"]), len(phenotype_to_genes["omim_alzheimer"])

    # Functional similarity matrix
    f = open(base_dir + "functional_similarity.dat", 'w')
    f2 = open(base_dir + "functional_similarity_matrix.dat", 'w')
    f2.write("%s\n" % " ".join(phenotypes))
    f.write("phenotype1 phenotype2 n_go_1 n_go_2 n_go_intersection n_go_union jaccard\n")
    #print phenotype_to_functions
    for i, pheno1 in enumerate(phenotypes):
	f2.write("%s " % pheno1)
	values = []
	for j, pheno2 in enumerate(phenotypes):
	    g1 = phenotype_to_functions[pheno1]
	    g2 = phenotype_to_functions[pheno2]
	    if len(g1|g2) != 0:
		jaccard = float(len(g1&g2))/len(g1|g2)
	    else:
		jaccard = 0
	    if i < j:
		f.write("%s %s %d %d %d %d %f\n" % (pheno1, pheno2, len(g1), len(g2), len(g1&g2), len(g1|g2), jaccard))
	    values.append(jaccard)
	f2.write("%s\n" % " ".join(map(str, values)))
    f.close()
    f2.close()

    # Gene similarity matrix
    f = open(base_dir + "gene_similarity_matrix.dat", 'w')
    f.write("%s\n" % " ".join(phenotypes))
    for i, pheno1 in enumerate(phenotypes):
	f.write("%s " % pheno1)
	values = []
	for j, pheno2 in enumerate(phenotypes):
	    g1 = phenotype_to_genes[pheno1]
	    g2 = phenotype_to_genes[pheno2]
	    if len(g1|g2) != 0:
		jaccard = float(len(g1&g2))/len(g1|g2)
	    else:
		jaccard = 0
	    values.append(jaccard)
	f.write("%s\n" % " ".join(map(str, values)))
    f.close()

    # Phenotype vs function matrix
    f = open(base_dir + "phenotype_vs_functions.dat", 'w')
    f.write("%s\n" % " ".join(phenotypes))
    all_go_ids = reduce(lambda x,y: x|y, phenotype_to_functions.values())
    for go_id in all_go_ids:
	f.write("%s " % "_".join(go_id_to_go_term[go_id].split(" ")))
	values = []
	for pheno in phenotypes:
	    if go_id in phenotype_to_functions[pheno]:
		values.append(1)
	    else:
		values.append(0)
	f.write("%s\n" % " ".join(map(str, values)))
    f.close()

    # Phenotype vs gene matrix
    f = open(base_dir + "phenotype_vs_genes.dat", 'w')
    f.write("%s\n" % " ".join(phenotypes))
    all_genes = reduce(lambda x,y: x|y, phenotype_to_genes.values())
    for gene in all_genes:
	f.write("%s " % gene)
	values = []
	for pheno in phenotypes:
	    if gene in phenotype_to_genes[pheno]:
		values.append(1)
	    else:
		values.append(0)
	f.write("%s\n" % " ".join(map(str, values)))
    f.close()

    #duplicated_genes = set(line.strip() for line in open("/home/emre/arastirma/data/disease/itan2010_gene_duplication/duplicated_genes.txt"))
    duplicated_genes = set(line.strip() for line in open("/home/emre/arastirma/data/disease/itan2010_gene_duplication/duplicated_genes_cheung.txt"))
    for pheno in phenotypes:
	print pheno, len(phenotype_to_genes[pheno]), len(phenotype_to_genes[pheno] & duplicated_genes)
    return 


def case_study_pruned_networks():
    from toolbox import network_utilities as gu
    from toolbox import functional_enrichment
    from toolbox import mcl_utilities as mcl
    from scipy.stats import hypergeom

    network_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_no_reliability/edge_scores.sif"
    user_entity_id_mapping_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_no_reliability/node_mapping.tsv.genesymbol.single"
    seeds_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_no_reliability/omim_breast_cancer/seed_scores.sif"
    auc_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/breast_cancer_pruned_p80.txt"
    module_file = DATA_DIR + "module/biana_no_tap-omim/mcl/modules.txt"

    # Get node mapping
    ueid_to_gene = get_ues_gene_mapping(user_entity_id_mapping_file)
    # Get seeds
    seeds = set([line.strip().split()[0] for line in open(seeds_file)])
    # Get neighborhood in the original network
    g_org = gu.create_network_from_sif_file(network_file, use_edge_data=False)
    g_neighborhood = gu.get_neighborhood_subgraph(g_org, seeds)
    #neighborhood_edges = set(g_neighborhood.edges()) # edge node order may be different for the same edge

    # Get indices of min/max networks
    aucs = []
    for line in open(auc_file):
	aucs.append(float(line.split()[1]))
    indices = zip(*sorted([ (auc, i) for i, auc in enumerate(aucs) ]))[1]
    print indices[:3], indices[-3:] # real index in file name is one higher

    # Get max neighborhood network
    g_maxs = [ ] 
    for i in indices[-2:]:
	network_file_pruned = DATA_DIR + "human_interactome_biana/pruned/omim_breast_cancer/80/sampled_graph.sif.%d" % (i+1)
	g = gu.create_network_from_sif_file(network_file_pruned, use_edge_data=False)
	g_maxs.append(g) 
	#g_maxs.append(g.subgraph(g_neighborhood.nodes()))

    # Get min neighborhood network
    g_mins = [ ] 
    for i in indices[:2]:
	network_file_pruned = DATA_DIR + "human_interactome_biana/pruned/omim_breast_cancer/80/sampled_graph.sif.%d" % (i+1)
	g = gu.create_network_from_sif_file(network_file_pruned, use_edge_data=False)
	g_mins.append(g)
	#g_mins.append(g.subgraph(g_neighborhood.nodes()))

    # Get common edges in min/max networks
    g_max = reduce(lambda x,y: gu.networkx.intersection(x, y), g_maxs)
    g_min = reduce(lambda x,y: gu.networkx.intersection(x, y), g_mins)

    # Get differential edges
    g_diff = gu.networkx.difference(g_max, g_min)
    print len(g_max.edges()), len(g_min.edges()), len(g_diff.edges())
    nodes = set()
    for node in g_diff.nodes():
	if node in ueid_to_gene:
	    nodes.add(node)
    g_sub = g_diff.subgraph(nodes)

    #!
    weak_edges = set()
    strong_edges = g_sub.edges()
    output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/breast_cancer_pruned_p80_diff.dot"
    gu.create_dot_network_file(g_sub, output_file, seeds, ueid_to_gene, weak_edges=weak_edges, draw_type="all")
    os.system("fdp -Tgif -O %s" % output_file) 
    return

    # Get seed GOs to check their coverage in top connected component
    phenotype_to_functions = get_go_function_counts() 
    seed_terms = phenotype_to_functions["omim_breast_cancer"][0]

    # Get all functions enriched in the network
    output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/network_genes.txt"
    f = open(output_file, 'w')
    [ f.write("%s\n" % gene) for gene in set(ueid_to_gene.values()) ]
    f.close()
    #functional_enrichment.check_functional_enrichment_of_human_gene_symbols(output_file, output_file+".funcassoc")
    network_go_terms = functional_enrichment.get_functional_enrichment(output_file + ".funcassoc", remove_parents=False, only_biological_processes=True)
    print 23928, len(network_go_terms)
    
    # Check the functions enriched in the largest connected component of each module
    for i, module in enumerate(gu.get_connected_components(g_sub, return_as_graph_list=False)):
	if len(module) < 10:
	    continue
	output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/M%d_.txt" % i
	# Draw diff network component
	weak_edges = set()
	g_sub_sub = g_sub.subgraph(module)
	gu.create_dot_network_file(g_sub_sub, output_file+".dot", seeds, ueid_to_gene, weak_edges=weak_edges, draw_type="all")
	os.system("fdp -Tgif -O %s" % output_file+".dot")
	# Get functions
	f = open(output_file, 'w')
	for node in module:
	    f.write("%s\n" % ueid_to_gene[node])
	f.close()
	#functional_enrichment.check_functional_enrichment_of_human_gene_symbols(output_file, output_file+".funcassoc")
	go_terms = functional_enrichment.get_functional_enrichment(output_file + ".funcassoc", remove_parents=False, only_biological_processes=True)
	print len(seed_terms), len(go_terms), len(seed_terms & go_terms), len(seed_terms & go_terms) / float(len(seed_terms))
	print "p_value:", sum(hypergeom.pmf(range(len(seed_terms & go_terms),len(go_terms)+1), len(network_go_terms), len(seed_terms), len(go_terms)))

    weak_edges = set()
    strong_edges = g_sub.edges()

    # Draw diff neighborhood network
    output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/breast_cancer_pruned_p80_diff.dot"
    gu.create_dot_network_file(g_sub, output_file, seeds, ueid_to_gene, weak_edges=weak_edges, draw_type="all")
    os.system("fdp -Tgif -O %s" % output_file) 
    return

    # Get modules of high scoring network
    modules = mcl.get_modules_from_file(module_file)
    
    # Create nested Cytoscape network file
    output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/breast_cancer_pruned_p80_modularized_diff.nnf"
    f = open(output_file, 'w')
    network_name = "breast_cancer_pruned_p80_"
    module_sets = []
    included_nodes = set()
    included_edges = set()
    # Output modules
    for i, module in enumerate(modules):
	m = set(module) & nodes
	if len(m) > 0:
	    module_sets.append(m)
	    f.write("%s M%d_\n" % (network_name, i))
	    for u, v in g_sub.edges(m):
		if u == v:
		    continue
		if u in m and v in m:
		    w = "pp"
		    if (u,v) in weak_edges or (v,u) in weak_edges:
			w = "weak"
		    elif (u,v) in strong_edges or (v,u) in strong_edges:
			w = "strong"
		    included_nodes.add(u)
		    included_nodes.add(v)
		    included_edges.add((u,v))
		    included_edges.add((v,u))
		    u = ueid_to_gene[u]
		    v = ueid_to_gene[v]
		    f.write("M%d_ %s %s %s\n" % (i, u, w, v))
	    for u in m:
		if u not in included_nodes:
		    included_nodes.add(u)
		    u = ueid_to_gene[u]
		    f.write("M%d_ %s\n" % (i, u))
    # Connect modules
    for i, module1 in enumerate(module_sets):
	for j, module2 in enumerate(module_sets):
	    if i<j:
		connected_weak = False
		connected_strong = False
		for u in module1:
		    for v in module2:
			if (u,v) in strong_edges:
			    connected_strong = True
			    break
			if (u,v) in weak_edges:
			    connected_weak = True
		if connected_strong:
		    f.write("%s M%d_ %s M%d_\n" % (network_name, i, "strong", j))
		elif connected_weak:
		    f.write("%s M%d_ %s M%d_\n" % (network_name, i, "weak", j))
    # Output the rest
    for u,v in g_sub.edges():
	if u == v:
	    continue
	if (u,v) in included_edges:
	    continue
	included_nodes.add(u)
	included_nodes.add(v)
	w = "pp"
	if (u,v) in weak_edges or (v,u) in weak_edges:
	    w = "weak"
	elif (u,v) in strong_edges or (v,u) in strong_edges:
	    w = "strong"
	u = ueid_to_gene[u]
	v = ueid_to_gene[v]
	f.write("%s %s %s %s\n" % (network_name, u,w,v))
    for node in nodes - included_nodes:
	f.write("%s %s\n" % (network_name, ueid_to_gene[node]))
    f.close()

    # Check the functions enriched in the largest connected component of each module
    for i, module in enumerate(module_sets):
	output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/M%d_.txt" % i
	g_sub_sub = g_sub.subgraph(module) 
	module = gu.get_connected_components(g_sub_sub, return_as_graph_list=False)[0] 
	f = open(output_file, 'w')
	for node in module:
	    f.write("%s\n" % ueid_to_gene[node])
	f.close()
	functional_enrichment.check_functional_enrichment_of_human_gene_symbols(output_file, output_file+".funcassoc")
    return

def case_study_pruned_networks_old():
    """
    cat /sbi/users/emre/data/netzcore/from_gaudi_2011/output_runs_on_random/biana_no_tap_no_reliability_pruned_p50_*/omim_breast_cancer/ns/r3i2/auc.txt > arastirma/netzcore/data/summary_runs_on_random/breast_cancer_pruned_p80.txt
    vi %s/"//g
    d<-read.table("breast_cancer_pruned_p50.txt")
    e<-d$V2
    f<-(e-mean(e))/sd(e)
    which(e %in% sort(e)[98:100])
    > 22 54 58
    which(e %in% sort(e)[1:3])
    > 7 39 55 99
    """
    from toolbox import network_utilities as gu
    from toolbox import functional_enrichment

    network_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_no_reliability/edge_scores.sif"
    user_entity_id_mapping_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_no_reliability/node_mapping.tsv.genesymbol.single"
    seeds_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_no_reliability/omim_breast_cancer/seed_scores.sif"
    network_file_pruned = DATA_DIR + "human_interactome_biana/pruned/omim_breast_cancer/80/sampled_graph.sif.58"
    network_file_pruned2 = DATA_DIR + "human_interactome_biana/pruned/omim_breast_cancer/80/sampled_graph.sif.7"
    module_file = DATA_DIR + "module/biana_no_tap-omim/mcl/modules.txt"
    #network_file_permuted = DATA_DIR + "human_interactome_biana/permuted/50/sampled_graph.sif.46"
    ueid_to_gene = get_ues_gene_mapping(user_entity_id_mapping_file)
    seeds = set([line.strip().split()[0] for line in open(seeds_file)])
    g = gu.create_network_from_sif_file(network_file, use_edge_data=False)
    g_neighborhood = gu.get_neighborhood_subgraph(g, seeds)
    neighborhood_edges = set(g_neighborhood.edges()) # edge node order may be different for the same edge
    #g_sub_pruned = gu.get_neighborhood_subgraph(g_pruned, seeds)
    #print len(g_sub.nodes()), len(g_sub.edges())
    #print len(g_sub_pruned.nodes()), len(g_sub_pruned.edges())
    g_pruned = gu.create_network_from_sif_file(network_file_pruned, use_edge_data=False)
    g_pruned2 = gu.create_network_from_sif_file(network_file_pruned2, use_edge_data=False)
    #weak_edges = set(g.edges()) - set(g_pruned.edges())
    strong_edges = set(g_pruned.edges())
    #strong_edges = set(g.edges()) - set(g_pruned2.edges())
    weak_edges = set(g_pruned2.edges())
    common_edges = weak_edges & strong_edges
    weak_edges -= common_edges
    strong_edges -= common_edges
    weak_edges &= neighborhood_edges
    strong_edges &= neighborhood_edges
    #print len(weak_edges), len(strong_edges), len(common_edges)
    g_sub = gu.create_graph()
    #g_sub.add_edges_from(weak_edges | strong_edges)

    #strong_edges = weak_edges # To check differential network from the other side (edges in min but not in max)

    g_sub.add_edges_from(strong_edges) 
    weak_edges = set() 

    # Run scoring on pruned networks
    if False:
	from toolbox import guild_utilities
	data_dir = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/" 
	executable_path = "scoreNetwork/scoreN"
	for network_type in ("pruned_max", "pruned_min"):
	    if network_type == "pruned_max":
		network_file = network_file_pruned
	    elif network_type == "pruned_min":
		network_file = network_file_pruned2
	    else:
		raise ValueError("Unknown network type!")
	    scoring_folder = data_dir + network_type + os.sep
	    # Create input files for scoring
	    guild_utilities.prepare_scoring(network_file, seeds_file, scoring_folder, non_seed_score=0.01, seed_score=1.0, edge_score=1.0, n_sample=1, delim=" ", name=None)
	    # Run GUILD and create output files, the case for Netcombo
	    guild_utilities.run_scoring(scoring_folder, executable_path, scoring_type="netscore", parameters={"n_iteration":2, "n_repetition":3}, qname=None, name=None, calculate_pvalue=False)

    # Get functions of high scoring portions in 3 networks 
    network_types = ("original", "pruned_max", "pruned_min")
    if False:
	import analyze_results
	association_scores_file_identifier_type = "genesymbol"
	node_mapping_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_no_reliability/node_mapping.tsv"
	node_mapping_file += "."+association_scores_file_identifier_type
	node_scores_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_no_reliability/omim_breast_cancer/node_scores.sif"
	for network_type in network_types:
	    if network_type == "original":
		output_scores_file = DATA_DIR + "output_runs_for_draft/biana_no_tap_no_reliability/omim_breast_cancer/ns/r3i2/node_scores.sif"
	    else: # network_type in ("pruned_max", "pruned_min"):
		output_scores_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/%s/output_scores.sif.netscore" % network_type
	    enrichment_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/%s.txt" % network_type
	    file_enrichment = open(enrichment_file, 'w')
	    analyze_results.check_functional_enrichment_at_given_cutoff(output_scores_file, node_scores_file, node_mapping_file, "5%", association_scores_file_identifier_type, file_enrichment.write, 0.01, exclude_seeds=False, specie = "Homo sapiens")
	    file_enrichment.close()

    # Get functions of the top scroing portions for all diseases
    if False: 
	from toolbox import OboParser

    	phenotype_to_functions = get_go_function_counts() 
    	seed_terms = phenotype_to_functions["omim_breast_cancer"][0]

	go = OboParser.getOboGraph("/home/emre/arastirma/celldiff/data/GO/gene_ontology.1_2.obo")
	all_terms = set()
	common_terms = None
	network_to_terms = {}
	#network_types = network_types[:2]
	for network_type in network_types:
	    enrichment_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/%s.txt" % network_type
	    go_terms = functional_enrichment.get_functional_enrichment(enrichment_file, remove_parents=False, only_biological_processes=True, only_slim=False)
	    network_to_terms[network_type] = set() | go_terms
	    print network_type, len(go_terms)
	    all_terms |= go_terms
	    if common_terms is None:
		common_terms = go_terms
	    else:
		common_terms &= go_terms
	all_terms |= seed_terms
	common_terms &= seed_terms
	#seed_terms = functional_enrichment.remove_parent_terms(seed_terms, go) 
	#all_terms = functional_enrichment.remove_parent_terms(all_terms, go)
	#common_terms = functional_enrichment.remove_parent_terms(common_terms, go)
	print len(all_terms), len(common_terms)

	for network_type in network_types:
	    for network_type2 in network_types:
		if network_type == network_type2:
		    continue
		print network_type, network_type2, len(all_terms & network_to_terms[network_type] & network_to_terms[network_type2])

	f = open(DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/functional_comparison.dat", 'w')
	f.write("seed terms\t%s\n" % "\t".join(network_types))
	for go_term in all_terms:
	    values = []
	    if go_term in seed_terms: 
		val = 1
	    else:
		val = 0
	    values.append(val)
	    for network_type in network_types:
		val = 0
		if go_term in network_to_terms[network_type]:
		    val= 1
		values.append(val)
	    f.write("%s\t%s\n" % (go.node[go_term]['n'], "\t".join(map(str, values))))
	f.close()

    if False:
	n = float(len(seeds))
	for graph in (g, g_pruned, g_pruned2):
	    count = 0
	    path_length = 0
	    n_path = 0
	    n_pair = 0.0
	    for i, seed1 in enumerate(seeds):
		for j, seed2 in enumerate(seeds):
		    if i<j:
			#count += len(find_all_paths(g, seed1, seed2, []))
			try:
			    paths = all_shortest_paths(graph, seed1, seed2)
			    for path in paths:
				count += 1
				path_length += len(path) - 1
			    n_pair += 1
			    n_path += len(path) - 1
			except:
			    continue
	    print n, n_pair, n_path/n_pair, count, count/n, count/n_pair, path_length/float(count)
	return

    # Check seed interaction counts on pruned max
    if False:
	output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/breast_cancer_pruned_p80_seeds.txt"
	f = open(output_file, 'w')
	for seed in seeds:
	    if seed not in ueid_to_gene:
		print seed
		continue
	    f.write("%s\n" % ueid_to_gene[seed])
	f.close()

	for network_type, graph in zip(network_types, [g, g_pruned, g_pruned2]):
	    output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/%s_seed_interaction_counts.txt" % network_type
	    g_neighborhood = gu.get_neighborhood_subgraph(graph, seeds)
	    f = open(output_file, 'w')
	    nodes = set()
	    for node in g_neighborhood.nodes():
		if node in ueid_to_gene:
		    nodes.add(node)
		    if node in seeds:
			f.write("%s\t%d\n" % (ueid_to_gene[node], g_neighborhood.degree(node)))
	    f.close()
	    #g_neighborhood = g_neighborhood.subgraph(nodes)    
	    #output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/breast_cancer_pruned_p0.dot"
	    #gu.create_dot_network_file(g_neighborhood, output_file, seeds, ueid_to_gene, draw_type="all")
	    #os.system("twopi -Tgif -O %s" % output_file)

    # Check modules in pruned max
    if False:
	nodes = set()
	for node in g_sub.nodes():
	    if node in ueid_to_gene:
		nodes.add(node)
	g_sub = g_sub.subgraph(nodes)

	output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/breast_cancer_pruned_p80_strong.dot"
	#gu.create_dot_network_file(g_sub, output_file, seeds, ueid_to_gene, weaks=weaks, draw_type="weak")
	#gu.create_dot_network_file(g_sub_pruned, output_file, seeds, ueid_to_gene, weaks=weaks, draw_type="weak")
	gu.create_dot_network_file(g_sub, output_file, seeds, ueid_to_gene, weak_edges=weak_edges, draw_type="all")
	os.system("twopi -Tgif -O %s" % output_file)
	
	from toolbox import mcl_utilities as mcl
	modules = mcl.get_modules_from_file(module_file)
	
	output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/breast_cancer_pruned_p80_modularized_strong.nnf"
	f = open(output_file, 'w')

	network_name = "breast_cancer_pruned_p80_"
	module_sets = []
	included_nodes = set()
	included_edges = set()
	# Output modules
	for i, module in enumerate(modules):
	    m = set(module) & nodes
	    if len(m) > 0:
		module_sets.append(m)
		f.write("%s M%d_\n" % (network_name, i))
		for u, v in g_sub.edges(m):
		    if u == v:
			continue
		    if u in m and v in m:
			w = "pp"
			if (u,v) in weak_edges or (v,u) in weak_edges:
			    w = "weak"
			elif (u,v) in strong_edges or (v,u) in strong_edges:
			    w = "strong"
			included_nodes.add(u)
			included_nodes.add(v)
			included_edges.add((u,v))
			included_edges.add((v,u))
			u = ueid_to_gene[u]
			v = ueid_to_gene[v]
			f.write("M%d_ %s %s %s\n" % (i, u, w, v))
		for u in m:
		    if u not in included_nodes:
			included_nodes.add(u)
			u = ueid_to_gene[u]
			f.write("M%d_ %s\n" % (i, u))
	# Connect modules
	for i, module1 in enumerate(module_sets):
	    for j, module2 in enumerate(module_sets):
		if i<j:
		    connected_weak = False
		    connected_strong = False
		    for u in module1:
			for v in module2:
			    if (u,v) in strong_edges:
				connected_strong = True
				break
			    if (u,v) in weak_edges:
				connected_weak = True
		    if connected_strong:
			f.write("%s M%d_ %s M%d_\n" % (network_name, i, "strong", j))
		    elif connected_weak:
			f.write("%s M%d_ %s M%d_\n" % (network_name, i, "weak", j))
	# Output the rest
	for u,v in g_sub.edges():
	    if u == v:
		continue
	    if (u,v) in included_edges:
		continue
	    included_nodes.add(u)
	    included_nodes.add(v)
	    w = "pp"
	    if (u,v) in weak_edges or (v,u) in weak_edges:
		w = "weak"
	    elif (u,v) in strong_edges or (v,u) in strong_edges:
		w = "strong"
	    u = ueid_to_gene[u]
	    v = ueid_to_gene[v]
	    f.write("%s %s %s %s\n" % (network_name, u,w,v))
	for node in nodes - included_nodes:
	    f.write("%s %s\n" % (network_name, ueid_to_gene[node]))
	f.close()

	output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/breast_cancer_pruned_p80_strong.sif"
	f = open(output_file, 'w')
	included_nodes = set()
	for u,v in g_sub.edges():
	    if u == v:
		continue
	    included_nodes.add(u)
	    included_nodes.add(v)
	    w = "pp"
	    if (u,v) in weak_edges or (v,u) in weak_edges:
		w = "weak"
	    elif (u,v) in strong_edges or (v,u) in strong_edges:
		w = "strong"
	    u = ueid_to_gene[u]
	    v = ueid_to_gene[v]
	    f.write("%s %s %s\n" % (u,w,v))
	for node in nodes - included_nodes:
	    f.write("%s\n" % ueid_to_gene[node])
	f.close()

	for i, module in enumerate(module_sets):
	    output_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/M%d_.txt" % i
	    #g_sub_sub = g_sub.subgraph(module) # For checking the functions enriched in the largest connected component of the module
	    #module = gu.get_connected_components(g_sub_sub, return_as_graph_list=False)[0] 
	    f = open(output_file, 'w')
	    for node in module:
		f.write("%s\n" % ueid_to_gene[node])
	    f.close()
	    functional_enrichment.check_functional_enrichment_of_human_gene_symbols(output_file, output_file+".funcassoc")

    if False:
	for i in range(4):
	    enrichment_file = DATA_DIR + "summary_runs_on_random/breast_cancer_pruned/" + "%s/M%d_.txt.funcassoc" % ("strong/func-all", i) # strong/func-all strong weak
	    go_terms = functional_enrichment.get_functional_enrichment(enrichment_file, remove_parents=False, only_biological_processes=True)
	    print "m%d<-c(\"%s\")" % (i, "\", \"".join(go_terms))
    return


def case_study_neighborhood():
    perturbation = "pruned" #"permuted" #
    network_file = DATA_DIR + "input/biana_no_tap_relevance/edge_scores.sif"
    seeds_file = DATA_DIR + "input/biana_no_tap_relevance/omim_alzheimer/seeds.txt" 
    neighborhood_network_file = DATA_DIR + "compare/biana_no_tap_relevance-omim_alzheimer-nn/seed_neighborhood.sif"
    #get_neighbors_of_nodes_in_network(network_file, seeds_file, neighborhood_network_file)

    user_entity_id_mapping_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_relevance/node_mapping.tsv.genesymbol.single"
    #seed_genes_file = DATA_DIR + "compare/biana_no_tap_relevance-omim_alzheimer-nd-top5/seeds.txt"
    all_genes_file = DATA_DIR + "compare/biana_no_tap_relevance-omim_alzheimer-nd-top5/all.txt"
    seed_genes_file = DATA_DIR + "omim/2009_Aug_27/alzheimer.txt" 
    aging_file = DATA_DIR + "uwaging/mutex_uwaging_genage_netage.txt"
    ad_file = DATA_DIR + "alzheimer_gold/gene_list.txt"

    f = open(DATA_DIR + "compare/biana_no_tap-omim_alzheimer-nn/%s/results.dat" % perturbation, 'w')
    f.write("\tpicked_good\ttotal\tgood\tpicked\n")

    comparison = get_corresponding_genes_of_ues(neighborhood_network_file, user_entity_id_mapping_file)
    seeds = set([gene.strip() for gene in open(seed_genes_file)])
    #all_genes2 = set([gene.strip() for gene in open(all_genes_file)])
    all_genes = get_corresponding_genes_of_ues(network_file, user_entity_id_mapping_file)
    ad_genes = set([gene.strip() for gene in open(ad_file)]) & all_genes
    aging_genes = set([gene.strip() for gene in open(aging_file)]) & all_genes
    print len(seeds&aging_genes), len(all_genes&seeds), len(aging_genes)
    #print len(seeds), len(all_genes), len(comparison)
    from scipy.stats import hypergeom
    picked = comparison-seeds
    total = all_genes-seeds
    print "comparison:", len(picked)
    print "all:", len(total)

    good = ad_genes-seeds
    matched_ad = (comparison&ad_genes)-seeds
    picked_good = matched_ad
    print "ad:", len(good)
    print "matched_ad:", len(matched_ad)
    print "p_value:", sum(hypergeom.pmf(range(len(picked_good),len(picked)+1), len(total), len(good), len(picked)))
    print len(picked_good), len(total), len(good), len(picked)
    f.write("%s\t%d\t%d\t%d\t%d\n" % ("p0", len(picked_good), len(total), len(good), len(picked)))

    # Numbers with random AD-related set
    f2 = open(DATA_DIR + "compare/biana_no_tap-omim_alzheimer-nn/%s/results_random.dat" % perturbation, 'w')
    f2.write("\tpercentage\trun\tpicked_good\ttotal\tgood\tpicked\n")
    from random import shuffle 
    genes = list(total) # all_genes-seeds
    random_genes = []
    for i in xrange(1, 101):
	shuffle(genes)
	ad_genes_random = set(genes[:len(good)])
	random_genes.append(ad_genes_random)
	good = ad_genes_random
	picked_good = comparison&ad_genes_random
	f2.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % ("p0-%d" % i, 0, i, len(picked_good), len(total), len(good), len(picked)))

    good = aging_genes-seeds
    matched_aging = (comparison&aging_genes)-seeds
    picked_good = matched_aging
    print "aging:", len(good)
    print "matched_aging:", len(matched_aging)
    print "p_value:", sum(hypergeom.pmf(range(len(picked_good),len(picked)+1), len(total), len(good), len(picked)))

    print "ad-aging:", len((ad_genes&aging_genes)-seeds)
    print "matched_ad-matched_aging:", len(matched_ad&matched_aging)
    
    # Check effect of permuting/pruning
    for p in xrange(10, 90, 10):
	comparisons = []
	n_picked_good, n_total, n_good, n_picked = [0.0] * 4
	#os.mkdir(DATA_DIR + "compare/biana_no_tap-omim_alzheimer-nn/%d" % p)
	for i in xrange(1,101):
	    if perturbation == "pruned":
		network_file = DATA_DIR + "human_interactome_biana/%s/omim_alzheimer/%d/sampled_graph.sif.%d" % (perturbation, p, i)
	    else:
		network_file = DATA_DIR + "human_interactome_biana/%s/%d/sampled_graph.sif.%d" % (perturbation, p, i)
	    neighborhood_network_file = DATA_DIR + "compare/biana_no_tap-omim_alzheimer-nn/%s/%d/seed_neighborhood.sif.%d" % (perturbation, p, i)
	    #get_neighbors_of_nodes_in_network(network_file, seeds_file, neighborhood_network_file)
	    comparison = get_corresponding_genes_of_ues(neighborhood_network_file, user_entity_id_mapping_file)
	    comparisons.append(comparison)
	    #all_genes = get_corresponding_genes_of_ues(network_file, user_entity_id_mapping_file) # all networks contain all genes (even though they are unconnected)
	    ad_genes = set([gene.strip() for gene in open(ad_file)]) & all_genes
	    picked = comparison-seeds
	    total = all_genes-seeds
	    good = ad_genes-seeds
	    matched_ad = (comparison&ad_genes)-seeds
	    picked_good = matched_ad
	    #print len(picked_good), len(total), len(good), len(picked)
	    n_picked_good += len(picked_good)
	    n_total += len(total)
	    n_good += len(good)
	    n_picked += len(picked)
	n_picked_good, n_total, n_good, n_picked = map(lambda x: int(round(x/100)), [n_picked_good, n_total, n_good, n_picked])
	print n_picked_good, n_total, n_good, n_picked
	print "p_value:", sum(hypergeom.pmf(range(n_picked_good,n_picked+1), n_total, n_good, n_picked))
	f.write("%s\t%d\t%d\t%d\t%d\n" % ("p%d"%p, n_picked_good, n_total, n_good, n_picked))

	# Numbers with random AD-related set
	genes = list(total) # all_genes-seeds
	for j in xrange(1, 101):
	    n_picked_good, n_total, n_good, n_picked = [0.0] * 4
	    for i in xrange(1,101):
		if perturbation == "pruned":
		    network_file = DATA_DIR + "human_interactome_biana/%s/omim_alzheimer/%d/sampled_graph.sif.%d" % (perturbation, p, i)
		else:
		    network_file = DATA_DIR + "human_interactome_biana/%s/%d/sampled_graph.sif.%d" % (perturbation, p, i)
		neighborhood_network_file = DATA_DIR + "compare/biana_no_tap-omim_alzheimer-nn/%s/%d/seed_neighborhood.sif.%d" % (perturbation, p, i)
		#comparison = get_corresponding_genes_of_ues(neighborhood_network_file, user_entity_id_mapping_file)
		comparison = comparisons[i-1]
		ad_genes_random = random_genes[j-1]
		picked_good = comparison&ad_genes_random
		good = ad_genes_random
		#all_genes = get_corresponding_genes_of_ues(network_file, user_entity_id_mapping_file)
		picked = comparison-seeds
		total = all_genes-seeds
		n_picked_good += len(picked_good)
		n_total += len(total)
		n_good += len(good)
		n_picked += len(picked)
	    n_picked_good, n_total, n_good, n_picked = map(lambda x: int(round(x/100)), [n_picked_good, n_total, n_good, n_picked])
	    f2.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % ("p%d-%d" % (p, j), p, j, n_picked_good, n_total, n_good, n_picked))

    f.close()
    f2.close()
    return

def get_go_function_counts():

    base_dir = "/home/emre/arastirma/netzcore/data/module/biana_no_tap-omim/"

    from toolbox import OboParser, functional_enrichment
    enrichment_file = base_dir + "enrichment.txt" # enrichment for all top 5% genes not only the ones in modules
    name_to_go_terms = functional_enrichment.get_functional_enrichment(enrichment_file, remove_parents=False, only_biological_processes=True)
    #for name, go_terms in name_to_go_terms.iteritems():
    #	print name, go_terms

    g=OboParser.getOboGraph("/home/emre/arastirma/celldiff/data/GO/gene_ontology.1_2.obo")

    i, mode, disease, name, seed_terms = None, None, None, None, None
    go_terms = None

    disease_to_go_terms = {} # seed_terms, module_terms, go_terms

    f = open(base_dir + "modules_ns.txt")
    #f_out = open(base_dir + "module_summary_ns-no_parent.dat", 'w')
    #f_out.write("ppi phenotype scoring n_seed_go n_seed_go_in_modules n_go_in_modules n_go\n")
    for line in f:
	if line.strip().endswith("Seed GO terms"):
	    if go_terms is not None: 
		#if name in name_to_go_terms: # Some diseases are excluded in the functional enrichment analysis
		#    all_go_terms = name_to_go_terms[name]
		#    count = str(len(all_go_terms))
		#    disease_to_go_terms["omim_"+disease.replace(" ", "_")] = (seed_terms, go_terms, all_go_terms)
		#else:
		#    count = "NA"
		#go_terms = functional_enrichment.remove_parent_terms(go_terms, g)
		#print disease, go_terms
		print i, i_seed, len(go_terms), len(go_terms & seed_terms), len(seed_terms)
		#f_out.write("biana_no_tap_no_reliability omim_%s ns %d %d %d %s\n" % ("_".join(disease.split()), len(seed_terms), len(go_terms & seed_terms), len(go_terms), count)) #i_seed, i))
		disease_to_go_terms[disease] = (seed_terms, go_terms)
	    mode = 1
	    disease = line.split(" - ")[1]
	    #name = "#biana_no_tap_no_reliability omim_%s ns" % disease.replace(" ", "_")
	    print disease
	    #f_out.write("biana_no_tap_no_reliability omim_%s no" % "_".join(disease.split()))
	    i = 0
	    seed_terms = set()
	elif line.startswith("bPPI") and line.strip().endswith("NetScore"):
	    #print i
	    #f_out.write(" %d NA NA NA\n" % i)
	    mode = 0
	    i = 0
	    i_seed = 0
	    #print seed_terms
	    #seed_terms = functional_enrichment.remove_parent_terms(seed_terms, g) 
	    #print disease, seed_terms
	    go_terms = set()

	#if line["#"]:
	#	continue
	words = line.split("\t")
	try:
	    float(words[0])
	except:
	    #print words[0]
	    continue
	go_term = words[5]
	if mode == 1:
	    if g.node[go_term]['t'] == "biological_process": # and 'a' in g.node[go_term]: #is bp and slim # ("molecular_function", "biological_process"):
		seed_terms.add(go_term)
	else:
	    if g.node[go_term]['t'] == "biological_process":
		go_terms.add(go_term)
	    if go_term in seed_terms:
		i_seed += 1

	if g.node[go_term]['t'] == "biological_process": # and 'a' in g.node[go_term]: #in ("molecular_function", "biological_process"):
	    i += 1

    #go_terms = functional_enrichment.remove_parent_terms(go_terms, g)
    #print disease, go_terms
    print i, i_seed, len(go_terms), len(go_terms & seed_terms), len(seed_terms)
    #if name in name_to_go_terms: # Some diseases are excluded in the functional enrichment analysis
    #	all_go_terms = name_to_go_terms[name]
    #	count = str(len(all_go_terms))
    #	disease_to_go_terms["omim_"+disease.replace(" ", "_")] = (seed_terms, go_terms, all_go_terms)
    #else:
    #	count = "NA"
    #f_out.write("biana_no_tap_no_reliability omim_%s ns %d %d %d %s\n" % ("_".join(disease.split()), len(seed_terms), len(go_terms & seed_terms), len(go_terms), count)) #i_seed, i))
    #f_out.close()
    disease_to_go_terms[disease] = (seed_terms, go_terms)
    f.close()

    import networkx
    f_out = open(base_dir + "module_summary_ns-seed_go.dat", 'w')
    f_out.write("phenotype go level\n")
    f_out2 = open(base_dir + "module_summary_ns-all_go.dat", 'w')
    f_out2.write("phenotype go level\n")
    root_id = "GO:0008150" # BP
    for disease, values in disease_to_go_terms.iteritems():
	name = "#biana_no_tap_no_reliability omim_%s ns" % disease.replace(" ", "_")
	if name not in name_to_go_terms: # Some diseases are excluded in the functional enrichment analysis
	    continue
	#print disease
	all_terms = name_to_go_terms[name]
	seed_terms, go_terms = values
	for go_id in seed_terms:
	    level = len(networkx.bidirectional_shortest_path(g, go_id, root_id))
	    f_out.write("omim_%s %s %d\n" % ("_".join(disease.split()), go_id, level)) 
	for go_id in all_terms:
	    level = len(networkx.bidirectional_shortest_path(g, go_id, root_id))
	    f_out2.write("omim_%s %s %d\n" % ("_".join(disease.split()), go_id, level)) 
    f_out.close()
    f_out2.close()

    disease_to_go_terms2 = {}
    f_out = open(base_dir + "module_summary_ns-only_bp.dat", 'w')
    f_out.write("ppi phenotype scoring n_seed_go n_seed_go_in_modules n_go_in_modules n_go\n")
    for disease, values in disease_to_go_terms.iteritems():
	seed_terms, go_terms = values
	count = "NA"
	name = "#biana_no_tap_no_reliability omim_%s ns" % disease.replace(" ", "_")
	if name in name_to_go_terms: # Some diseases are excluded in the functional enrichment analysis
	    all_go_terms = name_to_go_terms[name]
	else:
	    all_go_terms = set()
	#union_terms = seed_terms | go_terms | all_go_terms
	#union_terms = functional_enrichment.remove_parent_terms(union_terms, g)
	#seed_terms &= union_terms
	#go_terms &= union_terms
	#all_go_terms &= union_terms
	#if disease == "alzheimer":
	#    print len(seed_terms), seed_terms
	#    print len(go_terms), go_terms
	#    print len(all_go_terms), all_go_terms
	if name in name_to_go_terms: 
	    count = str(len(all_go_terms))
	    disease_to_go_terms2["omim_"+disease.replace(" ", "_")] = (seed_terms, go_terms, all_go_terms)
	f_out.write("biana_no_tap_no_reliability omim_%s ns %d %d %d %s\n" % ("_".join(disease.split()), len(seed_terms), len(go_terms & seed_terms), len(go_terms), count)) 
    f_out.close()

    return disease_to_go_terms2


def output_ns_related_modules():
    f = open("modules.txt")
    f_out = open("modules_ns.txt", 'w')

    mode = 0
    for line in f:
	if line.strip().endswith("Seed GO terms"):
	    mode = 1
	elif line.startswith("biana_no_tap_no_reliability") and line.strip().endswith("ns"):
	    mode = 1
	elif line.startswith("biana_no_tap_no_reliability"):
	    mode = 0

	if mode == 1:
	    f_out.write(line)

    f.close()
    f_out.close()
    return


def get_ues_gene_mapping(mapping_file):
    ueid_to_gene = {}
    f = open(mapping_file)
    f.readline()
    for line in f:
	ueid, gene = line.strip().split("\t")
	ueid_to_gene[ueid] = gene
    f.close()
    return ueid_to_gene 


def get_neighbors_of_nodes_in_network(network_file, node_file, output_file):
    from toolbox import network_utilities as gu
    g = gu.create_network_from_sif_file(network_file, use_edge_data=True)
    nodes = [ line.strip() for line in open(node_file) ]
    neighbors = []
    for node in nodes:
	neighbors.extend(g.neighbors(node))
    neighbors.extend(nodes)
    g_sub = g.subgraph(neighbors)
    f = open(output_file, 'w')
    for u,v,w in g_sub.edges(data=True):
	f.write("%s %s %s\n" % (u,w,v))
    f.close()
    return


if __name__ == "__main__":
    main()

