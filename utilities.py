import os

#DATA_DIR = "/home/emre/arastirma/netzcore/data/"
DATA_DIR = "../data/"

def main():
    #convert_ids_for_netscore_alzheimer_analysis()
    get_seed_gene_counts(DATA_DIR + "input/biana_no_tap_no_reliability/")#"input_runs_for_draft/biana_no_tap_no_reliability/")

    #file_name = DATA_DIR + "omim/alzheimer.txt"
    #out_file_name = "test.txt"
    #check_functional_enrichment(file_name, out_file_name)

    #get_common_non_empty_navlakha_files()

    #get_unique_ids_for_biana_ids("biana_ids.txt", "UniprotAccession", "biana_id_to_uniprot_ids.txt") # change type to unique
    #get_unique_ids_for_biana_ids("biana_ids.txt", "GeneID", "biana_id_to_gene_ids.txt")

    #check_pearson_correlation_between_seed_connectivity_and_performance()

    #convert_sif_file_to_R_matrix(DATA_DIR + "input_runs_for_draft/entrez/edge_scores.sif", "test_ppi.dat")

    #case_study_neighborhood()

    #get_go_terms()

    return

def get_go_terms():
    base_dir = "/home/emre/arastirma/netzcore/data/module/biana_no_tap-omim/"
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

def case_study_neighborhood():
    network_file = DATA_DIR + "input/biana_no_tap_relevance/edge_scores.sif"
    seeds_file = DATA_DIR + "input/biana_no_tap_relevance/omim_alzheimer/seeds.txt" 
    neighborhood_network_file = DATA_DIR + "compare/biana_no_tap_relevance-omim_alzheimer-nn/seed_neighborhood.sif"
    #get_neighbors_of_nodes_in_network(network_file, seeds_file, neighborhood_network_file)

    user_entity_id_mapping_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_relevance/node_mapping.tsv.genesymbol.single"
    #seed_genes_file = DATA_DIR + "compare/biana_no_tap_relevance-omim_alzheimer-nd-top5/seeds.txt"
    all_genes_file = DATA_DIR + "compare/biana_no_tap_relevance-omim_alzheimer-nd-top5/all.txt"
    seed_genes_file = DATA_DIR + "omim/alzheimer.txt" 
    aging_file = DATA_DIR + "uwaging/mutex_uwaging_genage_netage.txt"
    ad_file = DATA_DIR + "alzheimer_gold/gene_list.txt"

    f = open(DATA_DIR + "compare/biana_no_tap_relevance-omim_alzheimer-nn/results.dat", 'w')
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

    good = aging_genes-seeds
    matched_aging = (comparison&aging_genes)-seeds
    picked_good = matched_aging
    print "aging:", len(good)
    print "matched_aging:", len(matched_aging)
    print "p_value:", sum(hypergeom.pmf(range(len(picked_good),len(picked)+1), len(total), len(good), len(picked)))

    print "ad-aging:", len((ad_genes&aging_genes)-seeds)
    print "matched_ad-matched_aging:", len(matched_ad&matched_aging)

    # Check effect of pruning
    seeds_file = DATA_DIR + "input/biana_no_tap_relevance/omim_alzheimer/seeds.txt" 
    for p in xrange(10, 90, 10):
	n_picked_good, n_total, n_good, n_picked = [0.0] * 4
	#os.mkdir(DATA_DIR + "compare/biana_no_tap_relevance-omim_alzheimer-nn/%d" % p)
	for i in xrange(1,101):
	    #network_file = DATA_DIR + "human_interactome_biana/pruned/omim_alzheimer/%d/sampled_graph.sif.%d" % (p, i)
	    network_file = DATA_DIR + "human_interactome_biana/permuted/%d/sampled_graph.sif.%d" % (p, i)
	    neighborhood_network_file = DATA_DIR + "compare/biana_no_tap_relevance-omim_alzheimer-nn/%d/seed_neighborhood.sif.%d" % (p, i)
	    #get_neighbors_of_nodes_in_network(network_file, seeds_file, neighborhood_network_file)
	    comparison = get_corresponding_genes_of_ues(neighborhood_network_file, user_entity_id_mapping_file)
	    all_genes = get_corresponding_genes_of_ues(network_file, user_entity_id_mapping_file)
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
    f.close()
    return

def get_corresponding_genes_of_ues(network_file, mapping_file, network_file_delim=" "):
    ueid_to_gene = {}
    f = open(mapping_file)
    f.readline()
    for line in f:
	ueid, gene = line.strip().split("\t")
	ueid_to_gene[ueid] = gene
    f.close()
    comparison = set()
    for line in open(network_file):
	words = line.strip().split(network_file_delim)
	if len(words) == 3:
	    id1, score, id2 = words
	else:
	    id1 = words[0]
	    id2 = id1
	for ueid in [id1, id2]:
	    gene = ueid_to_gene[ueid]
	    if gene != "-":
		comparison.add(gene)
    return comparison

def get_neighbors_of_nodes_in_network(network_file, node_file, output_file):
    import biana.utilities.graph_utilities as gu
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

def convert_sif_file_to_R_matrix(file_name, out_file_name):
    id_pair_to_score = {}
    ids = set()
    for line in open(file_name):
	id1, score, id2 = line.strip().split()
	ids.add(id1)
	ids.add(id2)
	id_pair_to_score[(id1, id2)] = score
	id_pair_to_score[(id2, id1)] = score
    ids = sorted(list(ids))
    f = open(out_file_name, 'w')
    f.write(" %s\n" % " ".join(ids))
    for id1 in ids:
	f.write("%s" % id1)
	for id2 in ids:
	    if (id1, id2) in id_pair_to_score:
		f.write(" %s" % id_pair_to_score[(id1, id2)])
	    else:
		f.write(" 0")
	f.write("\n")
    f.close()
    return

def check_pearson_correlation_between_seed_connectivity_and_performance():
    from scipy import stats
    method_to_values = read_R_data_file(DATA_DIR + "summary/biana_no_tap7_vs_all5-top5/auc_ppis.dat")
    metric_to_values = read_R_data_file(DATA_DIR + "summary/biana_no_tap7_vs_all5-top5/seeds.dat")
    for method, values1 in method_to_values.iteritems():
	values1.sort()
	for metric, values2 in metric_to_values.iteritems():
	    values2.sort()
	    assert map(lambda x: x[0], values1) == map(lambda x: x[0], values2)
	    vals1 = map(lambda x: x[1], values1)
	    vals2 = map(lambda x: x[1], values2)
	    print method, metric, "%.2lf %e" % stats.pearsonr(vals1, vals2)
    return

def read_R_data_file(file_name, seperator="\t"):
    f = open(file_name)
    headers = f.readline().rstrip().split(seperator)
    header_to_values = {}
    for line in f:
	words = line.rstrip().split(seperator)
	rowname=None
	for header, value in zip(headers, words):
	    if header == "":
		rowname = value
		continue
	    header_to_values.setdefault(header, []).append((rowname, float(value)))
    return header_to_values

def get_unique_ids_for_biana_ids(file_name, attribute, out_file):
        biana_ids = [ line.rstrip() for line in open(file_name) ]
        fetch_unique_ids(biana_ids, attribute, out_file)

def fetch_unique_ids(biana_ids, attribute, out_file):
        import MySQLdb
        db = MySQLdb.connect("localhost", "", "", "test_biana")
        c=db.cursor()
        #c.execute("SELECT U.userEntityID, E.value FROM userEntityUnification_protocol_12 U, externalEntity%s E where U.externalEntityID=E.externalEntityID AND E.type='unique' AND U.userEntityID IN (%s)" % (attribute, ",".join(biana_ids)))
        c.execute("SELECT U.userEntityID, E.value FROM userEntityUnification_protocol_12 U, externalEntity%s E where U.externalEntityID=E.externalEntityID AND E.type='cross-reference' AND U.userEntityID IN (%s)" % (attribute, ",".join(biana_ids)))
        f = open(out_file, 'w')
        for row in c.fetchall():
                f.write("%s\n" % "\t".join(map(str, row)))
        f.close()
        return


def get_common_non_empty_navlakha_files():
    associations = get_non_empty_files_in_folder(DATA_DIR + "navlakha/associations/")
    candidates = get_non_empty_files_in_folder(DATA_DIR + "navlakha/candidates/") # Handle empty candidate files
    print len(associations), len(candidates), len(set(associations) & set(candidates))
    print '"%s"' % '", "'.join(sorted(set(associations) & set(candidates)))
    return

def get_non_empty_files_in_folder(base_dir):	
    import os, os.path
    files = []
    for name in os.listdir(base_dir):
	f = open(base_dir + name)
	if f.readline().strip() == "":
	    continue
	else:
	    files.append(name[:-4])
	    f.close()
    return files


def get_seed_gene_counts(input_dir):
    import os, os.path
    for name in os.listdir(input_dir):
	if not os.path.isdir(input_dir + name):
	    continue
	filename = input_dir + name + os.sep + "README"
	if not os.path.exists(filename):
	    continue
	print name, "\t",
	f = open(filename)
	n_seed = None
	for line in f.readlines():
	    if line.startswith("Covered genes (seed genes):"):
		n_seed = int(line.split(":")[1].split()[2])
		print n_seed
		break
	f.close()
    return

def convert_ids_for_netscore_alzheimer_analysis():
    node_mapping_file = DATA_DIR + "input_runs_for_draft/biana_no_tap_relevance/node_mapping.tsv.genesymbol"
    convert_multiple_entry_to_single(node_mapping_file)
    #node_association_file = DATA_DIR + "uwaging/aging.txt"
    #node_association_file = DATA_DIR + "genage/aging_candidates.txt"
    #node_association_file = DATA_DIR + "netage/longetivity.txt"
    #node_association_file = DATA_DIR + "netage/AD_genes.txt"
    #node_association_file = DATA_DIR + "alzheimer_gold/gene_list.txt"
    #node_association_file = DATA_DIR + "uwaging/uwaging_mutex_genage_netage.txt"
    #node_association_file = DATA_DIR + "uwaging/mutex_uwaging_genage_netage.txt"
    node_association_file1 = DATA_DIR + "uwaging/uwaging_genage_netage.txt"
    for node_association_file in [ node_association_file1 ]: #, node_association_file2, node_association_file3 ]:
    	convert_genesymbol_to_user_entity_id(node_mapping_file + ".single.reversed", node_association_file)


def convert_genesymbol_to_user_entity_id(node_mapping_file, node_association_file):
    from analyze_results import get_id_to_mapped_id_mapping
    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)
    f = open(node_association_file + ".user_entity_id", 'w')
    for line in open(node_association_file):
	gene = line.strip()
	if gene in id_to_mapped_ids:
	    f.write("%s\n" % id_to_mapped_ids[gene][0])
    f.close()
    return

def convert_multiple_entry_to_single(node_mapping_file):
    from analyze_results import get_id_to_mapped_id_mapping
    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)
    f = open(node_mapping_file + ".single", 'w')
    f.write("user entity id\tgene symbol\n")
    f2 = open(node_mapping_file + ".single.reversed", 'w')
    f2.write("gene symbol\tuser entity id\n")
    for id, mapped_ids in id_to_mapped_ids.iteritems():
	f.write("%s\t%s\n" % (id, mapped_ids[0]))
	f2.write("%s\t%s\n" % (mapped_ids[0], id))
    f.close()
    f2.close()
    return

def rename(base_dir, prefix, refix):
    import os
    #base_dir = "summary"
    #prefix = "biana_no_tap_no_reliability_pruned_p"
    #refix = "biana_no_tap_no_reliability_pruned_non_seed_interactions_p"
    #prefix = "biana_no_tap_pruned_p"
    #refix = "biana_no_tap_pruned_non_seed_interactions_p"

    for dir in os.listdir(base_dir):
	if dir.startswith(prefix):
	    new_dir = refix + dir[len(prefix):]
	    print "renaming", dir, new_dir
	    os.rename(os.path.join(base_dir, dir), os.path.join(base_dir, new_dir))
    return


def check_functional_enrichment(file_name, out_file_name):
    from analyze_results import check_functional_enrichment
    f=open(file_name)
    lines = f.readlines()
    a=map(lambda x: x.strip(), lines)
    check_functional_enrichment(a, None, "genesymbol", open(out_file_name, 'w').write, tex_format = False) #True)
    return

def prune_network():
    import prepare_data
    g=prepare_data.network_utilities.networkx.Graph()
    g.add_edges_from([(1,2),(2,3),(2,4),(2,6),(3,4),(4,5),(5,6)])
    g2=prepare_data.network_utilities.prune_graph_at_given_percentage(g,10)
    print g2.edges()
    return

def analyze_network():
    import biana.utilities.graph_utilities as gu
    g = gu.create_network_from_sif_file("/data/emre/toy_data/test_interactions_small.sif")
    degrees = g.degree(with_labels=True)
    node_to_values = gu.get_node_degree_related_values(g, ["v2","v3"])
    for v in g.nodes():
	print v, degrees[v], node_to_values[v]
    gu.create_R_analyze_network_script(g, ["v2","v3"])

#node_file = DATA_DIR + "output/biana/aneurysm/nb/i5/node_scores.sif"
#node_file_geneid = DATA_DIR + "output/biana/aneurysm/ns/r3i3/node_scores.sif.geneid"
#node_file_geneid = DATA_DIR + "output/goh/aneurysm/nb/i5/node_scores.sif"
#node_file_genesymbol = DATA_DIR + "output/biana/aneurysm/ns/r3i3/node_scores.sif.genesymbol"
#node_file_genesymbol = DATA_DIR + "output/goh/aneurysm/nb/i5/node_scores.sif.genesymbol"

#import prepare_data
#prepare_data.convert_file_using_new_id_mapping(node_file, "/home/emre/arastirma/netzcore/data/human_interactome_biana/human_nodes.tsv", "user entity id", "geneid", node_file_geneid)
#prepare_data.convert_file_using_new_id_mapping(node_file_geneid, "/home/emre/arastirma/netzcore/data/gene_info/genes.tsv", "geneid", "genesymbol", node_file_genesymbol)


#prepare_data.convert_file_using_new_id_mapping("/home/emre/arastirma/netzcore/data/input/biana/edge_scores.sif", "/home/emre/arastirma/netzcore/data/human_interactome_biana/human_nodes.tsv", "user entity id", "geneid", "test_mapping4.txt")

#from biana.utilities import biana_output_converter 
#print biana_output_converter.get_primary_uniprot_accessions("/home/emre/arastirma/netzcore/data/human_interactome_biana/human_nodes.tsv")

def get_detection_methods():
    create_new_session(sessionID="biana_session",dbname="test_biana",dbhost="127.0.0.1",dbuser="",dbpassword="",unification_protocol="(p3)(noself)(noprevious)uniprot_seqtax_geneid_scoppdb")
    objSession = available_sessions["biana_session"]
    temp = available_sessions["biana_session"].create_new_user_entity_set( identifier_description_list = [("geneid","390")], attribute_restriction_list=[], id_type="embedded", new_user_entity_set_id="User_Entity_Set_1", negative_attribute_restriction_list=[])
    temp = available_sessions["biana_session"].output_user_entity_set_details(user_entity_set_id="User_Entity_Set_1", output_format="xml", include_command_in_rows=True, attributes = ["uniprotaccession"], only_selected=False, output_1_value_per_attribute=False, output_only_native_values=False, output_only_unique_values=True)

    y2h = set(map(int, objSession.get_ontology(ontology_name="psimiobo", root_attribute_values = [18]).get_all_linked_attributes()))
    tap = set(map(int, objSession.get_ontology(ontology_name="psimiobo", root_attribute_values = [4]).get_all_linked_attributes()))
    tap.add(59)
    tap.add(109)
    experimental = set(map(int, objSession.get_ontology(ontology_name="psimiobo", root_attribute_values = [45]).get_all_linked_attributes()))
    predicted = set(map(int, objSession.get_ontology(ontology_name="psimiobo", root_attribute_values = [63]).get_all_linked_attributes()))

    print y2h, tap, experimental, predicted

if __name__ == "__main__":
    main()

