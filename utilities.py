import os

#DATA_DIR = "/home/emre/arastirma/netzcore/data/"
DATA_DIR = "../data/"

def main():
    #convert_ids_for_netscore_alzheimer_analysis()
    #get_seed_gene_counts(DATA_DIR + "input/biana_no_tap_no_reliability/")#"input_runs_for_draft/biana_no_tap_no_reliability/")

    #file_name = DATA_DIR + "omim/alzheimer.txt"
    #out_file_name = "test.txt"
    #check_functional_enrichment(file_name, out_file_name)

    #get_common_non_empty_navlakha_files()

    #get_unique_ids_for_biana_ids("biana_ids.txt", "UniprotAccession", "biana_id_to_uniprot_ids.txt") # change type to unique
    #get_unique_ids_for_biana_ids("biana_ids.txt", "GeneID", "biana_id_to_gene_ids.txt")

    #check_pearson_correlation_between_seed_connectivity_and_performance()

    #convert_sif_file_to_R_matrix(DATA_DIR + "input_runs_for_draft/entrez/edge_scores.sif", "test_ppi.dat")

    #case_study_high_scoring()

    # DATA_DIR + output/biana_no_tap_relevance/new_omim_diabetes/nc3/node_scores.sif
    #convert_ueid_scores_to_gene_scores("node_scores.sif", "output.txt")

    #get_omim_disease_similarity()
    #get_omim_disease_similarity_in_network()
    #get_omim_disease_similarity_in_network_extended()

    #get_drugs_by_targets()

    #get_average_age_of_disease_genes_in_network_extended()

    get_detection_methods()
    return

def get_average_age_of_disease_genes_in_network_extended():
    # Prerequist: run analyze for this phenotypes in score_with_original_seeds mode
    # That will convert and copy output files under the directory with proper names: extended/omim_xxx.txt
    from parse_omim import get_disease_genes
    dir_name = DATA_DIR + "omim/2009_Aug_27/extended/"
    protein_age_file = "/home/emre/arastirma/data/collaboration/macarena/edat_humanprots_ev62.sort"
    gene_to_ensembl_file = "/home/emre/arastirma/data/collaboration/macarena/mart_export.txt"
    age_category = [ "Eukarya", "Metazoans", "Vertebrates", "Mammals"] #, "Human-specific" ]
    output_file = DATA_DIR + "omim/2009_Aug_27/extended/age_category.dat"

    node_scores_file = "/home/emre/arastirma/data/collaboration/arcadi/uncommon-max/node_scores.sif.genesymbol.single"
    output_file = "/home/emre/arastirma/data/collaboration/arcadi/uncommon-max/age_category.dat"
    
    ensembl_to_category = dict([ line.strip().split("\t") for line in open(protein_age_file) ])
    gene_to_category = {}
    f = open(gene_to_ensembl_file)
    f.readline()
    for line in f:
	gene, ensembl = line.strip().split("\t")
	if ensembl not in ensembl_to_category:
	    continue
	category = ensembl_to_category[ensembl]
	if category == "Human-specific": 
	    category = "Mammals"
	if gene in gene_to_category and category != gene_to_category[gene]:
	    print "Warning: inconsistent category for", gene, category
	gene_to_category[gene] = category
    f.close()

    if node_scores_file is not None:
	f = open(output_file, 'w')
	f.write("gene score category\n")
	for line in open(node_scores_file):
	    gene, score = line.split()
	    if gene in gene_to_category:
		category = gene_to_category[gene]
	    else:
		category = "NA"
	    f.write("%s %s %s\n" % (gene, score, category))
	f.close()
	return

    f = open(output_file, 'w')
    f.write("%s\n" % " ".join(age_category))
    disease_to_genes = get_disease_genes(dir_name, top_percentage=5) 
    diseases = disease_to_genes.keys()
    diseases.sort()
    for disease in diseases:
	i = 0.0 
	category_to_count = dict([(category, 0) for category in age_category ])
	genes = disease_to_genes[disease]
	for gene in genes:
	    if gene in gene_to_category:
		i += 1
		category_to_count[gene_to_category[gene]] += 1
	f.write("omim_%s" % disease)
	category_counts = []
	for category, count in category_to_count.iteritems():
	    category_counts.append((count, category))
	category_counts.sort()
	category_counts.reverse()
	for count, category in category_counts:
	    f.write(" %f" % (count/i))
	f.write("\n")
    f.close()
    return

def get_drugs_by_targets():
    from parse_drugbank import get_drug_targets
    drugbank_file = "/home/emre/arastirma/data/disease/drugbank/drugbank.xml" 
    drug_to_targets, drug_to_description, drug_to_indication = get_drug_targets(drugbank_file)

    def output_drugs_by_targets(disease_to_genes, output_file):
	f = open(output_file, 'w')
	f.write("Drug\tDescription\tIndication\tTargets\tRepurposes\tAffected target ratio for repurposes\n")
	for drug, targets in drug_to_targets.iteritems():
	    description = drug_to_description[drug]
	    if description is not None:
		description = description.encode('ascii', 'replace').replace("\n", " ").replace("\r", "")
	    else: 
		description = ""
	    indication = drug_to_indication[drug]
	    if indication is not None:
		indication = indication.encode('ascii', 'replace').replace("\n", " ").replace("\r", "")
	    else: 
		indication = ""
	    diseases = []
	    for disease, genes in disease_to_genes.iteritems():
		ratio = len(targets & genes) / float(len(targets))
		if len(targets) > 1 and ratio > 0.65:
		    diseases.append((ratio, disease))
	    if len(diseases) == 0:
		continue
	    diseases.sort()
	    diseases.reverse()
	    #f.write("%s\t%s\t%s\t%s\t%s\n" % (drug, description, indication, " | ".join(targets), " | ".join(diseases)))
	    f.write("%s\t%s\t" % (drug, description))
	    f.write("%s\t%s\t" % (indication, " | ".join(targets)))
	    f.write("%s\t%s\n" % (" | ".join(zip(*diseases)[1]), " | ".join(map(str, zip(*diseases)[0]))))
	f.close()
	return

    from parse_omim import get_disease_genes
    dir_name = DATA_DIR + "omim/2009_Aug_27/" 
    disease_to_genes = get_disease_genes(dir_name) 
    output_file = "drug_to_repurpose.txt" 
    output_drugs_by_targets(disease_to_genes, output_file)

    dir_name = DATA_DIR + "omim/2009_Aug_27/extended/"
    disease_to_genes = get_disease_genes(dir_name, top_percentage=1)
    output_file = "drug_to_repurpose_extended.txt"
    output_drugs_by_targets(disease_to_genes, output_file)
    return

def get_omim_disease_similarity():
    from parse_omim import get_disease_similarity_matrix
    dir_name = DATA_DIR + "omim/2009_Aug_27/"
    get_disease_similarity_matrix(dir_name, dir_name + "similarity.dat")
    return

def get_omim_disease_similarity_in_network():
    from parse_omim import get_disease_similarity_matrix
    dir_name = DATA_DIR + "omim/2009_Aug_27/"
    network_file = DATA_DIR + "input/biana_no_tap_no_reliability/edge_scores.sif"
    node_mapping_file = DATA_DIR + "input/biana_no_tap_no_reliability/node_mapping.tsv.genesymbol"
    convert_multiple_entry_to_single(node_mapping_file)
    genes_in_network = get_corresponding_genes_of_ues(network_file, node_mapping_file + ".single")
    get_disease_similarity_matrix(dir_name, dir_name + "similarity_in_ppi.dat", genes_to_be_considered = genes_in_network)
    return

def get_omim_disease_similarity_in_network_extended():
    # Prerequist: run analyze for this phenotypes in score_with_original_seeds mode
    # That will convert and copy output files under the directory with proper names: extended/omim_xxx.txt
    from parse_omim import get_disease_similarity_matrix
    dir_name = DATA_DIR + "omim/2009_Aug_27/extended/"
    get_disease_similarity_matrix(dir_name, dir_name + "similarity_top5.dat", top_percentage=5)
    return

def convert_ueid_scores_to_gene_scores(node_scores_file, output_scores_file):
    node_mapping_file = DATA_DIR + "input/biana_no_tap_relevance/node_mapping.tsv.genesymbol.single"
    from analyze_results import get_id_to_mapped_id_mapping
    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)
    node_to_score = {}
    for line in open(node_scores_file):
	ueid, score = line.strip().split()
	node_to_score[ueid] = score
    gene_to_score = {}
    for node, score in node_to_score.iteritems():
	if node not in id_to_mapped_ids: # ueid to "-" mapping
	    continue
	gene = id_to_mapped_ids[node][0]
	if len(id_to_mapped_ids[node])>1:
	    print id_to_mapped_ids[node]
	if gene in gene_to_score:
	    if score > gene_to_score[gene]:
		gene_to_score[gene] = score
	else:
	    gene_to_score[gene] = score
    f = open(output_scores_file, 'w')
    for gene, score in gene_to_score.iteritems():
	f.write("%s\t%s\n" % (gene, score))
    f.close()
    values = gene_to_score.items()
    values.sort(lambda x,y: cmp(y[1],x[1]))
    f = open(output_scores_file+".ranked", 'w')
    for i, val in enumerate(values):
	f.write("%s\t%s\n" % (val[0], i+1))
    f.close()
    return

def case_study_high_scoring():
    base_dir = "/home/emre/arastirma/netzcore/doc/draft/alzheimer_case_study/"
    network_file = DATA_DIR + "input/biana_no_tap_relevance/edge_scores.sif"
    node_mapping_file = DATA_DIR + "input/biana_no_tap_relevance/node_mapping.tsv.genesymbol.single"
    seeds_file = DATA_DIR + "input/biana_no_tap_relevance/new_omim_alzheimer/seed_scores.sif" 
    candidates_file = DATA_DIR + "input/biana_no_tap_relevance/new_omim_alzheimer/candidates.txt" 
    ad_omim_file = DATA_DIR + "omim/2011_Nov_2/new_omim_alzheimer.txt" 
    aging_file = DATA_DIR + "uwaging/mutex_uwaging_genage_netage.txt"
    ad_file = DATA_DIR + "alzheimer_gold/gene_list.txt"

    f = open(base_dir + "weighted_bPPI.eda", 'w')
    f.write("weight\n")
    nodes_in_network = set()
    for line in open(network_file):
	words = line.strip().split()
	f.write("%s (%s) %s = %s\n" % (words[0], words[1], words[2], words[1]))
	nodes_in_network.add(words[0])
	nodes_in_network.add(words[2])
    f.close()

    f = open(base_dir + "seeds_in_network.txt", 'w')
    for line in open(seeds_file):
	f.write("%s\n" % line.strip().split()[0])
    f.close()

    convert_genesymbol_to_user_entity_id(node_mapping_file + ".reversed", ad_omim_file, base_dir + "seed_nodes.txt")
    convert_genesymbol_to_user_entity_id(node_mapping_file + ".reversed", aging_file, base_dir + "aging_nodes.txt")
    convert_genesymbol_to_user_entity_id(node_mapping_file + ".reversed", ad_file, base_dir + "ad_gold_nodes.txt")
    seed_nodes = set([gene.strip() for gene in open(base_dir + "seed_nodes.txt")]) & nodes_in_network
    aging_nodes = set([gene.strip() for gene in open(base_dir + "aging_nodes.txt")]) & nodes_in_network
    ad_nodes = set([gene.strip() for gene in open(base_dir + "ad_gold_nodes.txt")]) & nodes_in_network
    candidate_nodes = set([gene.strip() for gene in open(candidates_file)]) & nodes_in_network
    
    f = open(base_dir + "node_phenotypes.noa", 'w')
    f.write("phenotype\n")
    common = seed_nodes & aging_nodes & ad_nodes
    for node in common:
	f.write("%s = %s\n" % (node, "seed-ad-aging"))
    for node in (seed_nodes & ad_nodes)-common:
	f.write("%s = %s\n" % (node, "seed-ad"))
    for node in (seed_nodes & aging_nodes)-common:
	f.write("%s = %s\n" % (node, "seed-aging"))
    for node in (ad_nodes & aging_nodes)-common:
	f.write("%s = %s\n" % (node, "ad-aging"))
    for node in (seed_nodes - ad_nodes) - aging_nodes:
	f.write("%s = %s\n" % (node, "seed"))
    for node in (ad_nodes - seed_nodes) - aging_nodes:
	f.write("%s = %s\n" % (node, "ad"))
    for node in (aging_nodes - seed_nodes) - ad_nodes:
	f.write("%s = %s\n" % (node, "aging"))
    for node in candidate_nodes - (seed_nodes | ad_nodes | aging_nodes):
	f.write("%s = %s\n" % (node, "candidate"))
    f.close()

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


def get_corresponding_genes_of_ues(network_file, mapping_file, network_file_delim=" "):
    ueid_to_gene = get_ues_gene_mapping(mapping_file)
    comparison = set()
    for line in open(network_file):
	words = line.strip().split(network_file_delim)
	if len(words) == 3:
	    id1, score, id2 = words
	else:
	    id1 = words[0]
	    id2 = id1
	for ueid in [id1, id2]:
	    if ueid not in ueid_to_gene:
		continue 
	    gene = ueid_to_gene[ueid]
	    if gene != "-":
		comparison.add(gene)
    return comparison


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
    method_to_values = read_R_data_file(DATA_DIR + "summary/biana_no_tap_vs_all-wo_LI/auc_ppis.dat")
    metric_to_values = read_R_data_file(DATA_DIR + "summary/biana_no_tap_vs_all-wo_LI/seeds.dat")
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
        fetch_attributes_for_biana_ids(biana_ids, attribute, out_file, unique=True)

def fetch_attributes_for_biana_ids(biana_ids, attribute, out_file, unique=False):
        import MySQLdb
        #db = MySQLdb.connect("localhost", "jgarcia", "", "biana_server_release4")
        db = MySQLdb.connect("localhost", "", "", "test_biana")
        c=db.cursor()
        if unique:
                query = "SELECT U.userEntityID, E.value FROM userEntityUnification_protocol_12 U, externalEntity%s E where U.externalEntityID=E.externalEntityID AND E.type='unique' AND U.userEntityID IN (%s)" % (attribute, ",".join(biana_ids))
        else:
                query = "SELECT U.userEntityID, E.value FROM userEntityUnification_protocol_12 U, externalEntity%s E where U.externalEntityID=E.externalEntityID AND U.userEntityID IN (%s)" % (attribute, ",".join(biana_ids))
        c.execute(query)
        biana_id_to_values = {}
        for row in c.fetchall():
                biana_id_to_values.setdefault(row[0], set()).add(row[1])
        f = open(out_file, 'w')
        for biana_id, values in biana_id_to_values.iteritems():
                f.write("%s\t%s\n" % (biana_id, "|".join(map(str, values))))
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
    	convert_genesymbol_to_user_entity_id(node_mapping_file + ".single.reversed", node_association_file, node_association_file + ".user_entity_id")


def convert_genesymbol_to_user_entity_id(node_mapping_file, node_association_file, out_file):
    from analyze_results import get_id_to_mapped_id_mapping
    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)
    f = open(out_file, 'w')
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
    from toolbox import network_utilities as gu
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
    from biana import create_new_session, available_sessions
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

