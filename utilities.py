
DATA_DIR = "/home/emre/arastirma/netzcore/data/"

def main():
    #convert_ids_for_netscore_alzheimer_analysis()
    #get_seed_gene_counts(DATA_DIR + "input_runs_for_draft/biana_no_tap_no_reliability/")

    #file_name = DATA_DIR + "omim/alzheimer.txt"
    #out_file_name = "test.tex"
    #check_functional_enrichment(file_name, out_file_name)

    #get_common_non_empty_navlakha_files()

    #get_unique_ids_for_biana_ids("biana_ids.txt", "UniprotAccession", "biana_id_to_uniprot_ids.txt") # change type to unique
    #get_unique_ids_for_biana_ids("biana_ids.txt", "GeneID", "biana_id_to_gene_ids.txt")

    #check_pearson_correlation_between_seed_connectivity_and_performance()

    convert_sif_file_to_R_matrix(DATA_DIR + "input_runs_for_draft/entrez/edge_scores.sif", "test_ppi.dat")
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
    method_to_values = read_R_data_file(DATA_DIR + "summary/biana_no_tap-all/auc_ppis.dat")
    metric_to_values = read_R_data_file(DATA_DIR + "summary/biana_no_tap-all/seeds.dat")
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
    import sys

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

