import prepare_data

#node_file = "/home/emre/arastirma/netzcore/data/output/biana/aneurysm/nb/i5/node_scores.sif"
#node_file_geneid = "/home/emre/arastirma/netzcore/data/output/biana/aneurysm/ns/r3i3/node_scores.sif.geneid"
#node_file_geneid = "/home/emre/arastirma/netzcore/data/output/goh/aneurysm/nb/i5/node_scores.sif"
#node_file_genesymbol = "/home/emre/arastirma/netzcore/data/output/biana/aneurysm/ns/r3i3/node_scores.sif.genesymbol"
#node_file_genesymbol = "/home/emre/arastirma/netzcore/data/output/goh/aneurysm/nb/i5/node_scores.sif.genesymbol"

#prepare_data.convert_file_using_new_id_mapping(node_file, "/home/emre/arastirma/netzcore/data/human_interactome_biana/human_nodes.tsv", "user entity id", "geneid", node_file_geneid)
#prepare_data.convert_file_using_new_id_mapping(node_file_geneid, "/home/emre/arastirma/netzcore/data/gene_info/genes.tsv", "geneid", "genesymbol", node_file_genesymbol)


#prepare_data.convert_file_using_new_id_mapping("/home/emre/arastirma/netzcore/data/input/biana/edge_scores.sif", "/home/emre/arastirma/netzcore/data/human_interactome_biana/human_nodes.tsv", "user entity id", "geneid", "test_mapping4.txt")

#from biana.utilities import biana_output_converter 
#print biana_output_converter.get_primary_uniprot_accessions("/home/emre/arastirma/netzcore/data/human_interactome_biana/human_nodes.tsv")

import sys
try: from biana import *
except: sys.exit(10)

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

