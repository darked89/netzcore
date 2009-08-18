#########################################################################
# Methods for analyzing output files from scoring methods
#
# eg 13/08/2009
#########################################################################

from biana.utilities import graph_utilities as network_utilities


def calculate_seed_coverage_at_given_percentage(file_seed_scores, file_result, percentage):
    """
	Calculates number of seed nodes included in the given percentage of high ranking nodes in result file 
	Returns a tuple containing number of seed nodes and all nodes at that percentage
    """
    setDummy, setDummy, dictNode, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_seed_scores, store_edge_type = False)
    setDummy, setDummy, dictNodeResult, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_result, store_edge_type = False)
    
    result_scores = dictNodeResult.items()
    result_scores.sort(lambda x,y: cmp(x[1], y[1]))

    i = 0
    n = len(result_scores)*percentage/100
    n_seed = 0

    for id, score in result_scores:
	i+=1
	if i>n:
	    break
	# checking whether node was a seed (using initial non-seed score assumption)
	if dictNode[id] > 0.0001:
	    n_seed += 1
    
    return n_seed, n


