#########################################################################
# Methods for analyzing output files from scoring methods
#
# eg 13/08/2009
#########################################################################

from biana.utilities import graph_utilities as network_utilities


def calculate_performance_metric_counts(file_result, file_seed_test_scores, file_seed_scores, score_threshold, n_random_negative_folds = 1, default_score = 0):
    """
	Calculate and return TP, FP, TN, FN (FP and TN based on random selected non-seed nodes)
	file_result: output scores file
	file_seed_test_scores: seed nodes separated for test
	file_seed_scores: initial node/seed scores
	score_threshold: threshold for considering data as P, N
	n_random_negative_folds: Number of negative data selection folds (each fold contain same number of test nodes) to be averaged for FP and TN calculation. FP and TN are not necesserily integers when n_random_negative_folds > 1
    """
    setNodeResult, setDummy, dictNodeResult, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_result, store_edge_type = False)
    setNodeTest, setDummy, dictNodeTest, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_seed_test_scores, store_edge_type = False)
    setNode, setDummy, dictNode, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_seed_scores, store_edge_type = False)

    listNonSeed = [ id for id, score in dictNode.iteritems() if score > default_score ]

    (nTP, nFP, nFN, nTN) = (0.0, 0.0, 0.0, 0.0)
    for id, score in dictNodeResult.iteritems():
        if id in setNodeTest:
            if score >= score_threshold:
                nTP += 1
            else:
                nFN += 1

    from random import shuffle
    shuffle(listNonSeed)
    n_actual_folds = 0
    for i in range(n_random_negative_folds):
	if (i+1)*len(setNodeTest) > len(listNonSeed):
	    break
	setNegative = set( listNonSeed[i*len(setNodeTest):(i+1)*len(setNodeTest)] )
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


def calculatePerformance(nTP, nFP, nFN, nTN):
    acc = (nTP + nTN) / (nTP + nFP + nTN + nFN)
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


def calculate_seed_coverage_at_given_percentage(file_result, file_seed_scores, percentage, default_score=0):
    """
	Calculates number of seed nodes included in the given percentage of high ranking nodes in result file 
	Returns a tuple containing number of seed nodes and all nodes at that percentage
    """
    setDummy, setDummy, dictNodeResult, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_result, store_edge_type = False)
    setDummy, setDummy, dictNode, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_seed_scores, store_edge_type = False)
    
    result_scores = dictNodeResult.items()
    result_scores.sort(lambda x,y: cmp(y[1], x[1]))

    i = 0
    n = len(result_scores)*percentage/100
    n_seed = 0
    
    last_score = result_scores[0][1]
    for id, score in result_scores:
	i+=1
	if i>n:
	    if last_score != score:
		#print "i,n:", i, n
		#print "id,score,last:", id, score, last_score
		break
	# checking whether node was a seed (using initial non-seed score assumption)
	if dictNode.has_key(id) and dictNode[id] > default_score:
	    n_seed += 1
	last_score = score
    
    return n_seed, n, i


