#########################################################################
# Scoring methods
#
# eg 25/06/2009
#########################################################################

import networkx, cPickle
from biana.utilities import graph_utilities as network_utilities


def score_by_random_model(file_node_scores, file_edge_scores, file_output_scores, default_score = 0, max_score = 1):
    from random import uniform
    f = open(file_output_scores, "w")
    setNode, setDummy, dictDummy, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_node_scores, store_edge_type = False)
    for id in setNode:
	f.write("%s\t%f\n" % (id, uniform(default_score, max_score)))
    f.close()
    return


def create_network_from_edge_file(edge_file_weights):
    g = network_utilities.create_network_from_sif_file(network_file = edge_file_weights[:-3]+"sif", weighted=True)
    return g

def create_network_from_weight_and_score_files(edge_file_weights, edge_file_scores):
    g = network_utilities.create_network_from_sif_file(network_file = edge_file_weights[:-3]+"sif", weighted=True)
    setNode, setEdge, dictNode, dictEdge = network_utilities.get_nodes_and_edges_from_sif_file(file_name = edge_file_scores[:-3]+"sif", store_edge_type = True)
    for e, s in dictEdge.iteritems():
	u,v=e
	s = float(s)
	w = g.get_edge(u,v)
	s = s*100 + 1
	#if s == 0:
	#    s = 0.1
	w /= s
	#print w
	g.add_edge(u,v,w)
    return g

#def run_and_save_results(edge_file_weights, edge_file_scores, dump_file):
def run_and_save_results(edge_file_weights, dump_file):
    # Edge weights should be handled in data preperation phase
    #g = create_network_from_weight_and_score_files(edge_file_weights = edge_file_weights, edge_file_scores = edge_file_scores)
    g = create_network_from_weight_and_score_files(edge_file_weights = edge_file_weights)
    node_to_score = score_by_shortest_paths(g) 
    f = open(dump_file, "w")
    cPickle.dump(node_to_score, f)
    f.close()
    return node_to_score

def test_run(edge_file_weights):
    from time import clock
    t1 = clock()
    g = create_network_from_edge_file(edge_file_weights = edge_file_weights)
    print map(len, networkx.connected_components(g))
    t2 = clock()
    print "Time to load graph: ", t2-t1
    node_to_score = {}
    v = "25604"
    node_to_distance, node_to_path = networkx.single_source_dijkstra (g, v)
    node_to_score[v] = 10000/reduce(lambda x,y: x+y, node_to_distance.values())
    return node_to_score[v]

def run_and_assess_performance_of_folds(k, edge_file_weights, edge_file_scores_prefix, node_file_scores, node_file_scores_prefix, result_file_prefix):
    setNode, setDummy, dictNode, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_file_scores[:-3]+"sif")
    dictNodeToListScore ={}
    for i in xrange(k):
	setNodeTrain, setDummy, dictNodeTrain, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_file_scores_prefix+"_%i.sif" % i)
	seeds_test = set()
	seeds = set()
	for u, s in dictNode.iteritems():
	    if dictNodeTrain[u] != s:
		seeds_test.add(u)
	    if float(s) > 0:
		seeds.add(u)
	node_to_scores = run_and_save_results(edge_file_weights = edge_file_weights, edge_file_scores = edge_file_scores_prefix+"_%i.sif" % i, dump_file = result_file_prefix+"_%i.dump" % i)
	for u, s in node_to_scores.iteritems():
	    if u in seeds: label = 1
	    else: label = 0
	    if u not in seeds_test: label = 0
	    dictNodeToListScore.setdefault(u, []).append((s, label))
    #print dictNodeToListScore
    createROCRPredictionsData(dictNodeToListScore, result_file_prefix+"_predictions.txt", result_file_prefix+"_labels.txt")
    return

def score_by_shortest_paths(g):
    node_to_score = {}
    #i = 0
    for v in g.nodes_iter():
	node_to_distance, node_to_path = networkx.single_source_dijkstra (g, v)
	node_to_score[v] = 10000/reduce(lambda x,y: x+y, node_to_distance.values())
	#i +=1
	#if i > 5:
	#    break
    return node_to_score

def createROCRPredictionsData(dictNodeToListScore, fileNameOut_predictions, fileNameOut_labels):
    fileOut = open(fileNameOut_predictions, 'w')
    fileOut2 = open(fileNameOut_labels, 'w')
    firstTime = True
    for id, listTuple in dictNodeToListScore.iteritems():
        if firstTime:
            firstTime = False
            for i in xrange(len(listTuple)):
                fileOut.write("\tFold_" + str(i+1))  
                fileOut2.write("\tFold_" + str(i+1))  
            fileOut.write("\n")  
            fileOut2.write("\n")  
        fileOut.write(id)
        fileOut2.write(id)
        for (score, label) in listTuple:
            fileOut.write("\t" + str(score))
            fileOut2.write("\t" + str(label))
        fileOut.write("\n")
        fileOut2.write("\n")
    fileOut.close()
    fileOut2.close()
    return

if __name__ == "__main__":
    pass

