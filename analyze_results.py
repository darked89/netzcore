#########################################################################
#output_scores_file+"."+association_scores_file_identifier_type Methods for analyzing output files from scoring methods
#
# eg 13/08/2009
#########################################################################

from biana.utilities import graph_utilities as network_utilities
from funcassociate import client
from biana.utilities import TsvReader


def check_functional_enrichment_of_high_scoring_modules(network_file, module_detection_type, output_scores_file, node_scores_file, node_mapping_file, percentage, id_type, output_method, default_non_seed_score = 0, exclude_seeds = False, specie = "Homo sapiens", mode = "unordered"):
    """
	Check functional enrichment of highest scoring modules at given percentage
    """
    g = network_utilities.create_network_from_sif_file(network_file)

    dummy, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = output_scores_file, store_edge_type = False)
    ids, n, i = get_top_scoring_nodes_at_given_percentage(node_to_score, percentage)

    if module_detection_type != "greedy":
	raise ValueError("Only greedy all highest neighbor node inclusion is supported!")

    sub_graph = network_utilities.get_subgraph(g, ids)
    modules = network_utilities.get_connected_components(sub_graph, return_as_graph_list=True)
    #print "NetworkX way:"
    print len(modules), map(len, modules)

    #modules =  get_high_scoring_modules(g, node_to_score, ids)
    #print "Handcrafted way:"
    #print len(modules), map(len, modules)

    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)

    all_ids = reduce(lambda x,y: x+y, id_to_mapped_ids.values())
    for module in modules:
	if len(module) == 1:
	    continue
	if exclude_seeds:
	    ids_in_module = [ id for id in module if node_to_score[id] <= default_non_seed_score ]
	else:
	    ids_in_module = module
	selected_ids = reduce(lambda x,y: x+y, [ id_to_mapped_ids[id] for id in ids_in_module])
	output_method("Module: %d genes among %d\n" % (len(selected_ids), len(all_ids)))
	output_method("%s\n" % ", ".join(selected_ids))
	check_functional_enrichment(selected_ids, all_ids, id_type, output_method, specie = specie, mode = mode)

    return


def get_high_scoring_modules(g, node_to_score, ids):
    """
	Can use ids instead of node_to_score & min_score, helper function can be removed and made inline 
    """

    def get_high_scoring_neighbors(v, g, node_to_score, min_score):
	neighbors = set()
	for u in g.neighbors(v):
	    if node_to_score[u] >= min_score:
		neighbors.add(u)
	return neighbors

    min_score = node_to_score[ids[-1]]
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
	for u in current_module:
	    if u in ids_included_in_modules:
		continue
	    neighbors = get_high_scoring_neighbors(u, g, node_to_score, min_score)
	    ids_included_in_modules.add(u) 
	    #! buggy below modifying set being iterated
	    current_module |= neighbors
    return modules


def get_id_to_mapped_id_mapping(node_mapping_file):
    reader = TsvReader.TsvReader(node_mapping_file, inner_delim = ",")
    columns, id_to_mapped_ids = reader.read(fields_to_include = None, merge_inner_values = True)
   
    id_to_mapped_ids_formatted = {}
    for id, vals in id_to_mapped_ids.iteritems():
	vals = reduce(lambda x,y: x+y, vals)
	id_to_mapped_ids_formatted[id] = vals

    return id_to_mapped_ids_formatted


def get_top_scoring_node_ids_at_given_percentage(output_scores_file, node_scores_file, node_mapping_file, percentage, id_type, default_non_seed_score=0, exclude_seeds=False):
    """
	Get ids of highest scoring nodes at given percentage
    """
    dummy, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = output_scores_file, store_edge_type = False)
    ids, n, i = get_top_scoring_nodes_at_given_percentage(node_to_score, percentage)

    #reader = TsvReader.TsvReader(node_mapping_file, inner_delim = ",")
    #columns, id_to_mapped_ids = reader.read(fields_to_include = None, merge_inner_values = True)
    id_to_mapped_ids = get_id_to_mapped_id_mapping(node_mapping_file)
   
    selected_ids = []
    all_ids = []

    if exclude_seeds:
	dummy, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = node_scores_file, store_edge_type = False)

    for id in ids:
	if exclude_seeds:
	    if node_to_score[id] > default_non_seed_score:
		continue
	#vals = reduce(lambda x,y: x+y, id_to_mapped_ids[id])
	vals = id_to_mapped_ids[id]
	selected_ids.extend(vals)

    for val_list in id_to_mapped_ids.values():
    #	all_ids.extend(reduce(lambda x,y: x+y, val_list))
    	all_ids.extend(val_list)

    return selected_ids, all_ids


def check_functional_enrichment_at_given_percentage(output_scores_file, node_scores_file, node_mapping_file, percentage, id_type, output_method, default_non_seed_score=0, exclude_seeds=False, specie = "Home sapiens", mode="unordered"):
    """
	Check functional enrichment of highest scoring nodes at given percentage
    """
    selected_ids, all_ids = get_top_scoring_node_ids_at_given_percentage(output_scores_file, node_scores_file, node_mapping_file, percentage, id_type, default_non_seed_score, exclude_seeds)
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
    #f = open(output_file, "a")
    #output_method = f.write

    if id_type == "geneid":
	id_type = "entrezgene"
    elif id_type == "genesymbol":
	id_type = "hgnc_symbol"
    elif id_type == "uniprotaccession":
	id_type = "uniprot_accession"
    elif id_type == "uniprotentry":
	id_type = "uniprot_id"
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

    #output_method("OVERREPRESENTED ATTRIBUTES"\n)

    headers = ["N", "M", "X", "LOD", "P", "P_adj", "attrib ID", "attrib name"]
    if mode == "unordered":
	headers.pop(1)
    if tex_format:
	output_method("%s\\\\\n" % " & ".join(headers))
    else:
	output_method("%s\n" % "\t".join(headers))

    zero = "< %f" % (1.0/float(reps))

    for row in response["over"]:
	if mode == "unordered":
	    row.pop(1)
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

    #f.close()
    return


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

def create_R_script(fileName, absolute_dir, title):
    f = open(fileName, "w")
    f.write("library(ROCR)\n") 
    f.write("f<-function(perf) { d<-rep(0, times=length(perf@x.values[[1]]))\n")
    f.write("for(vals in perf@x.values) { d<-d+vals; }\n")
    f.write("for(i in 1:length(d)) { if(d[i] == Inf) { d[i]<-0; }; }\n")
    f.write("d<-d/length(perf@x.values); d<-max(d); return(d); }\n")
    f.write("extractvals<-function(perf) { x=c(); y=c(); for ( i in 1:length(perf@y.values) ) { x[i]<-perf@x.values[i]; y[i]<-perf@y.values[i] }; list(x,y); }\n")
    f.write("v<-read.table(\"%spredictions.txt\")\n" % absolute_dir) 
    f.write("l<-read.table(\"%slabels.txt\")\n" % absolute_dir)
    f.write("pred<-prediction(v, l)\n")
    #if image_type = "eps":
    f.write("postscript(\"%sperformance.eps\", width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = \"special\", title = \"%s\")\n" % (absolute_dir, title))
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
    setCandidates, setDummy, dictCandidates, dictDummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = candidates_file, store_edge_type = False)
    node_validation_data = [ (dictNodeResult[id], 1) for id in setNodeTest ]
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


def calculate_performance_metric_counts(dictNodeResult, setNodeTest, non_seeds, score_threshold, n_random_negative_folds = None, replicable=123):
    (nTP, nFP, nFN, nTN) = (0.0, 0.0, 0.0, 0.0)
    for id, score in dictNodeResult.iteritems():
        if id in setNodeTest:
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
    return calculate_performance_metric_counts(dictNodeResult, setNodeTest, non_seeds, score_threshold, n_random_negative_folds, replicable)


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


def get_top_scoring_nodes_at_given_percentage(node_to_score, percentage):
    """
	Returns highest scoring nodes at given percentage
    """

    result_scores = node_to_score.items()
    result_scores.sort(lambda x,y: cmp(y[1], x[1]))

    i = 0
    n = len(result_scores)*percentage/100
    ids = [] # the order is important
    
    last_score = result_scores[0][1]
    for id, score in result_scores:
	i+=1
	if i>n:
	    if last_score != score:
		break
	ids.append(id)
	last_score = score

    return ids, n, i


def calculate_seed_coverage_at_given_percentage(file_result, file_seed_scores, percentage, default_score=0):
    """
	Calculates number of seed nodes included in the given percentage of high ranking nodes in result file 
	Returns a tuple containing number of seed nodes and all nodes at that percentage
    """

    n_seed = 0

    dummy, dummy, node_to_score_initial, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_seed_scores, store_edge_type = False)
    # getting seeds (using initial non-seed score assumption)
    seeds = set([ id for id, score in node_to_score_initial.iteritems() if score > default_score ])

    dummy, dummy, node_to_score, dummy = network_utilities.get_nodes_and_edges_from_sif_file(file_name = file_result, store_edge_type = False)
    ids, n, i = get_top_scoring_nodes_at_given_percentage(node_to_score, percentage)
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


