#!python 
#:vim:tabstop=4

#########################################################################
# Netzcore Utility Library (compiled from previous scripts)
# eg 9/08/2008
#
# Input: See usage
#
# Output: See usage
#
# Usage: python netzcore.py (see help & usage info printed on the screen)
#
# Creates: 3 directories named input/ & output/ & randomNetworks/
# 
#########################################################################

from sets import *
import re
import sys, os
import random, math
import networkx
import getopt

INPUT_DIR_NAME = "input" 
OUTPUT_DIR_NAME = "output"
SAMPLE_DIR_NAME = "randomNetworks"

DEFAULT_SEED_SCORE = 1.0
DEFAULT_NON_SEED_SCORE = 0.0
DEFAULT_EDGE_SCORE = 1.0

def main():
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "t:s:p:f:i:n:k:e:g:bx:S:o:awyvhZ") #, ["help", "output="])
        #print opts, args
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    verbose = False
    seed_file_name = None
    node_file_name = None
    edge_file_name = None
    n_fold = None
    n_iteration = None
    error = None
    plot_type = "ROC" 
    boost_seed = False
    t_score = None
    n_sample = None
    writeSampled = False 
    readSampled = False 
    flag_store_edge_value = False
    flag_normalize_once = True
    output_prefix = ""

    for o, a in opts:
        #print o, a
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o == "-t":
            if a in ("prepare", "run", "analyze", "runFflow", "runKMajority", "other"):
                operation_type = a
            else:
                print "Unknown operation type"
                usage()
                sys.exit(2)
        elif o == "-s":
            seed_file_name = a
        elif o == "-p":
            node_file_name = a
            print "Ignoring node file anyway"
        elif o == "-i":
            edge_file_name = a
        elif o == "-f":
            if a in ("sif", "t3t", "t5t"):
                network_format = a
                if network_format != "sif":
                    print "Coming soon: t3t and t5t formats are not supported yet.."
            else:
                print "Unknown network format"
                usage()
                sys.exit(2)
            print "Using default sif format anyway"
            network_format = "sif"
        elif o == "-n":
            n_fold = int(a)
        elif o == "-k":
            n_iteration = int(a)
        elif o == "-e":
            error = float(a)
        elif o == "-g":
            if a in ("ROC", "PPV", "AUC"):
                plot_type = a
            else:
                print "Unrecognized graphical representation type"
                usage()
                sys.exit(2)
        elif o == "-b":
            boost_seed = True
        elif o == "-x":
            t_score = float(a)
        elif o == "-S":
            n_sample = int(a)
        elif o == "-o":
            output_prefix = a
        elif o == "-w":
            writeSampled = True
        elif o == "-y":
            readSampled = True
        elif o == "-a":
            flag_store_edge_value = True
        elif o == "-Z": # set to make normalization at each step
            normalize_once = False
        else:
            assert False, "unhandled option"

    if len(opts) < 3:
        usage()
        sys.exit(2)
    else:
        f = open("chances.txt", "a")
        f2 = open("all.txt", "a")
        if random.randint(0, 41) == 17:
            wellcome()
            f.write(".")
            f2.write(".")
        else:
            f2.write(".")
        f.close()
        f2.close()

    try:
        if not os.path.exists("./"+INPUT_DIR_NAME):
            os.mkdir("./"+INPUT_DIR_NAME)
        if not os.path.exists("./"+OUTPUT_DIR_NAME):
            os.mkdir("./"+OUTPUT_DIR_NAME)
        if not os.path.exists("./"+SAMPLE_DIR_NAME):
            os.mkdir("./"+SAMPLE_DIR_NAME)
    except exc, err:
        print "Warning: Input/Output directory can not be created"
        print str(err)

    if operation_type == "prepare":
        prepareInput(fileSeed=seed_file_name, fileEdge=edge_file_name, nFold=n_fold, flagStoreEdgeValue=flag_store_edge_value)
    elif operation_type == "run":
        runNetzcore(nFold = n_fold, nIteration = n_iteration, error = error, nSample = n_sample, writeSampled = writeSampled, readSampled = readSampled, flagNormalizeOnce = flag_normalize_once, verbose = verbose)
    elif operation_type == "analyze":
        analyzePerformance(nFold = n_fold, output_prefix = output_prefix, plotType = plot_type, seedBoost = boost_seed, plotAdd = False, tScore = t_score)
    elif operation_type == "runFflow":
        runFflow(nFold = n_fold, nIteration = n_iteration)
    elif operation_type == "runKMajority":
        runKMajority(nFold = n_fold, nIteration = n_iteration)
    elif operation_type == "other":
        #analyzeNetwork(edge_file_name, output_prefix)
        analyzeNetZscoreResults(seed_file_name = seed_file_name, fileName_network = edge_file_name, output_prefix = output_prefix, tScore = t_score)

    return

def usage():
    #print sys.argv[0], "-t prepare -n <n_fold> -s <seed_nodes_file>{one_node_id_per_line/in_sif_format} -i <interaction_network_file> [-f <network_format>{(sif|t3t|t5t) -p <node_score_file>{node_id_and_node_score_each_line}]"
    print "\n\t", sys.argv[0], "-t <operation_type>{prepare|run|analyze} <list_of_option_type_specific_arguments>{explained_below}"
    print "\n------- Option Specific Arguments"
    print "\nPREPARE: -t prepare -n <n_fold> -s <seed_nodes_file>{one_node_id_per_line/in_sif_format} -i <interaction_network_file_in_sif_format> [-a <flag_store_edge_value>{False} ]"
    print "\nRUN: -t run -n <n_fold> -k <n_iteration> [ -Z{default:False(normalize_once)/normalize_at_each_step} -w{default:False/use_the_same_sampled_networks_by_writing/reading_file} -e <convergence_criteron_error>{default:0.01} -S <n_graphs_to_be_sampled>{default:10} -v{verbose_mode_shows_messaging_steps_in_error_files ]" 
    print "\nRUN FUNCTIONAL FLOW: -t runFflow -n <n_fold> -k <n_iteration>" 
    print "\nRUN K-MAJORITY: -t runKMajority -n <n_fold> -k <n_iteration>" 
    print "\nANALYZE: -t analyze -n <n_fold> -g <performance_representation>{ROC|PPV|AUC} [ -b{if_specified_assumes_maxScore_for_seed_nodes} -x <threshold_score>{prints_basic_performance_metrics_for_this_threshold_score} ]" 
    print "OTHER: -t other -s <seed_nodes_file> -i <interaction_network_file_in_sif_format> [ -o <output_prefix>{prefix_to_be_added_to_the_begining_of_output_files} -x <threshold_score>{prints_basic_performance_metrics_for_this_threshold_score} ]"
    print "\n------- Brief Description of Option Types"
    print "PREPARE: create input files in the directory named input."
    print "RUN: run NetZcore program to score nodes with the files created by PREPARE. Scored nodes (in a 5 column network format: id1 score1 id2 score2 score_edge) are saved in the directory named output."
    print "RUN FUNCTIONAL FLOW: run functionalFlow program to score nodes with the files created by PREPARE. Scored nodes (in a 5 column network format: id1 score1 id2 score2 score_edge) are saved in the directory named output."
    print "ANALYZE: create performance files for ROCR and visualization script for R in output folder and calculate performance metrics."
    print "OTHER: Some temporary tasks incorporated just to use methods defined here."
    print 
    return

def wellcome():
    print "\nDear User,\nFor customer support please call + 34 690 018 028 (about 25cents/min depending on mutual operator agreements) or email to emre <dat> guney <at> upf <dut> edu\nHave fun!"
    print "-------\n"
    return

def prepareInput(fileSeed, fileEdge, nFold, flagStoreEdgeValue=False, fileNode=None):
    """
        fileNode: node score
        fileEdge: if not sif format edge scores (t3t) or both node and edge scores (t5t)
    """
    print "Prepare: ", fileSeed, fileEdge, nFold
    if fileSeed is None or fileEdge is None or nFold is None:
        print "Prepare: missing argument"
        usage()
        sys.exit(2)

    setSeed, setSeedEdge, dictSeed, dictDummy = getNodesAndEdgesFromSifFile(fileSeed, flagStoreEdgeValue)
    setInteractor, setInteractorEdge, dictInteractor, dictEdgeWeight = getNodesAndEdgesFromSifFile(fileEdge, flagStoreEdgeValue)
    #print setInteractor
    setInteractorReal = Set()
    for i,j in setInteractorEdge:
        setInteractorReal.add(i)
        setInteractorReal.add(j)
    setSeedExisting = setSeed & setInteractorReal
    print " nodes: %d \n nodes w/ interaction: %d \n seeds: %d \n seeds w/interaction: %d" % (len(setInteractor), len(setInteractorReal), len(setSeed), len(setSeedExisting))
    
    #writeSetIntoFile("seeds_in_network.txt", setSeedExisting)

    createCrossValidationProteinFiles(namePrefix = INPUT_DIR_NAME + "/proteins", setRoot = setSeedExisting, setInteractor = setInteractorReal, nFold = nFold, dictRoot = dictSeed)
    #createPercentageValidationProteinFiles(namePrefix = INPUT_DIR_NAME + "/proteins", setRoot = setSeedExisting, setInteractor = setInteractorReal, nFold = nFold, dictRoot = dictSeed)
    writeEdgeInT3tFormat(fileName = INPUT_DIR_NAME + "/interactions.txt", setEdge = setInteractorEdge, dictEdge = dictEdgeWeight) #, setRoot = setSeedExisting)

    #(dictInteractionToSetDB, dictDBToSetInteraction, setRoot) = readPianaNetworkIntoDictionaries(filePiana)
    #print len(dictInteractionToSetDB)
    #for i in dictDBToSetInteraction.iterkeys():
    #    print i
    #for i,j in dictInteractionToSetDB.iterkeys():
    #    print i,j
    #setRootOrg = readRootProteinsIntoSet(seed_file_name)
    #f = open("root_proteins.txt", "w")
    #for i in setRootOrg:
    #    f.write(i + "\n")
    #f.close()
    #print len(setRoot), len(setRootOrg)
    #print len(setRootOrg)
    #print setRootOrg - setRoot
    #for i in setRoot:
    #    print i

    #g = createNetwork(dictInteractionToSetDB)
    #print getNetworkRadius(g)

    #setInteractor = getInteractorSet(dictInteractionToSetDB)

    #dictDBToReliability = calculateDBReliability(dictDBToSetInteraction, setRoot)
    #print dictDBToReliability

    #createCytoscapeNetworkFromInteractionDictionary(INPUT_DIR_NAME + "interactions.sif", dictInteractionToSetDB)

    #createAffinityFileFromInteractionDictionary(INPUT_DIR_NAME + "interactions.txt", dictInteractionToSetDB, dictDBToReliability)

    #for root in setRoot:
    #    if random.random() < 1.0/nFold:
    #        setRootTrain.add(root)
    #    else:
    #        print root

    #createAffinityFileFromInteractionDictionary("affinity.txt", dictInteractionToSetDB)
    #createAbundanceAndLocalityFiles("local.txt", "abundance.txt", setRootTrain, setInteractor)
    #createCytoscapeNetworkFromInteractionDictionary("affinity.sif", dictInteractionToSetDB)

    return

def runNetzcore(nFold, nIteration, error, nSample, writeSampled, readSampled, flagNormalizeOnce, verbose):
    if flagNormalizeOnce:
        strAdditionalParameters = "-g -z" # normalize once
    else:
        strAdditionalParameters = "-s" # normalize each step
    if verbose:
        strAdditionalParameters += " -w" # verbose mode
    if error is not None:
        strAdditionalParameters += " -e %f" % error
    if nSample is not None:
        strAdditionalParameters += " -S %d" % nSample

    #tError = error 
    typeRandomization = 101
    for i in xrange(1, nFold+1):
        fileProtein = INPUT_DIR_NAME + "/proteins_%d.txt" % i
        fileInteraction = INPUT_DIR_NAME + "/interactions.txt" 
        if writeSampled:
            if i == 1 and not readSampled:
                strAdditionalParametersUpdated = strAdditionalParameters + " -X"
            else:
                strAdditionalParametersUpdated = strAdditionalParameters + " -Y"
        else:
            strAdditionalParametersUpdated = strAdditionalParameters

        if nSample == 0:
            print "./netzcore -p %s -i %s -n %d -N" % (fileProtein, fileInteraction, nIteration)
            os.system("./netzcore -p %s -i %s -n %d -N > %s/proteins_%d.out 2> %s/proteins_%d.err" % (fileProtein, fileInteraction, nIteration, OUTPUT_DIR_NAME, i, OUTPUT_DIR_NAME, i))
        else:
            print "./netzcore -p %s -i %s -n %d -R %d %s" % (fileProtein, fileInteraction, nIteration, typeRandomization, strAdditionalParametersUpdated)
            os.system("./netzcore -p %s -i %s -n %d -R %d %s > %s/proteins_%d.out 2> %s/proteins_%d.err" % (fileProtein, fileInteraction, nIteration, typeRandomization, strAdditionalParametersUpdated, OUTPUT_DIR_NAME, i, OUTPUT_DIR_NAME, i))
        #print "./netzcore -p %s -i %s -n %d -e %f -R %d %s" % (fileProtein, fileInteraction, nIteration, tError, nSample, typeRandomization, strAdditionalParameters)
        #os.system("./netzcore -p %s -i %s -n %d -e %f -R %d %s > %s/proteins_%d.out 2> %s/proteins_%d.err" % (fileProtein, fileInteraction, nIteration, tError, typeRandomization, strAdditionalParameters, OUTPUT_DIR_NAME, i, OUTPUT_DIR_NAME, i))
    return

def runFflow(nFold, nIteration):
    for i in xrange(1, nFold+1):
        fileProtein = INPUT_DIR_NAME + "/proteins_%d.txt" % i
        fileInteraction = INPUT_DIR_NAME + "/interactions.txt" 
        print "./functionalFlow %s %s %d" % (fileProtein, fileInteraction, nIteration)
        os.system("./functionalFlow %s %s %d > %s/proteins_%d.out 2> %s/proteins_%d.err" % (fileProtein, fileInteraction, nIteration, OUTPUT_DIR_NAME, i, OUTPUT_DIR_NAME, i))
    for i in xrange(1, nFold+1):
        fileName_root_test = INPUT_DIR_NAME + "/proteins_test_%d.txt" % i
        fileName_network = OUTPUT_DIR_NAME + "/proteins_%d.out" % i
        fileName_out = OUTPUT_DIR_NAME + "/joan_proteins_%d.txt" % i
        setRootTest = readRootProteinsIntoSet(fileName_root_test)
        dictNodeToScore, minScore, maxScore = parseNodesScoresFromProgramOutputIntoDictionary(fileName_network, True)
        for id in setRootTest:
            dictNodeToScore[id+"*"] = dictNodeToScore[id]
            del dictNodeToScore[id]
        f = open(fileName_out, "w")
        for id, score in dictNodeToScore.iteritems():
            f.write("%s\t%f\n" % (id, score))
        f.close()
    return

def runKMajority(nFold, nIteration):
    for i in xrange(1, nFold+1):
        fileProtein = INPUT_DIR_NAME + "/proteins_%d.txt" % i
        fileInteraction = INPUT_DIR_NAME + "/interactions.txt" 
        print "./kMajority %s %s %d e" % (fileProtein, fileInteraction, nIteration)
        os.system("./kMajority %s %s %d e > %s/proteins_%d.out 2> %s/proteins_%d.err" % (fileProtein, fileInteraction, nIteration, OUTPUT_DIR_NAME, i, OUTPUT_DIR_NAME, i))
    return

def analyzePerformance(nFold, output_prefix, plotType, seedBoost, plotAdd, tScore = None):
    dictNodeToListScore = {}
    setRootAll = Set()
    for i in xrange(1, nFold+1):
        fileName_root_test = INPUT_DIR_NAME + "/proteins_test_%d.txt" % i
        setRootTest = readRootProteinsIntoSet(fileName_root_test)
        setRootAll = setRootAll.union(setRootTest)
    for i in xrange(1, nFold+1):
        fileName_root_test = INPUT_DIR_NAME + "/proteins_test_%d.txt" % i
        fileName_network = OUTPUT_DIR_NAME + "/proteins_%d.out" % i
        setRootTest = readRootProteinsIntoSet(fileName_root_test)
        setRootTrain = setRootAll-setRootTest
        dictNodeToScore, minScore, maxScore = parseNodesScoresFromProgramOutputIntoDictionary(fileName_network, True)
        #dictNodeToScore, minScore, maxScore = parseNodesScoresFromProgramOutputIntoDictionary(fileName_network, False)
        if tScore is not None:
            (acc, sens, spec, ppv) = calculatePerformance(setRootTest, dictNodeToScore, tScore, setRootTrain)
            if ppv is None: ppv = float("-inf")
            if sens is None: sens = float("-inf")
            print "Fold: %d - For score threshold %.2f:\nacc: %.2f\tsens: %.2f\tspec: %.2f\tppv: %.2f\n" % (i, tScore, acc, sens, spec, ppv)
        for id, score in dictNodeToScore.iteritems():
            if id in setRootTest:
                label = 1
            else:
                if id in setRootTrain:
                    if seedBoost: # performance check including training nodes
                        label = 1
                        score = maxScore 
                    else:
                        label = 0
                        score = 0
                        #continue # wrong since we need a label for each fold
                else:
                    label = 0
            #print i, id, (score,label)
            dictNodeToListScore.setdefault(id, []).append((score, label))
    fileNameOut_predictions = OUTPUT_DIR_NAME + "/" + output_prefix + "predictions.txt"
    fileNameOut_labels = OUTPUT_DIR_NAME + "/" + output_prefix + "labels.txt"
    createROCRPredictionsData(dictNodeToListScore, fileNameOut_predictions, fileNameOut_labels)
    createRPerformanceScript(OUTPUT_DIR_NAME + "/drawROC.R", output_prefix, plotType, plotAdd)
    return

def analyzeNetwork(fileEdge, output_prefix=""):
    setInteractor, setInteractorEdge, dictInteractor, dictEdgeWeight = getNodesAndEdgesFromSifFile(fileEdge, False)
    g = networkx.XGraph()
    g.add_edges_from(setInteractorEdge)
    histogram = getNetworkDegreeHistogram(g)
    print histogram
    return
    file_output = OUTPUT_DIR_NAME + "/" + output_prefix + "drawDistribution.R"
    createRNetworkAnalysisScript(file_output, histogram)
    return

def createRNetworkAnalysisScript(fileName, histogram):
    f = open(fileName, "w")
    f.write("h<-c(%s)\n" % ",".join([str(i) for i in histogram])) 
    f.write("hist(h, length(h))\n") 
    f.write("\nlegend(\"topright\", c(\"Degree distribution\"))\n") 
    f.close()
    return

def analyzeNetZscoreResults(seed_file_name, fileName_network, output_prefix = "", tScore = None):
    seedBoost = False
    dictNodeToScore, minScore, maxScore = parseNodesScoresFromNetZscoreOutputIntoDictionary(fileName_network)
    setRoot, dummy, dictRoot, dummy = getNodesAndEdgesFromSifFile(seed_file_name, False)
    print len(setRoot)
    print len(dictNodeToScore)
    dictNodeToListScore = {}
    if tScore is not None:
        (acc, sens, spec, ppv) = calculatePerformance(setRoot, dictNodeToScore, tScore)
        if ppv is None: ppv = float("-inf")
        print "For score threshold %.2f:\nacc: %.2f\tsens: %.2f\tspec: %.2f\tppv: %.2f\n" % (tScore, acc, sens, spec, ppv)
    for id, score in dictNodeToScore.iteritems():
        if id in setRoot:
            label = 1
            if seedBoost: # performance check including training nodes
                score = maxScore 
        else:
            label = 0
        dictNodeToListScore.setdefault(id, []).append((score, label))
    fileNameOut_predictions = OUTPUT_DIR_NAME + "/" + output_prefix + "predictions.txt"
    fileNameOut_labels = OUTPUT_DIR_NAME + "/" + output_prefix + "labels.txt"
    createROCRPredictionsData(dictNodeToListScore, fileNameOut_predictions, fileNameOut_labels)
    createRPerformanceScript(OUTPUT_DIR_NAME + "/" + output_prefix + "drawROC.R", output_prefix, "PPV", False)
    return

def calculatePerformance(setTrue, dictNodeToScore, tScore, setIgnore=None):
    (nTP, nFP, nFN, nTN) = (0.0, 0.0, 0.0, 0.0)
    for id, score in dictNodeToScore.iteritems():
        if setIgnore is not None:
            if id in setIgnore:
                continue
        if id in setTrue:
            if score >= tScore:
                nTP += 1
            else:
                nFN += 1
        else:
            if score >= tScore:
                nFP += 1
            else:
                nTN += 1

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

    #print tScore, sens, 1-spec

    #if spec is not None:
    #    return (sens, (1-spec))
    #else:
    #    return (sens, None)

    return (acc, sens, spec, ppv)

def parseNodesScoresFromProgramOutputIntoDictionary(fileName_network, flagParseNumberBeforePar):
    maxScore = float("-inf")
    minScore = float("inf")
    dictNodeToScore = {}
    file = open(fileName_network)
    line = file.readline()
    i=0
    while line:
        words = line.split()
        id1 = words[0]
        id2 = words[2]
        if not dictNodeToScore.has_key(id1):
            if flagParseNumberBeforePar:
                #temp = words[1][:words[1].find('(')]
                temp = words[1]
                score1 = float(temp)
            else:
                temp = re.search("(\(.+\))" ,words[1])
                score1 = float(temp.group()[1:-1])
            #print score1
            dictNodeToScore[id1] = score1
        if not dictNodeToScore.has_key(id2):
            if flagParseNumberBeforePar:
                #temp = words[3][:words[3].find('(')]
                temp = words[3]
                score2 = float(temp)
            else:
                temp = re.search("(\(.+\))" ,words[3])
                score2 = float(temp.group()[1:-1])
            #print score2
            dictNodeToScore[id2] = score2
        if score1 > maxScore:
            maxScore = score1
        if score1 < minScore:
            minScore = score1
        if score2 > maxScore:
            maxScore = score2
        if score2 < minScore:
            minScore = score2
        line = file.readline()
    return dictNodeToScore, minScore, maxScore

def parseNodesScoresFromNetZscoreOutputIntoDictionary(fileName_network):
    exp = re.compile("\s+\d+\s+(\w+)\s+\[(.+)\]\s+=>\s+\[(.+)\]")
    exp_with_err = re.compile("\s+\d+\s+(\w+)\s+\[(.+)\]\s+=>\s+\[(.+)\]\s+\+/\-\s+\[(.+)\]")
    maxScore = float("-inf")
    minScore = float("inf")
    dictNodeToScore = {}
    file = open(fileName_network)
    line = file.readline()
    lineCount = 1
    try:
        while line:
            if line.startswith("NODES"):
                line = file.readline()
                lineCount += 1
                while not line.startswith("EDGES"):
                    result = exp_with_err.search(line) # not using err
                    if result is None:
                        result = exp.search(line)
                    id = result.group(1)
                    score_old = float(result.group(2))
                    score_new = float(result.group(3))
                    #if not dictNodeToScore.has_key(id):
                        #print score_old, score_new
                        #dictNodeToScore[id] = score_new
                    dictNodeToScore[id] = score_new
                    if score_new > maxScore:
                        maxScore = score_new
                    if score_new < minScore:
                        minScore = score_new
                    line = file.readline()
                    lineCount += 1
            else:
                line = file.readline()
                lineCount += 1
    except:
        print lineCount, line
    return dictNodeToScore, minScore, maxScore


def createRPerformanceScript(fileName, output_prefix, type, add = False):
    f = open(fileName, "w")
    f.write("library(ROCR)\n\n") 
    f.write("v<-read.table(\"%spredictions.txt\")\n" % output_prefix) 
    #f.write("#v<-v[,1]\n") 
    f.write("l<-read.table(\"%slabels.txt\")\n" % output_prefix)
    #f.write("#l<-l[,1]\n") 
    f.write("pred<-prediction(v, l)\n")
    f.write("postscript(\"%sperformance.eps\", width = 6, height = 6, horizontal = FALSE, onefile = FALSE, paper = \"special\")\n" % output_prefix)
    if type == "ROC":
        f.write("%sperfNetZcore<-performance(pred, \"tpr\", \"fpr\")\n" % output_prefix)
        f.write("#%sperfNetZcorePPV<-performance(pred, \"ppv\")\n" % output_prefix)
        f.write("#%sperfNetZcoreSens<-performance(pred, \"sens\")\n" % output_prefix)
        f.write("\nplot(%sperfNetZcore, lwd=2, col=2, xlab=\"TPR\", ylab=\"FPR\", main=\"%s\", plotCI.col=2, avg=\"vertical\", spread.estimate=\"stddev\", show.spread.at=seq(0,1,by=0.20))\n" % (output_prefix, output_prefix)) 
        f.write("#plot(%sperfNetZcorePPV, lwd=2, col=2, ylab=\"PPV/Sens\", main=\"%s\")\n" % (output_prefix, output_prefix))
        f.write("#plot(%sperfNetZcoreSens, lwd=3, col=3, add=TRUE)\n" % output_prefix) 
        f.write("\nlegend(\"topright\", c(\"NetZcore\"), lty=c(1), col=c(2))\n") 
        f.write("#legend(\"bottomright\", c(\"PPV\", \"Sens\"), lty=c(1,1), col=c(2,3))\n") 
    elif type == "PPV":
        f.write("%sperfNetZcorePPV<-performance(pred, \"ppv\")\n" % output_prefix)
        f.write("%sperfNetZcoreSens<-performance(pred, \"sens\")\n" % output_prefix)
        f.write("#%sperfNetZcore<-performance(pred, \"tpr\", \"fpr\")\n" % output_prefix)
        f.write("\n#plot(%sperfNetZcore, lwd=2, col=2, xlab=\"TPR\", ylab=\"FPR\", main=\"%s\", plotCI.col=2, avg=\"vertical\", spread.estimate=\"stddev\", show.spread.at=seq(0,1,by=0.20))\n" %  (output_prefix, output_prefix)) 
        f.write("plot(%sperfNetZcorePPV, lwd=2, col=2, ylab=\"PPV/Sens\", main=\"%s\")\n" % (output_prefix, output_prefix))        
        f.write("plot(%sperfNetZcoreSens, lwd=3, col=3, add=TRUE)\n" % output_prefix) 
        f.write("\n#legend(\"topright\", c(\"NetZcore\"), lty=c(1), col=c(2))\n") 
        f.write("legend(\"bottomright\", c(\"PPV\", \"Sens\"), lty=c(1,1), col=c(2,3))\n") 
    elif type == "AUC":
        f.write("%sperfNetZcore<-performance(pred, \"auc\")\n" % output_prefix)
        f.write("print(%sperfNetZcore)\n" % output_prefix)
        f.close()
        return
    else:
        print "Warning: Unrecognized graphical representation type"  
    f.write("dev.off()\n")
    #f.write("\n#v<-read.table(\"_predictions.txt\")\n") 
    #f.write("#l<-read.table(\"_labels.txt\")\n")
    #f.write("#pred<-prediction(v, l)\n")
    #f.write("#perfNetZcore_<-performance(pred, \"tpr\", \"fpr\")\n\n")
    #if add:
    #    f.write("plot(perfNetZcore, lwd=2, col=2, plotCI.col=2, avg=\"vertical\", spread.estimate=\"stddev\", show.spread.at=seq(0,1,by=0.20), add=TRUE)\n") 
    #else:
    #    f.write("plot(perfNetZcore, lwd=2, col=2, plotCI.col=2, avg=\"vertical\", spread.estimate=\"stddev\", show.spread.at=seq(0,1,by=0.20))\n") 
    #    #f.write("#plot(perfNetZcore, lwd=2, col=2, plotCI.col=2, avg=\"threshold\", spread.estimate=\"stddev\", show.spread.at=seq(0,1,by=0.20), colorize=T)\n")
    #    #f.write("#plot(perfNetZcore, lwd=2, col=2, plotCI.col=2, print.cutoffs.at=seq(0,50,by=5))\n")
    #f.write("\n#plot(perfNetZcore_, lwd=2, col=2, plotCI.col=2, avg=\"vertical\", spread.estimate=\"stddev\", show.spread.at=seq(0,1,by=0.20), add=TRUE)\n") 
    #f.write("\nlegend(\"topright\", c(\"NetZcore\"), lty=c(1), col=c(2))\n") 
    f.close()
    #os.system("cd output; R CMD BATCH %s" % fileName) # does not work due to inproper directory - file naming
    #os.system("cd output")
    #print os.system("pwd")
    os.chdir("output")
    #print os.getcwd()
    #os.system("R CMD BATCH %s" % fileName) 
    os.system("R CMD BATCH %s" % "*.R") 
    return

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

def writeEdgeInT3tFormat(fileName, setEdge, dictEdge = None): #, setRoot = None):
    f = open(fileName, "w")
    for id1, id2 in setEdge:
        #if setRoot is not None:
        #    if id1 in setRoot:
        #        id1 = "*" + id1
        #    if id2 in setRoot:
        #        id2 = "*" + id2
        if dictEdge is not None:
            f.write("%s\t%s\t%f\n" % (id1, id2, float(dictEdge[(id1, id2)])))
        else:
            f.write("%s\t%s\t%f\n" % (id1, id2, DEFAULT_EDGE_SCORE))
    f.close()
    return

def getNodesAndEdgesFromSifFile(fileName, flagStoreEdgeValue):
    setNode = Set()
    setEdge = Set()
    dictNode = {}
    dictEdge = {}
    f=open(fileName)
    for line in f.readlines():
        words = line[:-1].split()
        id1 = words[0]
        setNode.add(id1)
        if len(words) == 2:
            score = float(words[1])
            dictNode[id1] = score
        elif len(words) == 3: 
            id2 = words[2]
            setNode.add(id2)
            setEdge.add((id1, id2))
            if flagStoreEdgeValue:
                dictEdge[(id1, id2)] = words[1]
    f.close()
    #return dictNode, dictEdge
    if len(dictNode) == 0:
        dictNode = None
    if len(dictEdge) == 0:
        dictEdge = None
    return setNode, setEdge, dictNode, dictEdge


def calculateDBReliability(dictDBToSetInteraction, setRoot):
    dictDBToReliability = {}
    for db, setInteraction in dictDBToSetInteraction.iteritems():
        nIncluded = 0
        for id1, id2 in setInteraction:
            if id1 in setRoot and id2 in setRoot:
                nIncluded += 1
        dictDBToReliability[db] = float(nIncluded) / len(setInteraction)

    return dictDBToReliability 

def getInteractorSet(dictInteractionToSetDB):
    setInteractor = Set()
    for id1, id2 in dictInteractionToSetDB.iterkeys():
        setInteractor.add(id1)
        setInteractor.add(id2)
    return setInteractor

def writeSetIntoFile(fileName, setRoot):
    f=open(fileName,"w")
    for i in setRoot:
        f.write(i+"\n")
    f.close()
    return

def createCrossValidationProteinFiles(namePrefix, setRoot, setInteractor, nFold, dictRoot = None):
    #k = float(len(setRoot)) / nFold
    #kCeil = int(math.ceil(k))
    #if (kCeil - k) >0.5:
    #    nElementPerFold = kCeil - 1
    #else:
    #    nElementPerFold = kCeil
    if nFold == 1:
        createProteinFile(namePrefix + "_%d.txt" % 1, setRoot, setInteractor, dictRoot)
        createCrossValidationTestProteinFile(namePrefix + "_test_%d.txt" % 1, Set())
        return
    listNElementPerFold = []
    k = len(setRoot) / nFold
    r = len(setRoot) % nFold 
    for iFold in xrange(1,nFold+1):
        nElementPerFold = k
        #print iFold
        if iFold <= r:
            nElementPerFold += 1
        #print nElementPerFold,
        listNElementPerFold.append(nElementPerFold)
    #print "Distribution of seeds in each fold: ", listNElementPerFold
    setRootCopy = setRoot.copy()
    for iFold in xrange(1,nFold+1):
        setRootTest = Set()
        if iFold == nFold:
            setRootTest = setRootCopy
        else:
            for i in xrange(listNElementPerFold[iFold]): #nElementPerFold
                iId = random.randint(0, len(setRootCopy)-1)
                id = None
                j=0
                for e in setRootCopy:
                    if j == iId:
                        id = e
                        break
                    j += 1
                setRootCopy.remove(id)
                setRootTest.add(id)
        setRootTrain = setRoot - setRootTest
        print "Fold %d - total: %d train: %d test: %d" % (iFold, len(setRoot), len(setRootTrain), len(setRootTest))
        createCrossValidationTestProteinFile(namePrefix + "_test_%d.txt" % iFold, setRootTest)
        createProteinFile(namePrefix + "_%d.txt" % iFold, setRootTrain, setInteractor, dictRoot)
    return

def createCrossValidationTestProteinFile(fileName, setRootTest):
    file = open(fileName, "w")
    for id in setRootTest:
        file.write("%s\n" % id)
    file.close()
    return

def createProteinFile(fileNameOut, setRoot, setInteractor, dictRoot = None): 
    fileOut = open(fileNameOut, 'w')
    fileOutLocal = open(fileNameOut+".local", 'w')
    fileOutAbundance = open(fileNameOut+".abundance", 'w')
    for id in setInteractor:
        fileOutLocal.write("%s\t1\n" % (id))
        if id in setRoot:
            #id = "*" + id
            if dictRoot is not None:
                score = dictRoot[id]
            else:
                score = DEFAULT_SEED_SCORE
            fileOut.write("%s\t%f\n" % (id, score))
            fileOutAbundance.write("%s\t0.0\t0.0\t%f\n" % (id, score))
        else:
            fileOut.write("%s\t%f\n" % (id, DEFAULT_NON_SEED_SCORE))
            fileOutAbundance.write("%s\t0.0\t0.0\t%f\n" % (id, DEFAULT_NON_SEED_SCORE))
    fileOut.close()
    fileOutLocal.close()
    fileOutAbundance.close()
    return

def createPercentageValidationProteinFiles(namePrefix, setRoot, setInteractor, nFold, dictRoot = None):
    listNElementPerFold = []
    k = len(setRoot) / nFold
    r = len(setRoot) % nFold 
    for iFold in xrange(1,nFold+1):
        nElementPerFold = k
        #print iFold
        if iFold <= r:
            nElementPerFold += 1
        #print nElementPerFold,
        listNElementPerFold.append(nElementPerFold)
    #print "Distribution of seeds in each fold: ", listNElementPerFold
    setRootCopy = setRoot.copy()
    setRootTrain = Set()
    for iFold in xrange(1,nFold+1):
        if iFold == nFold:
            setRootTrain = setRoot
        else:
            for i in xrange(listNElementPerFold[iFold]): #nElementPerFold
                iId = random.randint(0, len(setRootCopy)-1)
                id = None
                j=0
                for e in setRootCopy:
                    if j == iId:
                        id = e
                        break
                    j += 1
                setRootCopy.remove(id)
                #print len(setRootCopy), len(setRoot)
                setRootTrain.add(id)
        setRootTest = setRoot - setRootTrain
        print "Fold %d - total: %d train: %d test: %d" % (iFold, len(setRoot), len(setRootTrain), len(setRootTest))
        createCrossValidationTestProteinFile(namePrefix + "_test_%d.txt" % iFold, setRootTest)
        createProteinFile(namePrefix + "_%d.txt" % iFold, setRootTrain, setInteractor, dictRoot)
    return


def createNetwork(dictInteractionToSetDB):
    g = networkx.XGraph()
    for id1, id2 in dictInteractionToSetDB.iterkeys():
        g.add_edge(id1,id2)
    return g

def getNetworkRadius(g):
    return networkx.radius(g)

def getNetworkDegreeHistogram(g):
    return networkx.degree_histogram(g)

def readRootProteinsIntoSet(fileNameIn_root):  
    # read "root_proteins"
    fileIn_root = open(fileNameIn_root)
    setRoot = Set()
    line = fileIn_root.readline()
    while line:
        setRoot.add(line[:-1])
        line = fileIn_root.readline()
    fileIn_root.close()
    return setRoot

def createAffinityFileFromInteractionDictionary(fileNameOut_train, dictInteractionToSetDB, dictDBToReliability):
    dictAffinity = {}
    file = open(fileNameOut_train, 'w')
    for (tupleInteraction, setDB) in dictInteractionToSetDB.iteritems():
        id1, id2 = tupleInteraction
        affinity = 1
        for db in setDB:
            affinity *= 1-dictDBToReliability[db] # float(len(setDB))/nDB
        affinity = 1 - affinity
        file.write("%s\t%s\t%.6f\n" % (id1, id2, affinity))
    file.close()
    return

def readPianaNetworkIntoDictionaries(fileNameInPianaNetwork, flagCompact=False):    
    INT_1 = 0
    INT_2 = 2
    DB = 4
    fileIn = open(fileNameInPianaNetwork)
    dictInteractionToSetDB = {}
    dictDBToSetInteraction = {}
    setInteractor = Set()
    setDB = Set()
    setRoot = Set()
    line = fileIn.readline()
    line = fileIn.readline()
    while line:
        id1 = None
        id2 = None
        if flagCompact:
            words = line.split()
            id1 = words[INT_1]
            id2 = words[INT_2]
            listDB = [words[DB]]
        else:
            listDB = []
            iterator = re.finditer("db=\w*|protein_1=\w*|protein_2=\w*|root_1=is-root|root_2=is-root", line)
            for pattern in iterator:
                pattern = pattern.group()
                if pattern.startswith("db="):
                    listDB.append(pattern[3:])
                elif pattern.startswith("protein_1="):
                    id1 = pattern[10:]
                elif pattern.startswith("protein_2="):
                    id2 = pattern[10:]
                elif pattern.startswith("root_1="):
                    setRoot.add(id1)
                elif pattern.startswith("root_2="):
                    setRoot.add(id2)
        if id1 is None or id2 is None:
            line = fileIn.readline()
            continue
        setInteractor.add(id1)
        setInteractor.add(id2)
        for db in listDB:
            setDB.add(db)
        if dictInteractionToSetDB.has_key((id1, id2)):
            for db in listDB:
                dictInteractionToSetDB[(id1, id2)].add(db)
        elif dictInteractionToSetDB.has_key((id2, id1)):
            for db in listDB:
                dictInteractionToSetDB[(id2, id1)].add(db)  
        else:
            dictInteractionToSetDB[(id1, id2)]=Set(listDB)
        for db in listDB:
            if dictDBToSetInteraction.has_key(db):
                if not ((id2, id1) in dictDBToSetInteraction[db]):
                    dictDBToSetInteraction[db].add((id1, id2)) 
            else:
                dictDBToSetInteraction[db] = Set([(id1, id2)])
        line = fileIn.readline()
    
    #print setDB
    nInteraction = 0
    nDB = 0.0 # len(setDB) 
    union = Set()
    for db, setInteraction in dictDBToSetInteraction.iteritems():
        nInteraction += len(setInteraction)
        nDB += 1
        print db, len(setInteraction)
        union &= setInteraction
    print len(union) 

    return (dictInteractionToSetDB, dictDBToSetInteraction, setRoot)


if __name__ == "__main__":
    main()

