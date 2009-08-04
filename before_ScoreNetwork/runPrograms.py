
from sets import *
import re, os

INFINITY = 9999999
type="fFlow" 
#type="kMajority"
#type="netZcore"
#fileAffinity = "/home/emre/arastirma/netscore/new_data/yeast_affinity_biggest_sub_netscore_data/yeast_affinity_biggest_sub_mod.dat"
#fileAbundance = "/home/emre/arastirma/netscore/new_data/yeast_affinity_biggest_sub_netscore_data/yeast_affinity_copy_scr.dat"
#dirBaseData = "/home/emre/arastirma/netscore/aneurism_data/"
dirBaseData = "/home/emre/arastirma/netzcore/src/input/"
dirBaseProgram = "/home/emre/workspace/"
listIteration = [5] #[20] #,3,4,5,6,20] #[1,2,4,5,20] #[20,2,4,5] #range(1,5)
tError = 0.001
listTypeRandomization = [101] #100,101,104] #[100, 101, 104] #[102] #[101, 100, 102, 104] 
#listStrAdditionalParameters = ["-s -a", "-s -g"]
nFold = 10
#fileName_root = dirBaseData + "root_proteins.txt" 
#fileName_root = dirBaseData + "../../data/aneurism_data/root_proteins.txt" 
#fileName_root = dirBaseData + "../../data/aneurism_data/uri_new/aneurysm_seeds.txt" 
fileName_root = dirBaseData + "../../data/jordi_mestres/biana_updated/roots_geneid_complete_human_network.sif" 

if type == "netZcore":
    namePrefix="netscore2"
if type == "kMajority":
    namePrefix="kMajority"
if type == "fFlow":
    namePrefix="functionalFlow"

def main():
    fileInteraction = dirBaseData + "interactions.txt"
    setRootAll = readRootProteinsIntoSet(fileName_root)
    for nIteration in listIteration:
        for typeRandomization in listTypeRandomization:
            strAdditionalParameters = "-z" # -r -q"
            for i in xrange(1, nFold+1):
                fileProtein = dirBaseData + "proteins_%d.txt" % i
                print namePrefix + "/main -p %s -i %s -n %d -e %f -R %d %s" % (fileProtein, fileInteraction, nIteration, tError, typeRandomization, strAdditionalParameters)
                #netZcore
                if type == "netZcore":
                   os.system(dirBaseProgram + namePrefix + "/main -p %s -i %s -n %d -e %f -R %d %s > %s/proteins_%d.out 2> %s/proteins_%d.err" % (fileProtein, fileInteraction, nIteration, tError, typeRandomization, strAdditionalParameters, dirBaseProgram + namePrefix, i, dirBaseProgram + namePrefix, i))
                #kMajority 
                if type == "kMajority":
                   os.system(dirBaseProgram + namePrefix + "/main %s %s %d e > %s/proteins_%d.out 2> %s/proteins_%d.err" % (fileProtein, fileInteraction, nIteration, dirBaseProgram + namePrefix, i, dirBaseProgram + namePrefix, i))
                #functionalFlow
                if type == "fFlow":
                   os.system(dirBaseProgram + namePrefix + "/main %s %s %d > %s/proteins_%d.out 2> %s/proteins_%d.err" % (fileProtein, fileInteraction, nIteration, dirBaseProgram + namePrefix, i, dirBaseProgram + namePrefix, i))
            dictNodeToListScore = {}
            for i in xrange(1, nFold+1):
                fileName_root_test = dirBaseData + "proteins_test_%d.txt" % i
                fileName_network = dirBaseProgram + namePrefix + "/proteins_%d.out" % i
                setRootTest = readRootProteinsIntoSet(fileName_root_test)
                setRootTrain = setRootAll-setRootTest
                if type == "netZcore":
                    #dictNodeToScore = parseNodesScoresFromProgramOutputIntoDictionary(fileName_network, False)
                    dictNodeToScore = parseNodesScoresFromProgramOutputIntoDictionary(fileName_network, True)
                if type == "kMajority":
                    dictNodeToScore = parseNodesScoresFromProgramOutputIntoDictionary(fileName_network, True)
                if type == "fFlow":
                    #dictNodeToScore = parseNodesScoresFromProgramOutputIntoDictionary(fileName_network, False)
                    dictNodeToScore = parseNodesScoresFromProgramOutputIntoDictionary(fileName_network, True)
                for id, score in dictNodeToScore.iteritems():
                    if id in setRootTrain:
                        pass
                        #label = 1
                        #score = INFINITY # to make things even with functionalFlow
                    if id in setRootTest:
                        label = 1
                    else:
                        label = 0
                    #print i, id, (score,label)
                    insertIntoDictionary(dictNodeToListScore, id, (score, label), "list")
            if type == "netZcore":
                fileNameOut_predictions = namePrefix + "_predictions_n%d_e%s_R%d_%s.txt.nobias" % (nIteration, repr(tError), typeRandomization, strAdditionalParameters.replace('-', '').replace(' ', '_'))
                fileNameOut_labels = namePrefix + "_labels_n%d_e%s_R%d_%s.txt.nobias" % (nIteration, repr(tError), typeRandomization, strAdditionalParameters.replace('-', '').replace(' ', '_'))
            if type == "kMajority":
                fileNameOut_predictions = namePrefix + "_predictions_n%d_e%s_e.txt" % (nIteration, repr(tError))
                fileNameOut_labels = namePrefix + "_labels_n%d_e%s_e.txt" % (nIteration, repr(tError))
            if type == "fFlow":
                fileNameOut_predictions = namePrefix + "_predictions_n%d_e%s.txt" % (nIteration, repr(tError))
                fileNameOut_labels = namePrefix + "_labels_n%d_e%s.txt" % (nIteration, repr(tError))
            createROCRPredictionsData(dictNodeToListScore, fileNameOut_predictions, fileNameOut_labels)

    return        

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

def parseNodesScoresFromProgramOutputIntoDictionary(fileName_network, flagParseNumberBeforePar):
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
        line = file.readline()
    return dictNodeToScore

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

def insertIntoDictionary(dict, key, value, type = None):
    """
    Inserts a key to a dictionary where values could be containers such as list or set
    Type could be one of: None, "list", "set"
    """
    if dict.has_key(key):
        #obj = dict[key]
        if type == "set": #hasattr(obj, "add"): #isinstance(obj, set):
            dict[key].add(value)
        elif type == "list": #hasattr(obj, "append"): #isinstance(obj, list):
            dict[key].append(value)
        else:
            dict[key] = value
            print "Warning: Key exists"
    else:
        if type == "set": #hasattr(obj, "add"): 
            dict[key] = set.Set( [value] )
        elif type == "list": #hasattr(obj, "append"):
            dict[key] = [ value ]
        else:
            dict[key] = value
    return 



main()

