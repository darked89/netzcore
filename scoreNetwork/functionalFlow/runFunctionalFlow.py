
from sets import *
import re, os

dirBase = "/home/emre/arastirma/netscore/aneurism_data/"
listIteration = [4] #range(1,5)
listTScore = range(-1,40, 1)
nFold = 10
fileName_root = dirBase + "root_proteins.txt" 

def main():
    fileInteraction = dirBase + "interactions.txt"
    setRootAll = readRootProteinsIntoSet(fileName_root)
    for nIteration in listIteration:
        dictTScoreToPerformanceSum = {}
        for i in xrange(1, nFold+1):
            fileProtein = dirBase + "proteins_%d.txt" % i
            print "./main %s %s %d" % (fileProtein, fileInteraction, nIteration)
            os.system("./main %s %s %d > proteins_%d.out 2> proteins_%d.err" % (fileProtein, fileInteraction, nIteration, i, i))
            fileName_root_test = dirBase + "proteins_test_%d.txt" % i
            fileName_network = "proteins_%d.out" % i
            setRootTest = readRootProteinsIntoSet(fileName_root_test)
            dictNodeToScore = parseNodesScoresFromProgramOutputIntoDictionary(fileName_network)
            for tScore in listTScore:
                if i==1:
                    dictTScoreToPerformanceSum[tScore] = (0, 0)
                (sens, OneMinusSpec) = calculatePerformance(setRootTest, dictNodeToScore, tScore, (setRootAll-setRootTest))
                (sensSum, OneMinusSpecSum) = dictTScoreToPerformanceSum[tScore]
                dictTScoreToPerformanceSum[tScore] = (sensSum + sens, OneMinusSpecSum + OneMinusSpec)

        for tScore, tuple in dictTScoreToPerformanceSum.iteritems():
            (sens, OneMinusSpec) = tuple
            print tScore, sens / nFold, OneMinusSpec / nFold
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

def parseNodesScoresFromProgramOutputIntoDictionary(fileName_network):
    dictNodeToScore = {}
    file = open(fileName_network)
    line = file.readline()
    i=0
    while line:
        words = line.split()
        id1 = words[0]
        id2 = words[2]
        if not dictNodeToScore.has_key(id1):
            temp = re.search("(\(.+\))" ,words[1])
            score1 = float(temp.group()[1:-1])
            #print score1
            dictNodeToScore[id1] = score1
        if not dictNodeToScore.has_key(id2):
            temp = re.search("(\(.+\))" ,words[3])
            score2 = float(temp.group()[1:-1])
            #print score2
            dictNodeToScore[id2] = score2
        line = file.readline()
    return dictNodeToScore

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

    sens = nTP / (nTP + nFN)
    spec = nTN / (nTN + nFP)
    ppv = nTP / (nTP + nFP)

    #print tScore, sens, 1-spec

    return (sens, (1-spec))


main()

