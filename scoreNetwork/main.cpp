
//#include "ScoreNetwork.hpp"
#include "Netshort.hpp"
#include "Netscore.hpp"
#include "Netrank.hpp"
#include "Netzcore.hpp"
#include "Netrandom.hpp"
#include "Netzscore.hpp"
#include "Netween.hpp"
#include "Netlink.hpp"

#include <iostream>
#include <ctime>

using namespace std;

void runScoreNetwork();
void runNetshort(string, string, string);
void runNetscore(string, string, string, unsigned int, unsigned int);
void runNetrank(string, string, string, unsigned int);
void runNetzcore(string, string, string, string, unsigned int, unsigned int, unsigned int);
void runNetrandom(string, string, string );
void runNetzscore(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration); 
void runNetz1score(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration); 
void runNetween(string fileNode, string fileEdge, string fileOutput, float threshold);
void runNetlink(string fileNode, string fileEdge, string fileOutput, float threshold);

void printHelp() {
    cout << "scoreN\n" 
	 << "-s <scoring_method>{NETSHORT:d|NETSCORE:s|NETRANK:r|NETZCORE:z|NETRANDOM:x|NETZSCORE:h|NETZ1SCORE:1|NETWEEN:w|NETLINK:l}\n"
	 << "-n <node_file>\n" 
	 << "-e <edge_file>\n" 
	 << "-o <output_file>\n"
	 << "-r <number_of_repetition>\n"
	 << "-i <number_of_iteration>\n"
	 << "-t <seed_score_threshold>\n"
	 << "-x <number_of_sampled_graphs>\n" 
	 << "-d <sampling_graph_directory>\n"
	 << endl;
}

int main(int argc, char **argv)
{
    static const char *optString = "n:e:i:s:o:r:x:d:t:h?";
    int opt = 0;
    string fileNode, fileEdge, fileOutput, dirSampling;
    enum ScoringMethod { NETSHORT = 'd', NETSCORE = 's', NETRANK = 'r', NETZCORE = 'z', NETRANDOM = 'x', NETZSCORE = 'h', NETZ1SCORE = '1', NETWEEN = 'w', NETLINK = 'l' };
    unsigned int nIteration = 1, nRepetition = 1, nSampled = 0;  
    char scoring = 's';
    float threshold = 0.01;

    if(argc < 3) {
	printHelp();
    }

    opt = getopt( argc, argv, optString );
    while( opt != -1 ) {
	istringstream iss;
        switch( opt ) {
            case 'n':
            	fileNode = string(optarg); 
                break;
            case 'e':
            	fileEdge = string(optarg);
                break;
            case 'i':
                iss.str(optarg);
            	iss >> nIteration;
                break;
            case 'r':
                iss.str(optarg);
            	iss >> nRepetition;
                break;
            case 's':
                iss.str(optarg);
                iss >> scoring;
                break;
            case 'o':
            	fileOutput = string(optarg);
                break;
	    case 'x':
                iss.str(optarg);
            	iss >> nSampled;
                break;
            case 'd':
            	dirSampling = string(optarg); 
                break;
            case 't':
                iss.str(optarg);
            	iss >> threshold;
                break;
	    case 'h':
	    case '?':
		printHelp();
                return 0;
            default:
                break;
        }
        opt = getopt( argc, argv, optString );
    }
        
    cout << "Arguments: scoring type " << scoring << ", nRepetition " << nRepetition << ", nIteration " << nIteration << ", nodeFile " << fileNode << ", edgeFile " << fileEdge << ", outputFile " << fileOutput << endl;
    clock_t t1 = clock();
    if(scoring == NETSHORT)
	runNetshort(fileNode, fileEdge, fileOutput);
    else if (scoring == NETSCORE)
	runNetscore(fileNode, fileEdge, fileOutput, nRepetition, nIteration);
    else if (scoring == NETRANK)
	runNetrank(fileNode, fileEdge, fileOutput, nIteration);
    else if (scoring == NETZCORE)
	runNetzcore(fileNode, fileEdge, fileOutput, dirSampling + string("sampled_graph.sif."), nSampled, nRepetition, nIteration); 
    else if (scoring == NETRANDOM)
	runNetrandom(fileNode, fileEdge, fileOutput); 
    else if (scoring == NETZSCORE)
	runNetzscore(fileNode, fileEdge, fileOutput, dirSampling + string("sampled_graph.sif."), nSampled, nRepetition, nIteration); 
    else if (scoring == NETZ1SCORE)
	runNetz1score(fileNode, fileEdge, fileOutput, dirSampling + string("sampled_graph.sif."), nSampled, nRepetition, nIteration); 
    else if (scoring == NETWEEN)
	runNetween(fileNode, fileEdge, fileOutput, threshold); 
    else if (scoring == NETLINK)
	runNetlink(fileNode, fileEdge, fileOutput, threshold); 
    else
	//runScoreNetwork();
	cerr << "Unrecognized scoring type!" << endl;
    clock_t t2 = clock();
    cout << "Time: " << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;
    return 0;
}

void runNetlink(string fileNode, string fileEdge, string fileOutput, float threshold) 
{
    // flagAccumulateToInitialNodeScore, flagVerbose
    Netlink sN(fileNode, fileEdge, fileOutput, threshold, true, false);
    sN.run(1, 1);
}

void runNetween(string fileNode, string fileEdge, string fileOutput, float threshold) 
{
    // seedScoreThreshold, flagAccumulateToInitialNodeScore, flagVerbose
    Netween sN(fileNode, fileEdge, fileOutput, threshold, true, false);
    sN.run();
}

void runNetz1score(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration) 
{
    // flagUseEdgeScore, flagAccumulateToInitialNodeScore, flagResetSeedScoresToInitial, flagVerbose
    Netzcore sNz(fileNode, fileEdge, fileOutput, prefixSampledGraphs, nSampled, true, true, false, false);
    Netscore sNs(fileOutput, true, true, false, false);
    sNs.setPNetwork(sNz.getPNetwork());

    // Normalize score in the begining once
    sNz.run(1, 1);
    sNs.run(nRepetition, nIteration);

    sNz.deletePNetwork();
}

void runNetzscore(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration) 
{
    // flagUseEdgeScore, flagAccumulateToInitialNodeScore, flagResetSeedScoresToInitial, flagVerbose
    //Netzscore sN(fileNode, fileEdge, fileOutput, prefixSampledGraphs, nSampled, true, false, false, true);
    Netzcore sNz(fileNode, fileEdge, fileOutput, prefixSampledGraphs, nSampled, true, false, false, false);
    Netscore sNs(fileOutput, true, true, false, false);
    sNs.setPNetwork(sNz.getPNetwork());

    // Normalize score at each repetition
    for(unsigned int r=0; r < nRepetition; ++r)
    {
	sNz.run(1, 1);
	sNs.run(1, nIteration);
    }
    sNz.deletePNetwork();
}

void runNetzcore(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration) 
{
    // flagUseEdgeScore, flagAccumulateToInitialNodeScore, flagResetSeedScoresToInitial, flagVerbose
    Netzcore sN(fileNode, fileEdge, fileOutput, prefixSampledGraphs, nSampled, true, true, false, false);
    sN.run(nRepetition, nIteration);
}

void runNetscore(string fileNode, string fileEdge, string fileOutput, unsigned int nRepetition, unsigned int nIteration) 
{
    // flagUseEdgeScore, flagAccumulateToInitialNodeScore, flagResetSeedScoresToInitial, flagVerbose
    Netscore sN(fileNode, fileEdge, fileOutput, true, true, false, false);
    sN.run(nRepetition, nIteration);
}

void runNetshort(string fileNode, string fileEdge, string fileOutput) 
{
    //clock_t t1 = clock();
    // flagAccumulateToInitialNodeScore
    Netshort sN(fileNode, fileEdge, fileOutput, true);
    //clock_t t2 = clock();
    //cout << "Time to load graph: " << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;
    sN.run();
    //cout << sN.testRun() << endl;
}

void runNetrank(string fileNode, string fileEdge, string fileOutput, unsigned int nIteration) 
{
    // flagAccumulateToInitialNodeScore
    Netrank sN(fileNode, fileEdge, fileOutput, true);
    sN.run(20); // number of pagerank iterations
}

void runNetrandom(string fileNode, string fileEdge, string fileOutput) 
{
    Netrandom sN(fileNode, fileEdge, fileOutput, true);
    sN.run();
}

/*
void runScoreNetwork() 
{
    //ScoreNetwork sN("../../data/toy_data/test_proteins_small.txt", "../../data/toy_data/test_interactions_small.txt", "../../data/out_test.txt");
    ScoreNetwork sN("../../data/toy_data/test_proteins_small.txt", "../../data/toy_data/test_interactions_small.txt", "../../data/out_test.txt", true, true, false, false);
    if(sN.getNetwork().getSize() ==  5)
	cout << "Correcto" << endl;
    else {
	cout << sN.getNetwork().getSize() << endl;
    }

    map<Vertex, float>::iterator it, itEnd;
    map<Vertex, float> mapB = sN.getNetwork().calculateBetweenness();

    for(it=mapB.begin(), itEnd=mapB.end(); it != itEnd; ++it)
	cout << sN.getNetwork().getVertexIndex(it->first) << ": " << it->second << endl;

    map<Vertex, float> mapC = sN.getNetwork().calculateClusteringCoefficient();

    for(it=mapC.begin(), itEnd=mapC.end(); it != itEnd; ++it)
	cout << sN.getNetwork().getVertexIndex(it->first) << ": " << it->second << endl;

    map<Vertex, float> mapS = sN.getNetwork().calculateShortestPath(sN.getNetwork().getVertex("v1"));

    for(it=mapS.begin(), itEnd=mapS.end(); it != itEnd; ++it)
	cout << sN.getNetwork().getVertexIndex(it->first) << ": " << it->second << endl;
}
*/

