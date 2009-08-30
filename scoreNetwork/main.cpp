
//#include "ScoreNetwork.hpp"
#include "Netshort.hpp"
#include "Netscore.hpp"
#include "Netrank.hpp"
#include "Netzcore.hpp"

#include <iostream>
#include <ctime>

using namespace std;

void runScoreNetwork();
void runNetshort(string, string, string);
void runNetscore(string, string, string, unsigned int, unsigned int);
void runNetrank(string, string, string);
void runNetzcore(string, string, string, string, unsigned int, unsigned int, unsigned int);

int main(int argc, char **argv)
{
    static const char *optString = "n:e:i:s:o:r:x:d:h?";
    int opt = 0;
    string fileNode, fileEdge, fileOutput, dirSampling;
    enum ScoringMethod { NETSHORT = 'd', NETSCORE = 's', NETRANK = 'r', NETZCORE = 'z' };
    unsigned int nIteration = 1, nRepetition = 1, nSampled = 0;  
    char scoring = 's';

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
	    case 'h':
	    case '?':
		cout << "-n <node_file> -e <edge_file> -s <scoring_method>{NETSHORT:d|NETSCORE:s|NETRANK:r|NETZCORE:z}" << 
                " -i <number_of_iteration> -o <output_file> -r <number_of_repetition> -x <number_of_sampled_graphs>" <<
	        " -d <sampling_graph_directory> " << endl;
                return 0;
            default:
                break;
        }
        opt = getopt( argc, argv, optString );
    }
        
    cout << "Arguments: scoring type " << scoring << ", nIteration " << nIteration << ", nodeFile " << fileNode << ", edgeFile " << fileEdge << ", outputFile " << fileOutput << endl;
    clock_t t1 = clock();
    if(scoring == NETSHORT)
	runNetshort(fileNode, fileEdge, fileOutput);
    else if (scoring == NETSCORE)
	runNetscore(fileNode, fileEdge, fileOutput, nRepetition, nIteration);
    else if (scoring == NETRANK)
	runNetrank(fileNode, fileEdge, fileOutput);
    else if (scoring == NETZCORE)
	runNetzcore(fileNode, fileEdge, fileOutput, dirSampling + string("sampled_graph.sif."), nSampled, nRepetition, nIteration); 
    else
	//runScoreNetwork();
	cerr << "Unrecognized scoring type!" << endl;
    clock_t t2 = clock();
    cout << "Time: " << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;
    return 0;
}

void runNetzcore(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration) 
{
    // flagUseEdgeScore, flagAccumulateToInitialNodeScore, flagResetSeedScoresToInitial, flagVerbose
    Netzcore sN(fileNode, fileEdge, fileOutput, prefixSampledGraphs, nSampled, true, false, false, false);
    sN.run(nRepetition, nIteration);
}

void runNetscore(string fileNode, string fileEdge, string fileOutput, unsigned int nRepetition, unsigned int nIteration) 
{
    // flagUseEdgeScore, flagAccumulateToInitialNodeScore, flagResetSeedScoresToInitial, flagVerbose
    Netscore sN(fileNode, fileEdge, fileOutput, true, false, false, false);
    sN.run(nRepetition, nIteration);
}

void runNetshort(string fileNode, string fileEdge, string fileOutput) 
{
    //clock_t t1 = clock();
    Netshort sN(fileNode, fileEdge, fileOutput);
    //clock_t t2 = clock();
    //cout << "Time to load graph: " << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;
    sN.run();
    //cout << sN.testRun() << endl;
}

void runNetrank(string fileNode, string fileEdge, string fileOutput) 
{
    Netrank sN(fileNode, fileEdge, fileOutput);
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

