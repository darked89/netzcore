
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
void runNetscore(string, string, string, unsigned int);
void runNetrank(string, string, string);
void runNetzcore(string, string, string, string, unsigned int);

int main(int argc, char **argv)
{
    static const char *optString = "n:e:i:s:o:h?";
    int opt = 0;
    string fileNode, fileEdge, fileOutput;
    enum ScoringMethod { NETSHORT = 'd', NETSCORE = 's', NETRANK = 'r', NETZCORE = 'z' };
    //ScoringMethod scoring = NETSCORE;
    unsigned int nIteration = 1; //, scoring = 1; 
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
            case 's':
                iss.str(optarg);
                iss >> scoring;
                break;
            case 'o':
            	fileOutput = string(optarg);
                break;
	    case 'h':
	    case '?':
		cout << "-n <node_file> -e <edge_file> -s <scoring_method>{NETSHORT:0|NETSCORE:1|NETRANK:2}" << 
                "-i <number_of_iteration> -o <output_file>" << endl;
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
	runNetscore(fileNode, fileEdge, fileOutput, nIteration);
    else if (scoring == NETRANK)
	runNetrank(fileNode, fileEdge, fileOutput);
    else if (scoring == NETZCORE)
	runNetzcore(fileNode, fileEdge, fileOutput, "../../data/sampled_graphs/sampled_graph.txt.", 100);
    else
	//runScoreNetwork();
	cerr << "Unrecognized scoring type!" << endl;
    clock_t t2 = clock();
    cout << "Time: " << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;
    return 0;
}

void runNetzcore(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled) //, unsigned int nIteration) 
{
    Netzcore sN(fileNode, fileEdge, fileOutput, prefixSampledGraphs, nSampled, false, false, false, false);
    sN.run(1,1);
}

void runNetscore(string fileNode, string fileEdge, string fileOutput, unsigned int nIteration) 
{
    //clock_t t1 = clock();
    //Netscore sN("../../data/input/node_scores.txt", "../../data/input/edge_weights.txt");
    //Netscore sN("../../data/toy_data/test_proteins_small.txt", "../../data/toy_data/test_interactions_small.txt", true, true, false, true);
    Netscore sN(fileNode, fileEdge, fileOutput, false, false, false, false);
    //clock_t t2 = clock();
    //cout << "Time to load graph: " << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;
    sN.run(1,nIteration);
    //cout << sN.run(1,1) << endl;
}

void runNetshort(string fileNode, string fileEdge, string fileOutput) 
{
    //clock_t t1 = clock();
    //Netshort sN("../../data/input/node_scores.txt", "../../data/input/edge_weights.txt");
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

