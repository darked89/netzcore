
//#include "Netscore.h"
#include "ScoreNetwork.h"
#include <unistd.h>
#include <sstream>

void scoreNetwork(string fileAffinity, string fileAbundance, int nIteration, float tError, bool flagVerbose=false, bool flagDebug=false) {
	ScoreNetwork sNetwork = ScoreNetwork(fileAffinity, fileAbundance, flagDebug, flagVerbose);
        sNetwork.run(nIteration, tError);
}

int main(int argc, char **argv)
{ 	
	static const char *optString = "p:i:t:T:R:n:e:k:K:l:L:x:y:u:v:S:h?drszaogqNwXY";
	int opt = 0;
	string fileProtein, fileInteraction;
	int type = 0, typeTransfer = IDENTITY, typeRandomization=TOPOLOGY_PRESERVING_SAME_DEGREE_INTERCHANGE, nIteration = 1, nSample = N_SAMPLE;
	int tDegreeLinkerNodeMin, tDegreeLinkerNodeMinIncrease, tDegreeLinkerEdgeMin, tDegreeLinkerEdgeMinIncrease; 
	float tError=0.01, scoreMaxAllowedNode, scoreMaxAllowedEdge, tScoreNeighborMin, tScoreNeighborMax;
	//bool flagDebug=false, flagSampleOnce = true, flagUseEdgeReliabilityScore = true, flagNormalizeOnce = false, flagAccumulate = true, flagResetSeedScoresToInitial = false, flagUpdateSampledOnce = false;
	bool flagVerbose=false, flagDebug=false, flagSampleOnce = false, flagUseEdgeReliabilityScore = false, flagNormalizeOnce = false, flagAccumulate = false, flagResetSeedScoresToInitial = false, flagUpdateSampledOnce = false, flagScaleZScoreWithNodeDegree = false, flagNetscore = false, flagSaveRandomNetworks = false, flagLoadRandomNetworks = false;
	
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
	    istringstream iss;
	    switch( opt ) {
		case 'p':
		    fileProtein = string(optarg); 
		    break;
		    
		case 'i':
		    fileInteraction = string(optarg);
		    break;
		 
		case 'n':
		    iss.str(optarg);
		    iss >> nIteration;
		    break;
		
		case 'e':
		    iss.str(optarg);
		    iss >> tError;
		    break;   

		case 't':
		    iss.str(optarg);
		    iss >> type;
		    break;
		    
		case 'T':
		    iss.str(optarg);
		    iss >> typeTransfer;
		    break;
		
		case 'R':
		    iss.str(optarg);
		    iss >> typeRandomization;
		    break;

		case 'S':
		    iss.str(optarg);
		    iss >> nSample;
		    break;

		case 'r':
		    flagUseEdgeReliabilityScore = true;
		    break;

		case 's':
		    flagSampleOnce = true;
		    break;

		case 'z':
		    flagNormalizeOnce = true;
		    flagSampleOnce = true;
		    break;

		case 'a':
		    flagAccumulate = true;
		    break;

		case 'o':
		    flagResetSeedScoresToInitial = true;
		    break;

		case 'g':
		    flagUpdateSampledOnce = true;
		    break;

		case 'q':
		    flagScaleZScoreWithNodeDegree = true;
		    break;

		case 'N':
		    flagNetscore = true;
		    break;

		case 'k':
		    iss.str(optarg);
		    iss >> tDegreeLinkerNodeMin;
		    break;

		case 'K':
		    iss.str(optarg);
		    iss >> tDegreeLinkerNodeMinIncrease;
		    break;

		case 'l':
		    iss.str(optarg);
		    iss >> tDegreeLinkerEdgeMin;
		    break;

		case 'L':
		    iss.str(optarg);
		    iss >> tDegreeLinkerEdgeMinIncrease;
		    break;

		case 'x':
		    iss.str(optarg);
		    iss >> scoreMaxAllowedNode;
		    break;

		case 'y':
		    iss.str(optarg);
		    iss >> scoreMaxAllowedEdge;
		    break;

		case 'u':
		    iss.str(optarg);
		    iss >> tScoreNeighborMin;
		    break;

		case 'v':
		    iss.str(optarg);
		    iss >> tScoreNeighborMax;
		    break;

		case 'X':
		    flagSaveRandomNetworks = true;
		    break;
			    
		case 'Y':
		    flagLoadRandomNetworks = true;
		    break;
			    
		case 'w':
		    flagVerbose = true;
		    break;
			    
		case 'd':
		    flagDebug = true;
		    break;
		case 'h':
		case '?':
		    cout << "-p <protein_file> -i <interaction_file> -t <scoring_method> " << 
		    "-T <transfer_type> -n <number_of_iterations> -e <error> " << 
		    "-k <linker_degree_node> -K <linker_degree_node_increase> " << 
		    "-l <linker_degree_edge> -L <linker_degree_edge_increase> " << 
		    "-x <node_score_max_allowed> -y <edge_score_max_allowed> " << 
		    "-u <neighbor_score_min> -v <neighbor_score_max> " << endl;
		    return 0;
		    
		default:
		    break;
	    }
	    opt = getopt( argc, argv, optString );
	}
	//inputFiles = argv + optind;
	//numInputFiles = argc - optind;
   
	/*
	//cout <<  "Jola!" << endl;
	//netscore(FILE_AFFINITY, FILE_ABUNDANCE);
	
	sscanf(argv[3], "%d", &type);
	sscanf(argv[4], "%d", &typeTransfer);
	sscanf(argv[5], "%d", &nIteration);
	sscanf(argv[6], "%f", &tError);	
	sscanf(argv[7], "%d", &tDegreeLinkerNodeMin);
	sscanf(argv[8], "%d", &tDegreeLinkerNodeMinIncrease);
	sscanf(argv[9], "%d", &tDegreeLinkerEdgeMin);
	sscanf(argv[10], "%d", &tDegreeLinkerEdgeMinIncrease);
	sscanf(argv[11], "%f", &scoreMaxAllowedNode);
	sscanf(argv[12], "%f", &scoreMaxAllowedEdge);
	sscanf(argv[13], "%f", &tScoreNeighborMin);
	sscanf(argv[14], "%f", &tScoreNeighborMax);

	//cout << " " << argv[1] << " " << argv[2] << " " << type << " " << typeTransfer << " " << nIteration << " " << tError << " " << tDegreeLinkerNodeMin << " " << tDegreeLinkerNodeMinIncrease << " " << tDegreeLinkerEdgeMin << " " << tDegreeLinkerEdgeMinIncrease << " " << scoreMaxAllowedNode << " " << scoreMaxAllowedEdge << " " << tScoreNeighborMin << " " <<  tScoreNeighborMax << endl;
	netscore(argv[1], argv[2], type, typeTransfer, nIteration, tError, tDegreeLinkerNodeMin, tDegreeLinkerNodeMinIncrease, tDegreeLinkerEdgeMin, tDegreeLinkerEdgeMinIncrease, scoreMaxAllowedNode, scoreMaxAllowedEdge, tScoreNeighborMin, tScoreNeighborMax);
	*/
	//	cout << type << " " << typeTransfer << " " << nIteration << " " << tError << " " << tDegreeLinkerNodeMin << " " << tDegreeLinkerNodeMinIncrease << " " << tDegreeLinkerEdgeMin << " " << tDegreeLinkerEdgeMinIncrease << " " << scoreMaxAllowedNode << " " << scoreMaxAllowedEdge << " " << tScoreNeighborMin << " " <<  tScoreNeighborMax << endl;
    //netscore(fileProtein, fileInteraction, type, typeTransfer, nIteration, tError, tDegreeLinkerNodeMin, tDegreeLinkerNodeMinIncrease, tDegreeLinkerEdgeMin, tDegreeLinkerEdgeMinIncrease, scoreMaxAllowedNode, scoreMaxAllowedEdge, tScoreNeighborMin, tScoreNeighborMax);
    clock_t t1 = clock();
    scoreNetwork(fileProtein, fileInteraction, 10, 0);//nIteration, tError);
    clock_t t2 = clock();
    cout << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;

    return 0;
}

