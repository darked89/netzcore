#include "ScoreNetwork.h"
#include <fstream>
#include <cstdlib> 
#include <ctime>
#include <sstream>

ScoreNetwork::ScoreNetwork() {
    flagDebug = false;
    flagVerbose = false;
    typeTransfer = IDENTITY;
    flagUseEdgeReliabilityScore = true;
    flagAccumulate = true;
    flagResetSeedScoresToInitial = false;
    minScore = INFINITY;
    maxScore = -INFINITY;
    nIteration = 0;
    iterationCounter = 0;
}

ScoreNetwork::ScoreNetwork(string fileProtein, string fileInteraction, int typeTransfer, bool flagUseEdgeReliabilityScore, bool flagAccumulate, bool flagResetSeedScoresToInitial, bool fDebug, bool fVerbose) {
    loadProteins(fileProtein); 
    loadInteractions(fileInteraction); 
    flagDebug = fDebug;
    flagVerbose = fVerbose;
    minZScore = INFINITY;
    maxZScore = -INFINITY;
    minScore = INFINITY;
    maxScore = -INFINITY;
    nIteration = 0;
    iterationCounter = 0;
}

ScoreNetwork::~ScoreNetwork() {
}

void ScoreNetwork::run(int nIteration, float tError) { 
    nIteration = nIteration;
    float error = INFINITY;
    for(iterationCounter = 0; iterationCounter<nIteration and error > tError; iterationCounter++) {
	updateNetwork();
    }
    /*
    int i = 0;
    if(flagDebug)
        printNetwork("Initial");
    /// filter unconnected nodes & self edges
    filterNetwork(1); 
    if(flagDebug)
    	printNetwork("After filtering");
    //if(flagDebug)
    //	cout << network << endl;
    /// sample network
    //typeRandomization = DEGREE_DISTRIBUTION_PRESERVING; //DEGREE_DISTRIBUTION_PRESERVING_EDGE_INTERCHANGE; //RANDOM; //TOPOLOGY_PRESERVING_SAME_DEGREE_INTERCHANGE;
    for(iterationCounter = 0; iterationCounter<nIteration and error > tError; iterationCounter++) {
        if(flagDebug)
            cout << endl << "Iteration: " << i << endl;  
        if(flagVerbose)
            cerr << endl << "Iteration: " << i << endl;  
	initializeIteration();
	updateNetwork();
	/// calculate error & update scores
	error = calculateErrorAndUpdateScores();
	finalizeIteration();
        if(flagDebug) 
            printNetwork("After updateNetwok");
        if(flagDebug)
            cout << "error: " << error << endl;
        cerr << "i: " << i << " error: " << error << endl;
    }	
    cerr << "Iteration " << i << ":\nerror: " << error << endl;
    //if(flagDebug)
    //	printSampledNetworks();
    if(flagDebug)
        printNetwork();
    */
}

// check which is called
void ScoreNetwork::initializeIteration() {
    cout << endl << "initialize of base" << endl;  
}

//virtual void ScoreNetwork::finalizeIteration() { }

// should it be in destructor
void ScoreNetwork::finalizeScoring() {
    scaleNodeScores();
}

void ScoreNetwork::loadProteins(string const &fileName) {
    fstream file;
    file.open(fileName.c_str(), ios::in);
    string name, dummy;
    float abundance = 0.0;
    if(!file) cerr << "Warning: file can not be opened" << endl; 
    //while(file >> name >> dummy >> dummy >> abundance) {
    while(file >> name >> abundance) {
        //cout << abundance << endl;
        if(abundance > 0.0) {
            mapSource[network.addVertex(name, abundance)] = NULL;
        } else {
            network.addVertex(name, abundance);
        }
    }
    file.close();
}

void ScoreNetwork::loadInteractions(string const &fileName) {
    fstream file;
    file.open(fileName.c_str(), ios::in);
    string name1, name2;
    float affinity = 0.0;
    if(!file) cerr << "Warning: file can not be opened" << endl;
    while(file >> name1 >> name2 >> affinity) {
        //cout << affinity << endl;
        //network.addEdge(name1, name2, affinity);
        network.addEdge(name1, name2, Data(affinity)); 
    }
    file.close();
}

/*
// filter network w.r.t. given parameters 
void ScoreNetwork::filterNetwork(int degreeNodeMin, int degreeNodeMax, float scoreNodeMin, float scoreNodeMax, float scoreEdgeMin, float scoreEdgeMax, bool flagRemoveSelfEdges) {	
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
	MapIntToEdge *pMapEdge;
	Vertex *pVertex;
	Edge *pEdge;
	int vDegree; //not uInt !! because of conversion during comparison
	float vScore, eScore;
	vector<int> listNodeToRemove;
	vector< pair<int, int> > listEdgeToRemove;
	for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
		pVertex = itNode->second;		
		pMapEdge = pVertex->pMapIdToEdge;  
		vScore = pVertex->data.score;
		vDegree = pVertex->degree; //pMapEdge->size();
		#ifdef CONTROL
			if(pVertex->degree != int(pMapEdge->size())) cerr << "Warning: degree inconsistency" << endl;
			if(pVertex->id != itNode->first) cerr << "Warning: vertex id inconsistency" << endl;
		#endif //CONTROL
		if((vScore < scoreNodeMin) or (vScore > scoreNodeMax) or (vDegree < degreeNodeMin) or (vDegree > degreeNodeMax) ) {
			//network.removeVertex(itNode->first);
			listNodeToRemove.push_back(itNode->first);
		} else {
			for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
				pEdge = itEdge->second;
				eScore = pEdge->data.score;
				#ifdef CONTROL	
				if(pEdge->idSource != itNode->first or pEdge->idTarget != itEdge->first) cout << "Warning: edge vertex ids inconsistency" << endl;
				#endif //CONTROL
				if(flagRemoveSelfEdges and (pEdge->idSource == pEdge->idTarget)) {
					listEdgeToRemove.push_back( pair<int, int>(itNode->first, pEdge->idTarget) );
				} else if(eScore < scoreEdgeMin or eScore > scoreEdgeMax) {
					//network.removeEdge(pEdge->idSource, pEdge->idTarget);
					listEdgeToRemove.push_back( pair<int, int>(itNode->first, pEdge->idTarget) );
				}
			}
		}
	}
	for(unsigned int i = 0; i < listNodeToRemove.size(); i++ ) {
		network.removeVertex(listNodeToRemove[i]);
	}
	for(unsigned int i = 0; i < listEdgeToRemove.size(); i++ ) {
		network.removeEdge(listEdgeToRemove[i].first, listEdgeToRemove[i].second);
	}
        cerr << "nNode: " << network.nNode << " nEdge: "<< network.nEdge << endl;
}


void ScoreNetwork::resetScores() {
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToEdge::iterator itEdge, itEndEdge;
    MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
    MapIntToEdge *pMapEdge;
    Edge *pEdge;
    cerr << "//!CAUTION!! - initial scores may not be intial (updated after zScoring)" << endl;
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        pVertex->data.score = pVertex->data.scoreInitial; 
        pVertex->data.scoreAccumulatedLast = 0.0;
        //pVertex->data.scoreDual = 0.0;
        pMapEdge = pVertex->pMapIdToEdge;
        for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            pEdge->data.score = pEdge->data.scoreInitial;
            //pEdge->data.scoreAccumulatedLast = 0.0;
        }
    }
}
*/

void ScoreNetwork::updateNetwork() { 
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToEdge::iterator itEdge, itEndEdge;
    MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
    MapIntToEdge *pMapEdge;
    Edge *pEdge;
    hash_map<int,void*> mapIdProcessed;
	
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) ;
    /* 
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
        updateNodeScore(itNode->first); 
    }
    if(flagDebug)	
            printNetwork("After updateNodeScore");
    */
    /* //! Commented since yet not dealing with edges
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        pMapEdge = pVertex->pMapIdToEdge;
        for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
                    updateEdgeScore(typeTransfer, pEdge->idSource, pEdge->idTarget);
            }
        }
        mapIdProcessed[itNode->first] = NULL;
    }
    if(flagDebug)
        printNetwork("After updateEdgeScore");
    */
}

void ScoreNetwork::updateNodeScore(int vId) {}

void ScoreNetwork::updateEdgeScore(int idSource, int idTarget) {}

float ScoreNetwork::transferScore(float score, float a, float b) {
    switch(typeTransfer) {
    case IDENTITY:
        return score;
    case POLYNOMIAL:
        return a*score+b; // a=0.5, b=0.1, 0.2 - 0.6
    case LOGARITHMIC:
        return float(-(1/10)*log(score)); // 0.0 - 1.0 "J"
    case EXPONENTIAL:
        return float(exp(0.75*score)/10); // 0.1 - 0.2 //float(1-exp(-0.75*score)); // 0.1 - 0.5
    default:
        cout << "Warning: unidentified type" << endl;
        return 0.0;
    }
}

void ScoreNetwork::scaleNodeScores() {
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        pVertex->data.score /= maxScore; 
    }
    return;
}


float ScoreNetwork::calculateErrorAndUpdateScores() {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex, *pNeighbor;
	MapIntToEdge *pMapEdge;
    Edge *pEdge, *pEdgeNeighbor;
    hash_map<int,void*> mapIdProcessed;
    unsigned int nNode = 0, nEdge = 0;
    float error = 0.0, sumErrorNode = 0.0, sumErrorEdge = 0.0;
    minScore = INFINITY;
    maxScore = -INFINITY;
    
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        error = pVertex->data.scoreUpdated - pVertex->data.score; 
        if(isnan(pVertex->data.scoreUpdated) or isinf(pVertex->data.scoreUpdated)) {
            cerr << "Warning: scoreUpdated (and possibly error) is either nan/inf for id " << pVertex->id << " sU: " << pVertex->data.scoreUpdated << " s: " << pVertex->data.score << " e: " << error << endl;
        }
        //cout << "e: " << error << " sU: " << pVertex->data.scoreUpdated << " s: " <<  pVertex->data.score << endl;
        sumErrorNode += sq(error);
        //cout << "e^2: " << sq(error) << " sumE: " << sumErrorNode << endl;
        nNode++;        
        if(flagResetSeedScoresToInitial) {
            if(mapSource.find(pVertex->id) != mapSource.end()) { // for fFlow compatibility
                pVertex->data.scoreUpdated = pVertex->data.scoreInitial; 
            }
        }
        pVertex->data.score = pVertex->data.scoreUpdated;
        if(pVertex->data.score < minScore) {
            minScore = pVertex->data.score;
        }
        if(pVertex->data.score > maxScore) {
            maxScore = pVertex->data.score;
        }
        //pVertex->data.scoreAccumulatedLast = pVertex->data.scoreAccumulatedLastUpdated;

        /* //! Commented since yet not dealing with edges
        pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
	    	pEdge = itEdge->second;
	    	if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
                error = pEdge->data.scoreUpdated - pEdge->data.score;
	            sumErrorEdge += sq(error);
	            nEdge++;
	            pEdge->data.score = pEdge->data.scoreUpdated;	            
                    //pEdge->data.scoreAccumulatedLast = pEdge->data.scoreAccumulatedLastUpdated;
	            pNeighbor = (*pMapNode)[pEdge->idTarget];
                    pEdgeNeighbor = (*(pNeighbor->pMapIdToEdge))[pEdge->idSource];
	            pEdgeNeighbor->data.score = pEdgeNeighbor->data.scoreUpdated;
	            //pEdgeNeighbor->data.scoreAccumulatedLast = pEdgeNeighbor->data.scoreAccumulatedLastUpdated;
	    	}
        }
        */
	    mapIdProcessed[itNode->first] = NULL;
    }
    cerr << "minScore: " << minScore << " maxScore: " << maxScore << endl;
	#ifdef CONTROL
		if(nNode != network.nNode) cerr << "Warning: nNode inconsistency" << nNode << " " << network.nNode << endl;
		//if(nEdge != network.nEdge) cerr << "Warning: nEdge inconsistency" << nEdge << " " << network.nEdge << endl;
	#endif //CONTROL
    //return max(float(sqrt( (sumErrorEdge/nEdge)*(sumErrorNode/nNode) )), float( (sqrt(sumErrorNode/nNode)+ sqrt(sumErrorEdge/nEdge)) /2)); //float(sqrt(sumErrorEdge/nEdge)); //float(sqrt(sumErrorNode/nNode)+sqrt(sumErrorEdge/nEdge))/2;
    return float(sqrt( sumErrorNode/nNode ));
}


//void ScoreNetwork::outputGraph(Graph const & g, string fileName) {
void ScoreNetwork::outputGraph(Graph *g, string fileName) {
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToEdge::iterator itEdge, itEndEdge;
    MapIntToVertex *pMapNode = g->pMapIdToNode; //g.pMapIdToNode;
    MapIntToInt *pMapId = g->pMapIdToIdNew; //g.pMapIdToIdNew;
    Vertex *pVertex;
    MapIntToEdge *pMapEdge;
    Edge *pEdge;
    hash_map<int,void*> mapIdProcessed;
    MapIntToString mapName = g->mapIdToName;
    ofstream oss;
    if(fileName != "") {
            oss.open(fileName.c_str(), ios::out);
            if(!oss.is_open()) {
                cerr << "Warning: file can not be opened " << fileName << endl;
            }
    }
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	pVertex = itNode->second;
    	pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
/*
                if(format == SIF_FORMAT) {
                    if(g.isSampled()) {
                        if(fileName != "") {
                                oss << mapName[((*pMapId)[pEdge->idSource])] << "\tpp\t" << mapName[((*pMapId)[pEdge->idTarget])] << endl;	            
                        } else {
                                cout << mapName[((*pMapId)[pEdge->idSource])] << "\tpp\t" << mapName[((*pMapId)[pEdge->idTarget])] << endl;	            
                        }
                    } else {
                        if(fileName != "") {
                                oss << mapName[pEdge->idSource] << "\tpp\t" << mapName[pEdge->idTarget] << endl;
                        } else {
                                cout << mapName[pEdge->idSource] << "\tpp\t" << mapName[pEdge->idTarget] << endl;
                        }
                    }
                } else if(format == DETAILED_FORMAT) {
                    if(g.isSampled()) {
                        if(fileName != "") {
                                oss << mapName[((*pMapId)[pEdge->idSource])] << " " << (*pMapNode)[pEdge->idSource]->data.scoreUpdated 
                                             << "(" << (*pMapNode)[pEdge->idSource]->data.score << ")\t" 
                                             << mapName[((*pMapId)[pEdge->idTarget])] << " " << (*pMapNode)[pEdge->idTarget]->data.scoreUpdated 
                                             << "(" << (*pMapNode)[pEdge->idTarget]->data.score << ")\t" 
                                 << "\t" << pEdge->data.scoreUpdated << "(" << pEdge->data.score << ")" << endl;	            
                        } else {
                                cout << mapName[((*pMapId)[pEdge->idSource])] << " " << (*pMapNode)[pEdge->idSource]->data.scoreUpdated 
                                             << "(" << (*pMapNode)[pEdge->idSource]->data.score << ")\t" 
                                             << mapName[((*pMapId)[pEdge->idTarget])] << " " << (*pMapNode)[pEdge->idTarget]->data.scoreUpdated 
                                             << "(" << (*pMapNode)[pEdge->idTarget]->data.score << ")\t" 
                                 << "\t" << pEdge->data.scoreUpdated << "(" << pEdge->data.score << ")" << endl;	            
                        }
                    } else {
                        if(fileName != "") {
                                    oss << mapName[pEdge->idSource] << " " << (*pMapNode)[pEdge->idSource]->data.scoreUpdated 
                                             << "(" << (*pMapNode)[pEdge->idSource]->data.score << ")\t" 
                                             << mapName[pEdge->idTarget] << " " << (*pMapNode)[pEdge->idTarget]->data.scoreUpdated 
                                             << "(" << (*pMapNode)[pEdge->idTarget]->data.score << ")\t" 
                                     << "\t" << pEdge->data.scoreUpdated << "(" << pEdge->data.score << ")" << endl;	
                        } else {
                                    cout << mapName[pEdge->idSource] << " " << (*pMapNode)[pEdge->idSource]->data.scoreUpdated 
                                             << "(" << (*pMapNode)[pEdge->idSource]->data.score << ")\t" 
                                             << mapName[pEdge->idTarget] << " " << (*pMapNode)[pEdge->idTarget]->data.scoreUpdated 
                                             << "(" << (*pMapNode)[pEdge->idTarget]->data.score << ")\t" 
                                     << "\t" << pEdge->data.scoreUpdated << "(" << pEdge->data.score << ")" << endl;	
                        }
                    }
                } else if(format == COMPACT_FORMAT) {
*/
                    //if(g.isSampled()) {
                    if(g->isSampled()) {
                        if(fileName != "") {
                                oss << mapName[((*pMapId)[pEdge->idSource])] << "\t" << (*pMapNode)[pEdge->idSource]->data.score << "\t"
                                    << mapName[((*pMapId)[pEdge->idTarget])] << "\t" << (*pMapNode)[pEdge->idTarget]->data.score << "\t"
                                    << pEdge->data.score << endl;	            
                        } else {
                                cout << mapName[((*pMapId)[pEdge->idSource])] << "\t" << (*pMapNode)[pEdge->idSource]->data.score << "\t"
                                     << mapName[((*pMapId)[pEdge->idTarget])] << "\t" << (*pMapNode)[pEdge->idTarget]->data.score << "\t" 
                                     << pEdge->data.score << endl;	            
                        }
                    } else {
                        if( (*pMapNode)[pEdge->idSource]->data.scoreUpdated != (*pMapNode)[pEdge->idSource]->data.score ) {
                            cerr << "Warning: updated score and score inconsistency" << endl;
                        }
                        if(fileName != "") {
                                    oss << mapName[pEdge->idSource] << "\t" << (*pMapNode)[pEdge->idSource]->data.score << "\t" 
                                        << mapName[pEdge->idTarget] << "\t" << (*pMapNode)[pEdge->idTarget]->data.score << "\t" 
                                        << pEdge->data.score << endl;	
                        } else {
                                    cout << mapName[pEdge->idSource] << "\t" << (*pMapNode)[pEdge->idSource]->data.score << "\t"
                                         << mapName[pEdge->idTarget] << "\t" << (*pMapNode)[pEdge->idTarget]->data.score << "\t" 
                                         << "\t" << pEdge->data.score << endl;	
                        }
                    }
/*
                } else {
                    cerr << "Warning: unknown output format" << endl;
                }
*/
            }
        }
	    mapIdProcessed[itNode->first] = NULL;
    }
    //cout << "----" << endl;
    return;
}

void ScoreNetwork::printNetwork(string explanation) {
    printGraph(network, explanation);
}

void ScoreNetwork::printGraph(Graph const & g, string explanation) {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = g.pMapIdToNode;
	MapIntToInt *pMapId = g.pMapIdToIdNew;
    Vertex *pVertex;
	MapIntToEdge *pMapEdge;
    Edge *pEdge;
    hash_map<int,void*> mapIdProcessed;
    MapIntToString mapName = g.mapIdToName;
    if(explanation != "") cout << "----" << explanation << ": " << endl;
    //cout << g.isSampled() << endl; 
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	pVertex = itNode->second;
    	pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
            	if(g.isSampled()) {
		            cout << mapName[((*pMapId)[pEdge->idSource])] << " " << (*pMapNode)[pEdge->idSource]->data.score << "\t" 
		 	       	 << mapName[((*pMapId)[pEdge->idTarget])] << " " << (*pMapNode)[pEdge->idTarget]->data.score << "\t" 
                                 << pEdge->data.score << endl;	            
            	} else {
	        		cout << mapName[pEdge->idSource] << " " << (*pMapNode)[pEdge->idSource]->data.score << "\t"
		 	             << mapName[pEdge->idTarget] << " " << (*pMapNode)[pEdge->idTarget]->data.score << "\t" 
		                     << pEdge->data.score << endl;	
	            }
            }
        }
	    mapIdProcessed[itNode->first] = NULL;
    }
    //cout << "----" << endl;
}


