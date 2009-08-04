#include "Netscore.h"
#include <fstream>
#include <cstdlib> 
#include <ctime>

Netscore::Netscore() {
}

Netscore::Netscore(string fileProtein, string fileInteraction, bool fDebug) {
	loadProteins(fileProtein); 
    loadInteractions(fileInteraction); 
    flagDebug = fDebug;
    if(flagDebug)
    	printNodes("Initial: ");
}

Netscore::~Netscore() {
}

void Netscore::loadProteins(string const &fileName) {
    fstream file;
    file.open(fileName.c_str(), ios::in);
    string name, dummy;
    float abundance = 0.0;
    if(!file) cout << "Warning: file can not be opened" << endl; 
    //while(file >> name >> dummy >> dummy >> abundance) {
    while(file >> name >> abundance) {
        //cout << abundance << endl;
    	network.addVertex(name, Data(abundance));
    }
    file.close();
}

void Netscore::loadInteractions(string const &fileName) {
    fstream file;
    file.open(fileName.c_str(), ios::in);
    string name1, name2;
    float affinity = 0.0;
    if(!file) cout << "Warning: file can not be opened" << endl;
    while(file >> name1 >> name2 >> affinity) {
        //cout << affinity << endl;
        network.addEdge(name1, name2, Data(affinity));
    }
    file.close();
}

void Netscore::run(int type, int typeTransfer, int nIteration, float tError, int tDegreeLinkerNodeMin, int tDegreeLinkerNodeMinIncrease, int tDegreeLinkerEdgeMin, int tDegreeLinkerEdgeMinIncrease, float scoreMaxAllowedNode, float scoreMaxAllowedEdge, float tScoreNeighborMin, float tScoreNeighborMax) {
	float error = INFINITY;
    // specifying criteria for inputing edges based on their scores
    //filterNetwork(-INT_INFINITY, INT_INFINITY, -INFINITY, INFINITY, 0.1, INFINITY);
	filterNetwork(1);
	if(flagDebug) {
		printNodes("Filtered: ");
		printNetwork("Before normalization");
	}
	
	scaleScores(); //(true); //(5.0, 5.0);
	
    if(flagDebug)
    	printNetwork("Before iteration");
	for(int i=0; i<nIteration and error > tError; i++) {
		if(flagDebug)
			cout << endl << "----Iteration: " << i << endl;
	    // specifying no criteria for outputing edges based on their scores  
	    updateNetwork(type, typeTransfer, tDegreeLinkerNodeMin, tDegreeLinkerNodeMinIncrease, tDegreeLinkerEdgeMin, tDegreeLinkerEdgeMinIncrease, scoreMaxAllowedNode, scoreMaxAllowedEdge, tScoreNeighborMin, tScoreNeighborMax);
	    if(flagDebug)
	    	printNetwork("After updateNetwok");
        error = calculateError(type);
        if(flagDebug)
        	cout << "error: " << error << endl;
	}		

	float scoreTemp = getEdgeScore("YFL039C", "YKL007W");
	if(scoreTemp != scoreMaxAllowedEdge) {
		cout << type << " " << typeTransfer << " " << nIteration << " " << tError << " " << tDegreeLinkerNodeMin << " " << tDegreeLinkerNodeMinIncrease << " " << tDegreeLinkerEdgeMin << " " << tDegreeLinkerEdgeMinIncrease << " " << scoreMaxAllowedNode << " " << scoreMaxAllowedEdge << " " << tScoreNeighborMin << " " <<  tScoreNeighborMax << endl;
		cout << calculateScoreDifference() << " " << getEdgeScore("YLR347C", "YOR212W") << " " << scoreTemp << endl;
	}
}

// filter network w.r.t. given parameters 
void Netscore::filterNetwork(int degreeNodeMin, int degreeNodeMax, float scoreNodeMin, float scoreNodeMax, float scoreEdgeMin, float scoreEdgeMax) {	
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
			if(pVertex->degree != int(pMapEdge->size())) cout << "Warning: degree inconsistency" << endl;
			if(pVertex->id != itNode->first) cout << "Warning: vertex id inconsistency" << endl;
		#endif //CONTROL
		if((vScore < scoreNodeMin) or (vScore > scoreNodeMax) or (vDegree < degreeNodeMin) or (vDegree > degreeNodeMax) ) {
			//network.removeVertex(itNode->first);
			listNodeToRemove.push_back(itNode->first);
		} else {
			for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
				pEdge = itEdge->second;
				eScore = pEdge->data.score;
				if(pEdge->idSource != itNode->first or pEdge->idTarget != itEdge->first) cout << "Warning: edge vertex ids inconsistency" << endl;
				if(eScore < scoreEdgeMin or eScore > scoreEdgeMax) {
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
}

void Netscore::calculateClusteringCoefficient() {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge, itEdgeNeighbor, itEndEdgeNeighbor;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
	MapIntToEdge *pMapEdge, *pMapEdgeNeighbor;
    Edge *pEdge;
    int nEdgeConnecting = 0;
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	pVertex = itNode->second;
    	pMapEdge = pVertex->pMapIdToEdge;
		#ifdef CONTROL
    		if((unsigned int)pVertex->degree != pMapEdge->size()) cout << "Warning: degree inconsistency" << endl;
		#endif //CONTROL
    	for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            pMapEdgeNeighbor = (*pMapNode)[pEdge->idTarget]->pMapIdToEdge;
            for(itEdgeNeighbor = pMapEdgeNeighbor->begin(), itEndEdgeNeighbor = pMapEdgeNeighbor->end(); itEdgeNeighbor != itEndEdgeNeighbor; ++itEdgeNeighbor) {
            	if(pMapEdge->find(itEdgeNeighbor->first) != itEndEdge) {
            		nEdgeConnecting++;
            	}
            }
        }
    	//pVertex->cClustering = 2 * nEdgeConnecting / combination2(pVertex->degree);
		#ifdef CONTROL
    		if(pVertex->id != itNode->first) cout << "node id inconsistency" << endl;
		#endif //CONTROL
        mapIdToCCoefficient[pVertex->id] = 2 * nEdgeConnecting / combination2(pVertex->degree);
    }
}

void Netscore::clusterNetwork() {
    hash_map<int, float>::iterator it, itEnd;
    hash_map<int, void*> mapIdProcessed;
    //hash_map<int, void*>::iterator itProcessed, itEndProcessed;
    vector< pair<float, int> > listCCoefficient;
    int vId;
    float cCoefficient;
    
    calculateClusteringCoefficient();
    
    for(it = mapIdToCCoefficient.begin(), itEnd = mapIdToCCoefficient.end(); it != itEnd; ++it) {
    	listCCoefficient.push_back(pair<float, int>(it->second, it->first));
    }
    sort(listCCoefficient.begin(), listCCoefficient.end());
    int idClusterCurrent = 1;
    for(int i = listCCoefficient.size()-1; i >= 0; i--) {
    	vId = listCCoefficient[i].second;
    	cCoefficient = listCCoefficient[i].first;
    	if(mapIdProcessed.find(vId) == mapIdProcessed.end()) {
    		mapIdProcessed[vId] = NULL;
    		mapIdToClusterId[vId] = idClusterCurrent;
    		mapClusterIdToId[idClusterCurrent].push_back(vId);
    		if(flagDebug) 
    			cout << vId << " in cluster " << idClusterCurrent << endl;
    		expandCluster(vId, cCoefficient, idClusterCurrent, &mapIdProcessed);
    		idClusterCurrent++;
    	}
	}
}

void Netscore::expandCluster(int vId, float cCoefficient, int idCluster, hash_map<int, void*> *pMapIdProcessed) {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
	MapIntToEdge *pMapEdge;
    Edge *pEdge;
	//for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
	pVertex = (*pMapNode)[vId]; //itNode->second;
	pMapEdge = pVertex->pMapIdToEdge;
    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
        pEdge = itEdge->second;
        if(pMapIdProcessed->find(pEdge->idTarget) == pMapIdProcessed->end()) {
        	if(mapIdToCCoefficient[pEdge->idTarget] <= 0.75*cCoefficient) {
        		(*pMapIdProcessed)[pEdge->idTarget] = NULL;
        		mapIdToClusterId[pEdge->idTarget] = idCluster;
        		mapClusterIdToId[idCluster].push_back(pEdge->idTarget);
        		if(flagDebug)
        			cout << pEdge->idTarget << " in cluster " << idCluster << endl;
        		expandCluster(pEdge->idTarget, mapIdToCCoefficient[pEdge->idTarget], idCluster, pMapIdProcessed);
        	}
        }
    }
	//}
}

void Netscore::calculateClusterBasedNodeZScores() {
	hash_map<int, vector<int> >::iterator it, itEnd;
	vector<int> listId;
	int * pHistogramDegree;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
	Vertex *pVertex;
	float sumDegree = 0.0, sumDegreeSquare = 0.0;
	clusterNetwork();
	for(it = mapClusterIdToId.begin(), itEnd = mapClusterIdToId.end(); it != itEnd; ++it) {
		listId = it->second;
		pHistogramDegree = new int[listId.size()];
		for(unsigned int i=0; i < listId.size(); i++) pHistogramDegree[i] = 0; 
		for(unsigned int i=0; i < listId.size(); i++) {
			pVertex = (*pMapNode)[listId[i]];
			pHistogramDegree[pVertex->degree]++;
		}
		for(unsigned int i=0; i < listId.size(); i++) {
			pVertex = (*pMapNode)[listId[i]];
			sumDegree += (pHistogramDegree[pVertex->degree] / listId.size()) * pVertex->degree;
			sumDegreeSquare += (pHistogramDegree[pVertex->degree] / listId.size()) * sq(pVertex->degree);
		}
		for(unsigned int i=0; i < listId.size(); i++) {
			pVertex = (*pMapNode)[listId[i]];
			pVertex->data.score *= (pVertex->degree - sumDegree) / float(sqrt(sumDegreeSquare - sq(sumDegree)));
		}
	}
}

void Netscore::sampleNetwork(int nSample, hash_map<int, vector<int> > *pMapIdToListDegree) {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
	vector<int> listIdNode;
	srand((unsigned)time(0));
	int i, j, r, r2, r3;
	for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
		listIdNode.push_back(itNode->first);
	}
	#ifdef CONTROL
		if(listIdNode.size() != network.nNode) cout << "Warning: nNode inconsistency" << endl;
	#endif //CONTROL
    for(i=0; i < nSample; i++){
    	for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    		((*pMapIdToListDegree)[itNode->first]).push_back(itNode->second->degree);
    	}
        r = rand() % listIdNode.size(); //(rand()%10)+1;
        for(j=0; j<r; j++) {
        	r2 = rand() % listIdNode.size();
        	while(true) {
        		r3 = rand() % listIdNode.size();
        		if(r2 == r3) continue;
        		if(listIdNode[r2] > 0) {
        			(*pMapIdToListDegree)[listIdNode[r2]][i]--;
        			(*pMapIdToListDegree)[listIdNode[r3]][i]++;
	        		break;
        		} else {
        			r2 = rand() % listIdNode.size();
        		}
        	} 
        }
    } 
}

void Netscore::calculateNodeZScores(int nSample) {
	hash_map<int, vector<int> >::iterator itDegree, itEndDegree;
	hash_map<int, int>::iterator itHistogram, itEndHistogram;
	hash_map<int, vector<int> > mapIdToListDegree;
	hash_map<int, int> mapHistogramDegree;
	vector<int> listDegree;	
	//int maxDegree = -INT_INFINITY;
	unsigned int i = 0;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
	Vertex *pVertex;
	float sumDegree = 0.0, sumDegreeSquare = 0.0;
	
	sampleNetwork(nSample, &mapIdToListDegree);
	
	for(itDegree = mapIdToListDegree.begin(), itEndDegree = mapIdToListDegree.end(); itDegree != itEndDegree; ++itDegree) {
		listDegree = itDegree->second;
		for(i=0; i<listDegree.size(); i++) {
			//if(listDegree[i]>maxDegree)
			//	maxDegree = listDegree[i];
			if(mapHistogramDegree.find(listDegree[i]) == mapHistogramDegree.end()) 
				mapHistogramDegree[listDegree[i]] = 1;
			else
				mapHistogramDegree[listDegree[i]]++;
		}		
		//maxDegree = abs(maxDegree);
		sumDegree = 0.0;
		sumDegreeSquare = 0.0;
		for(itHistogram = mapHistogramDegree.begin(), itEndHistogram = mapHistogramDegree.end(); itHistogram != itEndHistogram; ++itHistogram) {	    
			sumDegree += (itHistogram->second / listDegree.size()) * itHistogram->first;
			sumDegreeSquare += (itHistogram->second / listDegree.size()) * sq(itHistogram->first);
		}
		pVertex = (*pMapNode)[itDegree->first];
		pVertex->data.score *= (pVertex->degree - sumDegree) / float(sqrt(sumDegreeSquare - sq(sumDegree)));
	}
}

void Netscore::scaleScores(bool flagUseZScoring) { //float scoreAllowedMaxNode, float scoreAllowedMaxEdge, bool flagUseZScoring) {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
	MapIntToEdge *pMapEdge;
    Edge *pEdge;
    float scoreMaxNode = -INFINITY, scoreMaxEdge = -INFINITY, scoreMinNode = INFINITY;
    
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	pVertex = itNode->second;
    	//if(pVertex->data.score > scoreAllowedMaxNode) {
    	//	pVertex->data.score = scoreAllowedMaxNode;
    	//}
    	//if(pVertex->data.score > scoreMaxNode) {
        //	scoreMaxNode = pVertex->data.score;
        //}
    	pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            //if(pEdge->data.score > scoreAllowedMaxEdge) {
        	//	pEdge->data.score = scoreAllowedMaxEdge;
        	//}
            if(pEdge->data.score > scoreMaxEdge) {
            	scoreMaxEdge = pEdge->data.score;
            }
        }
    }
    
    if(flagUseZScoring) {
    	calculateNodeZScores(100);
    }
    
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	pVertex = itNode->second;
    	if(pVertex->data.score > scoreMaxNode) {
        	scoreMaxNode = pVertex->data.score;
        } else if(pVertex->data.score < scoreMinNode) {
        	scoreMinNode = pVertex->data.score;
        }
    }
    
    if(flagUseZScoring) {
    	if(scoreMinNode > 0) cout << "Warning: positive min value after z-scoring" << endl;
    }
    
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	pVertex = itNode->second;
    	if(flagUseZScoring) {
    		pVertex->data.score = (pVertex->data.score + float(fabs(scoreMinNode)) + MIN_SCORE_AFTER_SHIFT) / (scoreMaxNode-scoreMinNode);
    	} else {
    		pVertex->data.score /= scoreMaxNode;
    	}
    	pVertex->data.scoreInitial = pVertex->data.score; 
    	pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            pEdge->data.score /= scoreMaxEdge;
            pEdge->data.scoreInitial = pEdge->data.score; // !! BO code renginering
	    }
    }
}

void Netscore::resetScores() {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
	MapIntToEdge *pMapEdge;
    Edge *pEdge;
	for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        pVertex->data.score = pVertex->data.scoreInitial; 
        pVertex->data.scoreAccumulatedLast = 0.0;
        pVertex->data.scoreDual = 0.0;
        pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            pEdge->data.score = pEdge->data.scoreInitial;
            pEdge->data.scoreAccumulatedLast = 0.0;
        }
    }
}

void Netscore::updateNetwork(int type, int typeTransfer, int tDegreeLinkerNodeMin, int tDegreeLinkerNodeMinIncrease, int tDegreeLinkerEdgeMin, int tDegreeLinkerEdgeMinIncrease, float scoreMaxAllowedNode, float scoreMaxAllowedEdge, float tScoreNeighborMin, float tScoreNeighborMax) {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
	MapIntToEdge *pMapEdge;
    Edge *pEdge;
    hash_map<int,void*> mapIdProcessed;
	
	for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
		updateNodeScore(type, typeTransfer, itNode->first, tDegreeLinkerNodeMin, tDegreeLinkerNodeMinIncrease, scoreMaxAllowedNode, tScoreNeighborMin, tScoreNeighborMax);
	}
	if(flagDebug)	
		printNetwork("After updateNodeScore");
	for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
		mapIdProcessed[itNode->first] = NULL;
        pVertex = itNode->second;
        pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
	    	pEdge = itEdge->second;
	    	if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
	    		updateEdgeScore(type, typeTransfer, pEdge->idSource, pEdge->idTarget, tDegreeLinkerEdgeMin, tDegreeLinkerEdgeMinIncrease, scoreMaxAllowedEdge, tScoreNeighborMin, tScoreNeighborMax);
	    	}
        }
    }
	if(flagDebug)
		printNetwork("After updateEdgeScore");
	if(type == STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION) {
		for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
			updateDualScore(itNode->first, tDegreeLinkerNodeMinIncrease, scoreMaxAllowedNode);
		}
		if(flagDebug)	
			printNetwork("After updateDualScore");
	}		
}

float Netscore::transferScore(int typeTransfer, float score, float a, float b) {
    switch(typeTransfer) {
    case IDENTITY:
        return score;
    case POLYNOMIAL:
        return a*score+b; // a=0.5, b=0.1, 0.2 - 0.6
    case LOGARITHMIC:
        return float((-1/10)*log(score)); // 0.0 - 1.0 "J"
    case EXPONENTIAL:
        return float(exp(0.75*score)/10); // 0.1 - 0.2 //float(1-exp(-0.75*score)); // 0.1 - 0.5
    //case RADIAL GAUSSIAN etc:
    //    return score;
    default:
        cout << "Warning: unidentified type" << endl;
        return 0.0;
    }
}

float Netscore::calculateScoreDifference() {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
	MapIntToEdge *pMapEdge;
    Edge *pEdge;
    hash_map<int,void*> mapIdProcessed;
    float sumDifferenceNode = 0.0, sumDifferenceEdge = 0.0;
    
	for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
		mapIdProcessed[itNode->first] = NULL;
        pVertex = itNode->second;
        sumDifferenceNode += (pVertex->data.score - pVertex->data.scoreInitial);
        pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
	    	pEdge = itEdge->second;
	    	if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
	    		sumDifferenceEdge += (pEdge->data.score - pEdge->data.scoreInitial);
	    	}
        }
    }
	return sumDifferenceEdge;
}

float Netscore::calculateError(int type) {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex, *pNeighbor;
	MapIntToEdge *pMapEdge;
    Edge *pEdge;
    hash_map<int,void*> mapIdProcessed;
    unsigned int nNode = 0, nEdge = 0;
    float error = 0.0, sumErrorNode = 0.0, sumErrorEdge = 0.0;
    
	for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
		mapIdProcessed[itNode->first] = NULL;
        pVertex = itNode->second;
        //error = pVertex->data.scoreInitial - pVertex->data.score; // !! BO code rengineering
        error = pVertex->data.scoreUpdated - pVertex->data.score; // decideErrorScoreNode(type, pVertex);
        sumErrorNode += sq(error);
        nNode++;        
        pVertex->data.scoreAccumulatedLast = pVertex->data.scoreUpdated - pVertex->data.score;
        pVertex->data.score = pVertex->data.scoreUpdated;
        pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
	    	pEdge = itEdge->second;
	    	if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
	            //error = pEdge->data.scoreUpdated - pEdge->data.scoreInitial; // !! BO code rengineering
                error = pEdge->data.scoreUpdated - pEdge->data.score; // decideErrorScoreEdge(type, pEdge);
	            sumErrorEdge += sq(error);
	            nEdge++;
	            pEdge->data.scoreAccumulatedLast = pEdge->data.scoreUpdated - pEdge->data.score; 
	            pEdge->data.score = pEdge->data.scoreUpdated;	            
	            pNeighbor = (*pMapNode)[pEdge->idTarget];
	           (*(pNeighbor->pMapIdToEdge))[pEdge->idSource]->data.score = pEdge->data.scoreUpdated;
	           (*(pNeighbor->pMapIdToEdge))[pEdge->idSource]->data.scoreAccumulatedLast = pEdge->data.scoreAccumulatedLast;
	            //cout << "e: " << error << " s: " << sumErrorEdge << "\t";
	    	}
        }
    }
	#ifdef CONTROL
		if(nNode != network.nNode) cout << "Warning: nNode inconsistency" << nNode << " " << network.nNode << endl;
		if(nEdge != network.nEdge) cout << "Warning: nEdge inconsistency" << nEdge << " " << network.nEdge << endl;
	#endif //CONTROL
    return max(float(sqrt( (sumErrorEdge/nEdge)*(sumErrorNode/nNode) )), float( (sqrt(sumErrorNode/nNode)+ sqrt(sumErrorEdge/nEdge)) /2)); //float(sqrt(sumErrorEdge/nEdge)); //float(sqrt(sumErrorNode/nNode)+sqrt(sumErrorEdge/nEdge))/2;
}

/*
float Netscore::decideErrorScoreNode(int type, Vertex *pVertex) {
    switch(type) {
    case STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: { // !! BO code rengineering
        return pVertex->data.scoreUpdated - pVertex->data.score; //pVertex->data.scoreAccumulatedLast; // pVertex->data.scoreUpdated - pEdge->data.scoreInitial;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: {
        return pVertex->data.scoreUpdated - pVertex->data.score;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_CHECK_DEGREE_DIVISION: {
        return pVertex->data.scoreUpdated - pVertex->data.score;
    }
    case STEPWISE_CONTRIBUTION_LINKER_PARAMETER_CHECK_DEGREE_DIVISION: {
        return pVertex->data.scoreUpdated - pVertex->data.score;
    }
    default:
        cout << "Warning: unidentified type" << endl;
        return 0.0;
    }
}

float Netscore::decideErrorScoreEdge(int type, Edge *pEdge) {
    switch(type) {
    case STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: { // !! BO code rengineering
    	return pEdge->data.scoreUpdated - pEdge->data.score; 
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: {
    	return pEdge->data.scoreUpdated - pEdge->data.score;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_CHECK_DEGREE_DIVISION: {
    	return pEdge->data.scoreUpdated - pEdge->data.score;
    }
    case STEPWISE_CONTRIBUTION_LINKER_PARAMETER_CHECK_DEGREE_DIVISION: {
        return pEdge->data.scoreUpdated - pEdge->data.score;
    }
    default:
        cout << "Warning: unidentified type" << endl;
        return 0.0;
    }
}
*/

void Netscore::updateNodeScore(int type, int typeTransfer, int vId, int tDegreeLinkerNodeMin, int tDegreeLinkerNodeMinIncrease, float scoreAllowedMaxNode, float tScoreNeighborMin, float tScoreNeighborMax) {
    MapIntToVertex *pMapNode = network.pMapIdToNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	Vertex *pVertex = (*pMapNode)[vId], *pNeighbor;
	MapIntToEdge *pMapEdge = pVertex->pMapIdToEdge;
	Edge *pEdge;
	int nCount = 0;
    float scoreSum = 0.0, scoreAccumulated = 0.0;
    //scoreSum = pVertex->data.scoreInitial;
    //scoreSum = decideInitialScoreNode(type, pVertex);
	for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
		pEdge = itEdge->second;
		pNeighbor = (*pMapNode)[itEdge->first];
		//scoreSum += pEdge->data.score * pNeighbor->data.score / tDegreeLinkerNodeMin;	
		scoreAccumulated += transferScore(typeTransfer, decideStepScoreNode(type, pEdge, pNeighbor, tDegreeLinkerNodeMin, tScoreNeighborMin, tScoreNeighborMax, &nCount));
	}
	//scoreSum += scoreAccumulated;
    scoreSum = decideScoreNode(type, pVertex, scoreAccumulated, nCount, tDegreeLinkerNodeMin); // decideAccumulationScoreNode(type, scoreSum, nCount);
	pVertex->data.scoreUpdated = (scoreSum <= scoreAllowedMaxNode) ? scoreSum : scoreAllowedMaxNode;
	if(STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION) {
		pVertex->data.scoreDual = pVertex->data.scoreUpdated - pVertex->data.score; 
	}
}

/*
float Netscore::decideInitialScoreNode(int type, Vertex *pVertex) {
    switch(type) {
    case STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: { // !! BO code rengineering
        return pVertex->data.score - pVertex->data.scoreAccumulatedLast;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: {
    	return pVertex->data.score - pVertex->data.scoreAccumulatedLast;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_CHECK_DEGREE_DIVISION: {
	    return pVertex->data.score;
    }
    case STEPWISE_CONTRIBUTION_LINKER_PARAMETER_CHECK_DEGREE_DIVISION: {
	    return 0.0;
    }
    default:
        cout << "Warning: unidentified type" << endl;
        return 0.0;
    }
}
*/

float Netscore::decideStepScoreNode(int type, Edge *pEdge, Vertex *pNeighbor, int tDegreeLinkerNodeMin, float tScoreNeighborMin, float tScoreNeighborMax, int *nCount) {
    switch(type) {
    case STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: { // !! BO code rengineering
        return pEdge->data.score * pNeighbor->data.score / tDegreeLinkerNodeMin;	
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: {
        return pEdge->data.score * pNeighbor->data.score / tDegreeLinkerNodeMin;	
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_CHECK: {
        if(pNeighbor->data.score >= tScoreNeighborMin and pNeighbor->data.score <= tScoreNeighborMax) {
            (*nCount)++;
            return pEdge->data.score * pNeighbor->data.score;
        }
    }
    case STEPWISE_CONTRIBUTION_LINKER_PARAMETER_CHECK: {
        if(pNeighbor->data.score >= tScoreNeighborMin and pNeighbor->data.score <= tScoreNeighborMax) {
            (*nCount)++;
            return pEdge->data.score * pNeighbor->data.score;
        }
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_SCORE_CHECK_DEGREE_DIVISION: {
        if(pNeighbor->data.score >= tScoreNeighborMin and pNeighbor->data.score <= tScoreNeighborMax) {
            (*nCount)++;
            return pEdge->data.score * pNeighbor->data.score;
        }
    }
    case STEPWISE_CONTRIBUTION_LINKER_SCORE_CHECK_DEGREE_DIVISION: {
        if(pNeighbor->data.score >= tScoreNeighborMin and pNeighbor->data.score <= tScoreNeighborMax) {
            (*nCount)++;
            return pEdge->data.score * pNeighbor->data.score;
        }
    }    
    case ACCUMULATING_CONTRIBUTION_DEGREE_DIVISION: {
        (*nCount)++;
        return pEdge->data.score * pNeighbor->data.score;
    }
    default:
        cout << "Warning: unidentified type" << endl;
        return 0.0;
    }
}

float Netscore::decideScoreNode(int type, Vertex *pVertex, float scoreAccumulated, int nCount, int tDegreeLinkerNodeMin) { // decideAccumulationScoreNode(int type, float sumAccumulation, int nCount) {
    switch(type) {
    case STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: { // !! BO code rengineering
    	//float scoreDualLastAccumulated = pVertex->data.scoreDual;
    	//pVertex->data.scoreDual = scoreAccumulated;
	    return (pVertex->data.score - pVertex->data.scoreDual) + scoreAccumulated; //- scoreDualLastAccumulated;
    	//return (pVertex->data.score - scoreDualLastAccumulated) + scoreAccumulated;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: {
    	return (pVertex->data.score - pVertex->data.scoreAccumulatedLast) + scoreAccumulated;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_CHECK: {
    	if(nCount >= tDegreeLinkerNodeMin) {
    		return pVertex->data.score + scoreAccumulated;
    	} else {
    		return pVertex->data.score;
    	}
    }
    case STEPWISE_CONTRIBUTION_LINKER_PARAMETER_CHECK: {
    	if(pVertex->data.scoreUpdated == 0.0) { // if first iteration
    		if(nCount >= tDegreeLinkerNodeMin) {
	    		return pVertex->data.score + scoreAccumulated;
	    	} else {
	    		return pVertex->data.score;
	    	}
    	} else {
			if(nCount >= tDegreeLinkerNodeMin) {
				return scoreAccumulated;
			} else {
				return 0;
			}
    	}
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_SCORE_CHECK_DEGREE_DIVISION: {
	    return pVertex->data.score + (scoreAccumulated / nCount);
    }
    case STEPWISE_CONTRIBUTION_LINKER_SCORE_CHECK_DEGREE_DIVISION: {
    	if(pVertex->data.scoreUpdated == 0.0) { // if first iteration
    		return pVertex->data.score + scoreAccumulated / nCount;
    	}
	    return scoreAccumulated / nCount;
    }
    case ACCUMULATING_CONTRIBUTION_DEGREE_DIVISION: {
    	return pVertex->data.score + (scoreAccumulated / nCount);
    }
    default:
        cout << "Warning: unidentified type" << endl;
        return 0.0;
    }
}

void Netscore::updateEdgeScore(int type, int typeTransfer, int idSource, int idTarget, int tDegreeLinkerEdgeMin, int tDegreeLinkerEdgeMinIncrease, float scoreAllowedMaxEdge, float tScoreNeighborMin, float tScoreNeighborMax) {
    MapIntToVertex *pMapNode = network.pMapIdToNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	Vertex *pSource = (*pMapNode)[idSource], *pTarget = (*pMapNode)[idTarget];
	MapIntToEdge *pMapEdgeSource = pSource->pMapIdToEdge, *pMapEdgeTarget = pTarget->pMapIdToEdge;
	Edge *pEdgeSource = (*pMapEdgeSource)[idTarget], *pEdgeTarget = (*pMapEdgeTarget)[idSource], *pEdge;
    int eCount = 0;
    float scoreSum = 0.0, scoreAccumulated = 0.0;
	if(pEdgeSource->data.score != pEdgeTarget->data.score) {
		cout << "Warning: edge score inconsistency" << endl; 
	}
	// float scoreSum = pEdgeSource->data.scoreInitial; // !! BO code rengineering
    //scoreSum = decideInitialScoreEdge(type, pEdgeSource);
	for(itEdge = pMapEdgeSource->begin(), itEndEdge = pMapEdgeSource->end(); itEdge != itEndEdge; ++itEdge) {
        if(itEdge->first == idTarget) {
            continue;
        }
		pEdge = itEdge->second;
		// scoreSum += transferScore(pSource->data.score * pEdge->data.score / tDegreeLinkerEdgeMin);
		scoreAccumulated += transferScore(typeTransfer, decideStepScoreEdge(type, pSource, pEdge, tDegreeLinkerEdgeMin, tScoreNeighborMin, tScoreNeighborMax, &eCount)); // scoreSum += transferScore(IDENTITY, decideStepScoreEdge(type, pSource, pEdge, tDegreeLinkerEdgeMin, tScoreNeighborMin, tScoreNeighborMax, &eCount));
	}
	for(itEdge = pMapEdgeTarget->begin(), itEndEdge = pMapEdgeTarget->end(); itEdge != itEndEdge; ++itEdge) {
        if(itEdge->first == idSource) {
            continue;
        }
		pEdge = itEdge->second;
		//scoreSum += transferScore(pTarget->data.score * pEdge->data.score / tDegreeLinkerEdgeMin); 
		scoreAccumulated += transferScore(typeTransfer, decideStepScoreEdge(type, pTarget, pEdge, tDegreeLinkerEdgeMin, tScoreNeighborMin, tScoreNeighborMax, &eCount)); // scoreSum += transferScore(IDENTITY, decideStepScoreEdge(type, pTarget, pEdge, tDegreeLinkerEdgeMin, tScoreNeighborMin, tScoreNeighborMax, &eCount));
	}
    scoreSum = decideScoreEdge(type, pEdgeSource, scoreAccumulated, eCount, tDegreeLinkerEdgeMin); //decideAccumulationScoreEdge(type, scoreSum, eCount);
	if(scoreSum <= scoreAllowedMaxEdge) {
		pEdgeSource->data.scoreUpdated = scoreSum;
		pEdgeTarget->data.scoreUpdated = scoreSum; 
	} else {
		pEdgeSource->data.scoreUpdated = scoreAllowedMaxEdge;
		pEdgeTarget->data.scoreUpdated = scoreAllowedMaxEdge;
	}
	if(STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION) {
		float scoreTemp = pEdgeSource->data.scoreUpdated + (pEdgeSource->data.scoreUpdated / tDegreeLinkerEdgeMinIncrease);
		if(scoreTemp <= scoreAllowedMaxEdge) {
			pEdgeSource->data.scoreUpdated = scoreTemp;
			pEdgeTarget->data.scoreUpdated = scoreTemp; 
		} else {
			pEdgeSource->data.scoreUpdated = scoreAllowedMaxEdge;
			pEdgeTarget->data.scoreUpdated = scoreAllowedMaxEdge;
		}
	}
}

/*
float Netscore::decideInitialScoreEdge(int type, Edge *pEdge) {
    switch(type) {
    case STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: { // !! BO code rengineering
        return pEdge->data.score - pEdge->data.scoreAccumulatedLast; // pEdge->data.scoreInitial;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: {
	    return pEdge->data.score - pEdge->data.scoreAccumulatedLast; // pEdge->data.score;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_CHECK_DEGREE_DIVISION: {
	    return pEdge->data.score;
    }
    case STEPWISE_CONTRIBUTION_LINKER_PARAMETER_CHECK_DEGREE_DIVISION: {
	    return 0.0;
    }
    default:
        cout << "Warning: unidentified type" << endl;
        return 0.0;
    }
}
*/

float Netscore::decideStepScoreEdge(int type, Vertex *pNeighbor, Edge *pEdge, int tDegreeLinkerEdgeMin, float tScoreNeighborMin, float tScoreNeighborMax, int *eCount) {
    switch(type) {
    case STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: { // !! BO code rengineering
        return pNeighbor->data.score * pEdge->data.score / tDegreeLinkerEdgeMin;	
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: {
        return pNeighbor->data.score * pEdge->data.score / tDegreeLinkerEdgeMin;	
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_CHECK: {
        if(pEdge->data.score >= tScoreNeighborMin and pEdge->data.score <= tScoreNeighborMax) {
            (*eCount)++;
            return pNeighbor->data.score * pEdge->data.score;
        }
    }
    case STEPWISE_CONTRIBUTION_LINKER_PARAMETER_CHECK: {
        if(pEdge->data.score >= tScoreNeighborMin and pEdge->data.score <= tScoreNeighborMax) {
            (*eCount)++;
            return pNeighbor->data.score * pEdge->data.score;
        }
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_SCORE_CHECK_DEGREE_DIVISION: {
        if(pEdge->data.score >= tScoreNeighborMin and pEdge->data.score <= tScoreNeighborMax) {
            (*eCount)++;
            return pNeighbor->data.score * pEdge->data.score;
        }
    }
    case STEPWISE_CONTRIBUTION_LINKER_SCORE_CHECK_DEGREE_DIVISION: {
        if(pEdge->data.score >= tScoreNeighborMin and pEdge->data.score <= tScoreNeighborMax) {
            (*eCount)++;
            return pNeighbor->data.score * pEdge->data.score;
        }
    }
    case ACCUMULATING_CONTRIBUTION_DEGREE_DIVISION: {
        (*eCount)++;
        return pNeighbor->data.score * pEdge->data.score;
    }
    default:
        cout << "Warning: unidentified type" << endl;
        return 0.0;
    }
}

float Netscore::decideScoreEdge(int type, Edge *pEdge, float scoreAccumulated, int eCount, int tDegreeLinkerEdgeMin) {
    switch(type) {
    case STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: { // !! BO code rengineering
	    return (pEdge->data.score - pEdge->data.scoreAccumulatedLast) + scoreAccumulated;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION: {
	    return (pEdge->data.score - pEdge->data.scoreAccumulatedLast) + scoreAccumulated;
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_CHECK: {
    	if(eCount >= tDegreeLinkerEdgeMin) {
    		return pEdge->data.score + scoreAccumulated;
    	} else {
    		return pEdge->data.score;
    	}
    }
    case STEPWISE_CONTRIBUTION_LINKER_PARAMETER_CHECK: {
    	if(pEdge->data.scoreUpdated) { // if first iteration
    		if(eCount >= tDegreeLinkerEdgeMin) {
	    		return pEdge->data.score + scoreAccumulated;
	    	} else {
	    		return pEdge->data.score;
	    	}
    	} else {
	    	if(eCount >= tDegreeLinkerEdgeMin) {
	    		return scoreAccumulated;
	    	} else {
	    		return 0;
	    	}
    	}
    }
    case ACCUMULATING_CONTRIBUTION_LINKER_SCORE_CHECK_DEGREE_DIVISION: {
	    return pEdge->data.score + (scoreAccumulated / eCount);
    }
    case STEPWISE_CONTRIBUTION_LINKER_SCORE_CHECK_DEGREE_DIVISION: {
    	if(pEdge->data.scoreUpdated) { // if first iteration
    		return pEdge->data.score + (scoreAccumulated / eCount);
    	}
	    return scoreAccumulated / eCount;
    }
    case ACCUMULATING_CONTRIBUTION_DEGREE_DIVISION: {
    	return pEdge->data.score + (scoreAccumulated / eCount);
    }
    default:
        cout << "Warning: unidentified type" << endl;
        return 0.0;
    }
}

void Netscore::updateDualScore(int vId, int tDegreeLinkerNodeMinIncrease, float scoreAllowedMaxNode) {
    MapIntToVertex *pMapNode = network.pMapIdToNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	Vertex *pVertex = (*pMapNode)[vId];
	float scoreSum = 0.0, scoreAccumulation = 0.0;
	pVertex->data.score = pVertex->data.scoreUpdated; // error will be wrong
	scoreAccumulation = pVertex->data.scoreUpdated * combination2(pVertex->degree) / 2 * tDegreeLinkerNodeMinIncrease;
    scoreSum = (pVertex->data.scoreUpdated - pVertex->data.scoreAccumulatedLast) + scoreAccumulation;
	pVertex->data.scoreUpdated = (scoreSum <= scoreAllowedMaxNode) ? scoreSum : scoreAllowedMaxNode;
}

int Netscore::combination2(int x) {
	int result=1;
	if(x>2) {
		for(int i=x; i >= x-2; i--) {
			result *= i;
		}
		result /= 2;
	}
	return result;
}

void Netscore::printNetwork(string explanation) {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
	MapIntToEdge *pMapEdge;
    Edge *pEdge;
    hash_map<int,void*> mapIdProcessed;
    
    if(explanation != "") cout << "----" << explanation << ": " << endl;
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	mapIdProcessed[itNode->first] = NULL;
    	pVertex = itNode->second;
    	pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
	            /*
	            cout << network.mapIdToName[pEdge->idSource] << " " << (*pMapNode)[pEdge->idSource]->data.scoreUpdated 
	            	 << "(" << (*pMapNode)[pEdge->idSource]->data.score << " / " << (*pMapNode)[pEdge->idSource]->data.scoreInitial << ")" << "\t" 
	            	 << network.mapIdToName[pEdge->idTarget] << " " << (*pMapNode)[pEdge->idTarget]->data.scoreUpdated 
	            	 << "(" << (*pMapNode)[pEdge->idTarget]->data.score << " / " << (*pMapNode)[pEdge->idTarget]->data.scoreInitial << ")" << "\t" 
	                 << "\t" << pEdge->data.scoreUpdated << "(" << pEdge->data.score << " / " << pEdge->data.scoreInitial << ")" << endl;
	            */ 
	            cout << network.mapIdToName[pEdge->idSource] << " " << (*pMapNode)[pEdge->idSource]->data.scoreUpdated 
	 	        	 << "(" << (*pMapNode)[pEdge->idSource]->data.score << ")\t" 
	 	        	 << network.mapIdToName[pEdge->idTarget] << " " << (*pMapNode)[pEdge->idTarget]->data.scoreUpdated 
		        	 << "(" << (*pMapNode)[pEdge->idTarget]->data.score << ")\t" 
                     << "\t" << pEdge->data.scoreUpdated << "(" << pEdge->data.score << ")" << endl;	            
            }
        }
    }
}

void Netscore::printNodes(string explanationPrefix) {
	cout << explanationPrefix << network << endl;
}

void Netscore::printEdges() {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
	MapIntToEdge *pMapEdge;
    Edge *pEdge;
    hash_map<int,void*> mapIdProcessed;
    
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	mapIdProcessed[itNode->first] = NULL;
    	pVertex = itNode->second;
    	pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
	            cout << network.mapIdToName[pEdge->idSource] << "\t" << network.mapIdToName[pEdge->idTarget] 
	                 << "\t" << pEdge->data.scoreUpdated << "(" << pEdge->data.score << " / " << pEdge->data.scoreInitial << ")" << endl; 
            }
        }
    }
}

float Netscore::getEdgeScore(string vNameSource, string vNameTarget) {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = network.pMapIdToNode;
	MapIntToEdge *pMapEdge;
	MapIntToString mapName = network.mapIdToName;
	Vertex *pVertex;
	Edge *pEdge; 
	for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	pVertex = itNode->second;
    	if(mapName[pVertex->id] != vNameSource)
    		continue;
    	pMapEdge = pVertex->pMapIdToEdge;
    	for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            if(mapName[pEdge->idTarget] == vNameTarget) {
            	//if(pEdge->data.score != scoreMaxAllowedEdge) {
            		//cout << type << " " << typeTransfer << " " << nIteration << " " << tError << " " << tDegreeLinkerNodeMin << " " << tDegreeLinkerNodeMinIncrease << " " << tDegreeLinkerEdgeMin << " " << tDegreeLinkerEdgeMinIncrease << " " << scoreMaxAllowedNode << " " << scoreMaxAllowedEdge << " " << tScoreNeighborMin << " " <<  tScoreNeighborMax << endl;
            		//cout << pEdge->data.score << endl;
            		//flagPrint = true;
	            //}
    		}
    	}
	}
	return pEdge->data.score;
}
