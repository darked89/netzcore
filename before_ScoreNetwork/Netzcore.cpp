#include "Netzcore.h"
#include <fstream>
#include <cstdlib> 
#include <ctime>
#include <sstream>

Netzcore::Netzcore() {
    flagDebug = false;
    flagVerbose = false;
    minZScore = INFINITY;
    maxZScore = -INFINITY;
    minScore = INFINITY;
    maxScore = -INFINITY;
}

Netzcore::Netzcore(string fileProtein, string fileInteraction, bool fDebug, bool fVerbose) {
    loadProteins(fileProtein); 
    loadInteractions(fileInteraction); 
    flagDebug = fDebug;
    flagVerbose = fVerbose;
    minZScore = INFINITY;
    maxZScore = -INFINITY;
    minScore = INFINITY;
    maxScore = -INFINITY;
}

Netzcore::~Netzcore() {
}

void Netzcore::run(int nIteration, float tError, int typeRandomization, int typeTransfer, int nSample, bool flagSampleOnce, bool flagUseEdgeReliabilityScore, bool flagNormalizeOnce, bool flagAccumulate, bool flagResetSeedScoresToInitial, bool flagUpdateSampledOnce, bool flagScaleZScoreWithNodeDegree, bool flagNetscore, bool flagSaveRandomNetworks, bool flagLoadRandomNetworks) { 
    float error = INFINITY;
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
    for(i=0; i<nIteration and error > tError; i++) {
        if(flagDebug)
            cout << endl << "Iteration: " << i << endl;  
        if(flagVerbose)
            cerr << endl << "Iteration: " << i << endl;  
        /// handle first iteration seperately - sample and/or update
        if(i == 0) {
            /// either flagSampleOnce true or false
            if(!flagNetscore) {
                if(flagLoadRandomNetworks) {
                    cerr << "loading.." << endl;
                    readSampledNetworks(nSample);
                    updateSampledNetworkNodeScores(); 
	    	    //printNetwork("AfterUpdateSampled");
                    //printSampledNetworks();
                } else {
                    cerr << "sampling.." << endl;
                    sampleNetwork(nSample, typeRandomization);
                    updateSampledNetworkNodeScores(); 
	    	    //printNetwork("AfterSampleUpdate");
                    //printSampledNetworks();
                    if(flagSaveRandomNetworks) {
                        cerr << "saving.." << endl;
                        writeSampledNetworks();
                    }
                }
            }
            if(flagDebug)
            	printSampledNetworks();
            if(flagNetscore) {
    	        updateNetwork(typeTransfer, i, flagUseEdgeReliabilityScore, true, flagAccumulate, flagScaleZScoreWithNodeDegree); // flagNoZScoring:true dont use Z scoring
            } else {
            	updateNetwork(typeTransfer, i, flagUseEdgeReliabilityScore, false, flagAccumulate, flagScaleZScoreWithNodeDegree); // flagNoZScoring:false use Z scoring
            }
            /// calculate error & update scores
            error = calculateErrorAndUpdateScores(flagResetSeedScoresToInitial);
            /// update scores of sampled networks if sampling only once (flagSampleOnce is true) and independent of flagUpdateSampledOnce
            if(flagSampleOnce and !flagNormalizeOnce and !flagNetscore) {
                updateSampledNetworkNodeScores();
            }
        } else {
            if(!flagSampleOnce and !flagNetscore) {
                sampleNetwork(nSample, typeRandomization);
            }
            if(flagNetscore) {
    	        updateNetwork(typeTransfer, i, flagUseEdgeReliabilityScore, true, flagAccumulate, flagScaleZScoreWithNodeDegree); // flagNoZScoring:true dont use Z scoring
            } else {
                updateNetwork(typeTransfer, i, flagUseEdgeReliabilityScore, flagNormalizeOnce, flagAccumulate, flagScaleZScoreWithNodeDegree); /// flagNoZScoring = flagNormalizeOnce -> dont use Z scoring if normalize once is true after 1st iteration
            }
            /// calculate error & update scores
            error = calculateErrorAndUpdateScores(flagResetSeedScoresToInitial);
            /// update scores of sampled networks if sampling only once (flagSampleOnce is true) and flagUpdateSampledOnce is false
            if(flagSampleOnce and !flagUpdateSampledOnce and !flagNetscore) {
                updateSampledNetworkNodeScores();
            }
        } 
        if(flagDebug) 
            printNetwork("After updateNetwok");
        if(flagDebug)
            cout << "error: " << error << endl;
        cerr << "i: " << i << " error: " << error << endl;
        // delete sampled networks if sampling each step (flagSampleOnce is false)
        if((!flagSampleOnce and !flagNetscore) or flagNormalizeOnce) {
            deleteSampledNetworks();
        } 
    }	
    cerr << "nIteration: " << i << "\nerror: " << error << endl;
    scaleNodeScores();
    //if(flagDebug)
    //	printSampledNetworks();
    if(flagDebug)
        printNetwork();
    if(flagSampleOnce) {
        deleteSampledNetworks();
    }
}

void Netzcore::loadProteins(string const &fileName) {
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

void Netzcore::loadInteractions(string const &fileName) {
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

// filter network w.r.t. given parameters 
void Netzcore::filterNetwork(int degreeNodeMin, int degreeNodeMax, float scoreNodeMin, float scoreNodeMax, float scoreEdgeMin, float scoreEdgeMax, bool flagRemoveSelfEdges) {	
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


void Netzcore::resetScores() {
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

void Netzcore::sampleNetwork(int nSample, int typeRandomization) {
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToEdge::iterator itEdge, itEndEdge; //, itEdgeNeighbor, itEdgeNeighbor2;
    MapIntToVertex *pMapNode;
    MapIntToEdge *pMapEdge;//, *pMapEdgeNew, *pMapEdgeTemp;
    Vertex *pVertex, *pVertexNeighbor; //, *pVertexNew, ;
    //Edge *pEdge; //, *pEdgeNeighbor, *pEdgeNeighbor2;
    Graph *pNetworkSampled;
    //Data dataTemp;
    int i, j, k, nPerturbation, iNode, iNeighbor, iNodeNew, iNeighborNew, idNode, idNodeNew, idNeighbor, idNeighborNew;
    vector<int> listIdNode, listIdNodeTemp, listIdNodeTempNeighbor;
    hash_map<int, vector<int> *> mapDegreeToListIdNode;
    hash_map<int, hash_map<int, void*> *> mapDegreeToMapIdNode;
    vector<Data> listDataEdge;
    //hash_map<int, int> mapDegreeToCount;
    ostringstream oss;
    srand((unsigned)time(0));
    /// initialize helper storages
    listIdNode = network.getNodeIdList();
    if(typeRandomization == TOPOLOGY_PRESERVING_SAME_DEGREE_INTERCHANGE) {
        mapDegreeToListIdNode = network.getMapDegreeToNodeIdList();
        //hash_map<int, vector<int> *>::iterator itDegree;
        //for(itDegree = mapDegreeToListIdNode.begin(); itDegree != mapDegreeToListIdNode.end(); ++itDegree) {
        //    cout << itDegree->first << ": ";
        //    listIdNodeTemp = *(itDegree->second);
        //    for(i=0; i<(int)listIdNodeTemp.size(); i++) 
        //        cout << listIdNodeTemp[i] << " ";
        //    cout << endl;
        //}
    } else if(typeRandomization == RANDOM) {
        listDataEdge = network.getEdgeDataList();
        //for(i=0; i<(int)listDataEdge.size(); i++) 
        //    cout << listDataEdge[i].score << " ";
        //cout << endl;
        #ifdef CONTROL
            if((unsigned int)network.nEdge != listDataEdge.size()) cout << "Warning: number of edge inconsistency" << endl;
        #endif //CONTROL
    } else if(typeRandomization == DEGREE_DISTRIBUTION_PRESERVING) {
        //mapDegreeToCount = network.getDegreeHistogram(); // not needed since network is not created from scratch
        //mapDegreeToMapIdNode = network.getMapDegreeToNodeIdMap(); // problems for some degree switch to vector
        /// should be done for each sample
        //mapDegreeToListIdNode = network.getMapDegreeToNodeIdList();
	} else if(typeRandomization == DEGREE_DISTRIBUTION_PRESERVING_EDGE_INTERCHANGE) {
        mapDegreeToListIdNode = network.getMapDegreeToNodeIdList();
    }

    for(i=0; i < nSample; i++){
        //cout << i << " .... " << endl;
        if(typeRandomization == RANDOM) {
    	    pNetworkSampled = new Graph(network, true);
            for(k=0; (unsigned int)k < network.nEdge; k++) {
                // guarantee that graph will not be disconnected - if nNode close to nEdge may yield sparse graphs
                if(k < (int)listIdNode.size()) {
                    iNode = k;
                } else {
        	        iNode = rand() % listIdNode.size();
                }
        	    idNode = listIdNode[iNode];
                do {
                    do {
                        iNodeNew = rand() % listIdNode.size();
                    } while(iNodeNew == iNode);
                    idNodeNew = listIdNode[iNodeNew];
                    //cout << "add " << idNode << " " << idNodeNew << endl;
                }
                while(pNetworkSampled->addEdge(idNode, idNodeNew, listDataEdge[k]) == false);
            }
        } else {
            pNetworkSampled = new Graph(network);
        }
        if(typeRandomization == TOPOLOGY_PRESERVING or typeRandomization == TOPOLOGY_PRESERVING_SAME_DEGREE_INTERCHANGE) {
        	pNetworkSampled->initializeMapNewId();
        } else if(typeRandomization == DEGREE_DISTRIBUTION_PRESERVING) {
            mapDegreeToListIdNode = network.getMapDegreeToNodeIdList();
        }
    	pMapNode = pNetworkSampled->pMapIdToNode;
        nPerturbation = (rand() % min(network.nNode, network.nEdge)); // previously max
        if(nPerturbation < N_PERTURBATION_MIN) {
        	nPerturbation = N_PERTURBATION_MIN;
        }
        for(j=0; j<nPerturbation; j++) {
            oss.clear();
            oss << j;
            //printGraph(networkSampled, oss.str());
            iNode = rand() % listIdNode.size();
            idNode = listIdNode[iNode];
            pVertex = (*pMapNode)[idNode];
            if(typeRandomization == TOPOLOGY_PRESERVING_SAME_DEGREE_INTERCHANGE) {
                listIdNodeTemp = *(mapDegreeToListIdNode[pVertex->degree]);
                while(listIdNodeTemp.size() == 1) {
                    iNode = rand() % listIdNode.size();
                    idNode = listIdNode[iNode];
                    pVertex = (*pMapNode)[idNode];
                    listIdNodeTemp = *(mapDegreeToListIdNode[pVertex->degree]);
                }
                do {
                    iNodeNew = rand() % listIdNodeTemp.size();
                    //cout << iNodeNew << " " << listIdNoteTemp.size() << endl;
                } while(iNodeNew == iNode);
                idNodeNew = listIdNodeTemp[iNodeNew];
                pNetworkSampled->swapInMapNewId(idNode, idNodeNew);
            }
            else if(typeRandomization == RANDOM) {
                // created from scratch above
                break;
            }
            else if(typeRandomization == DEGREE_DISTRIBUTION_PRESERVING) {
                vector<int> * pListIdNodeTemp;
                /// select a node that has more than one edge 
                while(pVertex->degree == 1) {
                    iNode = rand() % listIdNode.size();
                    idNode = listIdNode[iNode];
                    pVertex = (*pMapNode)[idNode];
                }
            	pMapEdge = pVertex->pMapIdToEdge;
				#ifdef CONTROL
            		if((unsigned int)pVertex->degree != pMapEdge->size()) cout << "Warning: degree inconsistency" << endl;
				#endif //CONTROL
        		iNeighbor = rand() % pVertex->degree;
                Edge *pEdge = NULL;
                Data dataEdge;
        		for(k=0, itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge, k++) {
        			if(k == iNeighbor) {
                        pEdge = itEdge->second;
                        dataEdge = pEdge->data;
        				break;
                    }    
        		}
				#ifdef CONTROL
        			if(itEdge == itEndEdge) cerr << "Warning: iteration inconsistency" << endl;
                    if(pVertex->id == itEdge->first) cerr << "Warning: unexpected self edge" << endl;
				#endif //CONTROL
                idNeighbor = itEdge->first;
                pVertexNeighbor = (*pMapNode)[idNeighbor];
                //cout << "rm " << idNode << " " << idNeighbor << endl;
    			pNetworkSampled->removeEdge(pVertex->id, pVertexNeighbor->id); //(idNode, idNeighbor);
                //cout << pVertex->degree << " " << pVertexNeighbor->degree << endl;

                /// update degree to node list map 
                pListIdNodeTemp = mapDegreeToListIdNode[pVertex->degree+1];
                for(k=0; k<(int)pListIdNodeTemp->size(); k++) {
                    if(pListIdNodeTemp->at(k) == pVertex->id) {
                        pListIdNodeTemp->erase(pListIdNodeTemp->begin()+k);
                        break;
                    }
                }
                if(mapDegreeToListIdNode.find(pVertex->degree) == mapDegreeToListIdNode.end()) {
                    mapDegreeToListIdNode[pVertex->degree] = new vector<int>;
                }
                pListIdNodeTemp = mapDegreeToListIdNode[pVertex->degree];
                pListIdNodeTemp->push_back(pVertex->id);
                /// neighbor might have a degree one less for which there is no entry in map
                pListIdNodeTemp = mapDegreeToListIdNode[pVertexNeighbor->degree+1];
                for(k=0; k<(int)pListIdNodeTemp->size(); k++) {
                    if(pListIdNodeTemp->at(k) == pVertexNeighbor->id) {
                        pListIdNodeTemp->erase(pListIdNodeTemp->begin()+k);
                        break;
                    }
                }
                if(mapDegreeToListIdNode.find(pVertexNeighbor->degree) == mapDegreeToListIdNode.end()) {
                    mapDegreeToListIdNode[pVertexNeighbor->degree] = new vector<int>;
                }
                pListIdNodeTemp = mapDegreeToListIdNode[pVertexNeighbor->degree];
                pListIdNodeTemp->push_back(pVertexNeighbor->id);
                /// select new source and target for the edge
                int nLoop = 0;
                do {
                    if(nLoop > N_TRIAL_MAX) {
                        //cout << "max # of trial exceeded" << endl;
                        idNodeNew = pVertex->id; // for the assignment below
                        idNeighborNew = pVertexNeighbor->id;
                        pNetworkSampled->addEdge(idNodeNew, idNeighborNew, dataEdge);
                        break;
                    }
                    nLoop++;
                    pListIdNodeTemp = mapDegreeToListIdNode[pVertex->degree];
                    iNodeNew = rand() % pListIdNodeTemp->size();
                    for(k=0; k<(int)pListIdNodeTemp->size(); k++) {
                        if(k==iNodeNew) {
                            /// try to select a different node
                            if(pListIdNodeTemp->at(k) != pVertex->id) {
                                idNodeNew = pListIdNodeTemp->at(k);
                            } else {
                                if(k!=0) {
                                    idNodeNew = pListIdNodeTemp->at(k-1);
                                } else if(k != (int)pListIdNodeTemp->size()-1){
                                    idNodeNew = pListIdNodeTemp->at(k+1);
                                } else {
                                    idNodeNew = pVertex->id;
                                }
                            }
                            break;
                        }
                    }
                    pListIdNodeTemp = mapDegreeToListIdNode[pVertexNeighbor->degree];
                    iNeighborNew = rand() % pListIdNodeTemp->size();
                    for(k=0; k < (int)pListIdNodeTemp->size(); k++) {
                        if(k==iNeighborNew) {
                            /// try to select a different node
                            if(pListIdNodeTemp->at(k) != pVertexNeighbor->id) {
                                idNeighborNew = pListIdNodeTemp->at(k);
                            } else {
                                if(k!=0) { 
                                    idNeighborNew = pListIdNodeTemp->at(k-1);
                                } else if(k != (int)pListIdNodeTemp->size()-1){
                                    idNeighborNew = pListIdNodeTemp->at(k+1);
                                } else {
                                    idNeighborNew = pVertexNeighbor->id;
                                }
                            }
                            break;
                        }
                    }
                    //cout << "add " << idNodeNew << " " << idNeighborNew << endl;
                } while(idNodeNew == idNeighborNew or pNetworkSampled->addEdge(idNodeNew, idNeighborNew, dataEdge) == false); 
                pVertex = (*pMapNode)[idNodeNew];
                pVertexNeighbor = (*pMapNode)[idNeighborNew];
                //cout << pVertex->degree << " " << pVertexNeighbor->degree << endl;
                /// update degree to node list map 
                pListIdNodeTemp = mapDegreeToListIdNode[pVertex->degree-1];
                for(k=0; k<(int)pListIdNodeTemp->size(); k++) {
                    if(pListIdNodeTemp->at(k) == pVertex->id) {
                        pListIdNodeTemp->erase(pListIdNodeTemp->begin()+k);
                        break;
                    }
                }
                if(mapDegreeToListIdNode.find(pVertex->degree) == mapDegreeToListIdNode.end()) {
                    mapDegreeToListIdNode[pVertex->degree] = new vector<int>;
                }
                pListIdNodeTemp = mapDegreeToListIdNode[pVertex->degree];
                pListIdNodeTemp->push_back(pVertex->id);
                /// update neighbor
                pListIdNodeTemp = mapDegreeToListIdNode[pVertexNeighbor->degree-1];
                for(k=0; k<(int)pListIdNodeTemp->size(); k++) {
                    if(pListIdNodeTemp->at(k) == pVertexNeighbor->id) {
                        pListIdNodeTemp->erase(pListIdNodeTemp->begin()+k);
                        break;
                    }
                }
                if(mapDegreeToListIdNode.find(pVertexNeighbor->degree) == mapDegreeToListIdNode.end()) {
                    mapDegreeToListIdNode[pVertexNeighbor->degree] = new vector<int>;
                }
                pListIdNodeTemp = mapDegreeToListIdNode[pVertexNeighbor->degree];
                pListIdNodeTemp->push_back(pVertexNeighbor->id);
            }
            else if(typeRandomization == DEGREE_DISTRIBUTION_PRESERVING_EDGE_INTERCHANGE) { // not necesserily topology preserving
                ///check if there is at least one more node with the same degree
                bool flagFound = false, flagAddSuccess = false;
                vector<int> listIdNeighborAvailable;
                Vertex *pVertexNew, *pVertexNeighborNew;
                int nLoop = 0, nLoopInner=0;
                do {
                    if(nLoop > N_TRIAL_MAX) {
                        cerr << "max # of trial exceeded" << endl;
                        break;
                    }
                    nLoop++;
                    flagFound = false;
                    do {
                        iNode = rand() % listIdNode.size();
                        idNode = listIdNode[iNode];
                        pVertex = (*pMapNode)[idNode];
                        listIdNodeTemp = *(mapDegreeToListIdNode[pVertex->degree]);
                        if(listIdNodeTemp.size()<2) {
                            iNode = rand() % listIdNode.size();
                            idNode = listIdNode[iNode];
                            pVertex = (*pMapNode)[idNode];
                            continue;
                        }
                        pMapEdge = pVertex->pMapIdToEdge;
                        #ifdef CONTROL
                            if((unsigned int)pVertex->degree != pMapEdge->size()) cout << "Warning: degree inconsistency" << endl;
                        #endif //CONTROL
                        /// check neighbors that are in degree buckets that have more than one node
                        listIdNeighborAvailable.clear();
                        for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
                            pVertexNeighbor = (*pMapNode)[itEdge->first];
                            listIdNodeTemp = (*mapDegreeToListIdNode[pVertexNeighbor->degree]);
                            if(listIdNodeTemp.size()>1) {
                                listIdNeighborAvailable.push_back(itEdge->first);
                            }
                        }
                        if(listIdNeighborAvailable.size()<2) {
                            continue;
                        } else {
                            flagFound = true;
                        }
                    } while(flagFound==false);
                    iNeighbor = rand() % listIdNeighborAvailable.size();
                    idNeighbor = listIdNeighborAvailable[iNeighbor];
                    Edge *pEdge = NULL;
                    Data dataEdge, dataEdgeNeighbor;
                    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
                        if(itEdge->first == idNeighbor) {
                            pEdge = itEdge->second;
                            dataEdge = pEdge->data;
                            break;
                        }    
                    }
                    #ifdef CONTROL
                        if(itEdge == itEndEdge) cerr << "Warning: iteration inconsistency" << endl;
                        if(pVertex->id == itEdge->first) cerr << "Warning: unexpected self edge" << endl;
                    #endif //CONTROL
                    pVertexNeighbor = (*pMapNode)[idNeighbor];
                    listIdNodeTemp = *(mapDegreeToListIdNode[pVertex->degree]);
                    /// an exceptional case where these are the only two for that degree
                    if(!(pVertex->degree == pVertexNeighbor->degree and listIdNodeTemp.size()==2)) {
                        /// find another edge (pair of source and target) with same degrees 
                        //listIdNodeTemp = (*mapDegreeToListIdNode[pVertex->degree]);
                        nLoopInner=0;
                        do {
                            if(nLoopInner > N_TRIAL_MAX) {
                                break;
                            }
                            nLoopInner++;
                            do {
                                iNodeNew = rand() % listIdNodeTemp.size();
                                idNodeNew = listIdNodeTemp[iNodeNew];
                            } while(idNodeNew == idNode);
                            pVertexNew = (*pMapNode)[idNodeNew];
                            pMapEdge = pVertexNew->pMapIdToEdge;
                            /// check neighbors that has the same degree with idNeighbor
                            listIdNeighborAvailable.clear();
                            for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
                                pVertexNeighborNew = (*pMapNode)[itEdge->first];
                                if(pVertexNeighborNew->degree == pVertexNeighbor->degree) {
                                    listIdNeighborAvailable.push_back(itEdge->first);
                                }
                                pEdge = itEdge->second;
                                dataEdgeNeighbor = pEdge->data;
                            }
                        } while(listIdNeighborAvailable.size() < 1);
                        if(nLoopInner > N_TRIAL_MAX) {
                            nLoop++;
                            continue;
                        }
                        do {
                            iNeighborNew = rand() % listIdNeighborAvailable.size();
                            idNeighborNew = listIdNeighborAvailable[iNeighborNew];
                        } while(idNeighborNew == idNeighbor);
                        pVertexNeighborNew = (*pMapNode)[idNeighborNew];
                        /// some exceptions again
                        // rm 6 2 / rm 2 6 - add 6 6 add 2 2
                        if(pVertex->id == pVertexNeighborNew->id and pVertexNew->id == pVertexNeighbor->id) {
                            nLoop++;
                            continue;
                        }
                        pNetworkSampled->removeEdge(pVertex->id, pVertexNeighbor->id); //(idNode, idNeighbor);
                        //cout << "rm " << idNode << " " << idNeighbor << endl;
                        pNetworkSampled->removeEdge(pVertexNew->id, pVertexNeighborNew->id); 
                        //cout << "rm " << idNodeNew << " " << idNeighborNew << endl;
                        if(pNetworkSampled->addEdge(pVertex->id, pVertexNeighborNew->id, dataEdge) == false) {
                            /// restore removed edges
                            pNetworkSampled->addEdge(pVertex->id, pVertexNeighbor->id, dataEdge);
                            pNetworkSampled->addEdge(pVertexNew->id, pVertexNeighborNew->id, dataEdgeNeighbor);
                        } else if(pNetworkSampled->addEdge(pVertexNew->id, pVertexNeighbor->id, dataEdgeNeighbor) == false) {
                            pNetworkSampled->addEdge(pVertex->id, pVertexNeighbor->id, dataEdge);
                            pNetworkSampled->addEdge(pVertexNew->id, pVertexNeighborNew->id, dataEdgeNeighbor);
                        } else {
                            //cout << "add " << pVertex->id << " " << pVertexNeighborNew->id << endl;
                            //cout << "add " << idNode << " - " << pVertex->id << " / " << idNeighborNew << " - " << pVertexNeighborNew->id << endl;
                            //cout << "add " << pVertexNew->id << " " << pVertexNeighbor->id <<endl;
                            //cout << "add " << idNodeNew << " - " << pVertexNew->id << " / " << idNeighbor << " - " << pVertexNeighbor->id <<endl;
                            flagAddSuccess = true;
                        }
                    }
                } while(flagAddSuccess == false);
                /*
                int nLoop = 0;
                do {
                    if(pVertex->degree == 0 or mapDegreeToListIdNode.find(pVertex->degree) == mapDegreeToListIdNode.end()) {
                        idNodeNew = pVertex->id;
                        // select one one with one minus degree of idNeighbor
                        if(mapDegreeToListIdNode.find(pVertexNeighbor->degree) == mapDegreeToListIdNode.end()) {
                            idNeighborNew = pVertexNeighbor->id;
                        } else {
                            listIdNodeTemp = *(mapDegreeToListIdNode[pVertexNeighbor->degree]);
                            if(listIdNodeTemp.size() == 1) {
                                iNeighborNew = 0;
                            } else {
                                do {
                                    iNeighborNew = rand() % listIdNodeTemp.size();
                                } while(iNeighborNew == iNeighbor);
                            }
                            idNeighborNew = listIdNodeTemp[iNeighborNew];
                            if(idNeighborNew == idNodeNew) {
                                idNeighborNew = pVertexNeighbor->id;
                            }
                        }
                    } else if(pVertexNeighbor->degree == 0 or mapDegreeToListIdNode.find(pVertexNeighbor->degree) == mapDegreeToListIdNode.end()) {
                        idNodeNew = pVertexNeighbor->id;
                        // select one node with the one minus degree of idNode
                        if(mapDegreeToListIdNode.find(pVertex->degree) == mapDegreeToListIdNode.end()) {
                            idNeighborNew = pVertex->id;
                        } else {
                            listIdNodeTemp = *(mapDegreeToListIdNode[pVertex->degree]);
                            if(listIdNodeTemp.size() == 1) {
                                iNeighborNew = 0;
                            } else {
                                do {
                                    iNeighborNew = rand() % listIdNodeTemp.size();
                                } while(iNeighborNew == iNode);
                            }
                            idNeighborNew = listIdNodeTemp[iNeighborNew];
                            if(idNeighborNew == idNodeNew) {
                                idNeighborNew = pVertex->id;
                            }
                        }
                    } else {
                        // select one node with the one minus degree of idNode
                        listIdNodeTemp = *(mapDegreeToListIdNode[pVertex->degree]);
                        if(listIdNodeTemp.size() == 1) {
                            iNodeNew = 0;
                        } else {
                            do {
                                iNodeNew = rand() % listIdNodeTemp.size();
                            } while(iNodeNew == iNode);
                        }
                        idNodeNew = listIdNodeTemp[iNodeNew];
                        // select one one with one minus degree of idNeighbor
                        listIdNodeTemp = *(mapDegreeToListIdNode[pVertexNeighbor->degree]);
                        if(listIdNodeTemp.size() == 1) {
                            iNeighborNew = 0;
                        } else {
                            do {
                                iNeighborNew = rand() % listIdNodeTemp.size();
                            } while(iNeighborNew == iNeighbor);
                        }
                        idNeighborNew = listIdNodeTemp[iNeighborNew];
                        if(idNeighborNew == idNodeNew) {
                            idNodeNew = pVertex->id;
                            idNeighborNew = pVertexNeighbor->id;
                        }
                    }
                    nLoop++;
                    if(nLoop > N_TRIAL_MAX) {
                        cout << "max # of trial exceeded" << endl;
                        break;
                    }
                    cout << "add " << idNodeNew << " " << idNeighborNew << endl;
                } while(pNetworkSampled->addEdge(idNodeNew, idNeighborNew, dataEdge) == false); 
                if(nLoop > N_TRIAL_MAX) {
                    pNetworkSampled->addEdge(idNode, idNeighbor, pEdge->data); 
                }
                */
            }
            else if(typeRandomization == TOPOLOGY_PRESERVING) {
        		do {
            		iNodeNew = rand() % listIdNode.size();
            	} while(iNodeNew == iNode);
        		idNodeNew = listIdNode[iNodeNew];
        		//cout << "changing " << (pNetworkSampled->mapIdToName)[idNode] << " (" << idNode << ") " << (pNetworkSampled->mapIdToName)[idNodeNew] << " (" << idNodeNew << ")" << endl;
        		pNetworkSampled->swapInMapNewId(idNode, idNodeNew);
        		//pNetworkSampled->convertMapNewIdToMapId();
        		//printGraph(*pNetworkSampled, oss.str());
           	} else if(typeRandomization == RANDOM_PRESERVING_ONE_DEGREED_NODES){
	        	while(pVertex->degree < 2) {
	        		iNode = rand() % listIdNode.size();
		        	idNode = listIdNode[iNode];
		        	pVertex = (*pMapNode)[idNode];
	        	}
            	pMapEdge = pVertex->pMapIdToEdge;
				#ifdef CONTROL
            		if((unsigned int)pVertex->degree != pMapEdge->size()) cout << "Warning: degree inconsistency" << endl;
				#endif //CONTROL
        		iNeighbor = rand() % pVertex->degree;
                Edge *pEdge = NULL;
        		for(k=0, itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge, k++) {
        			if(k == iNeighbor) {
                        pEdge = itEdge->second;
        				break;
                    }
        		}
				#ifdef CONTROL
        			if(itEdge == itEndEdge) cerr << "Warning: iteration inconsistency" << endl;
				#endif //CONTROL
    			pNetworkSampled->removeEdge(pVertex->id, itEdge->first);
    			//cout << "remove: " << pVertex->id << " " << itEdge->first << endl;
        		do {
	            	do {
	            		iNodeNew = rand() % listIdNode.size();
	            	} while(iNodeNew == iNode);
	            	do {
	            		iNeighborNew = rand() % listIdNode.size();
	            	} while(iNeighborNew == iNodeNew);
	            	//cout << "add: " << listIdNode[iNodeNew] << " " << listIdNode[iNeighborNew] << endl;
        		} while(pNetworkSampled->addEdge(listIdNode[iNodeNew], listIdNode[iNeighborNew], pEdge->data) == false); 
        	}        	 
        }
        if(pNetworkSampled->isSampled()) { //if(flagPreserveTopology) {
        	//cout << "here" << endl;
        	pNetworkSampled->convertMapNewIdToMapId();
        }
        listNetworkSampled.push_back(pNetworkSampled);
    } 
    // clean temporaly used helper storages
    if(typeRandomization == TOPOLOGY_PRESERVING_SAME_DEGREE_INTERCHANGE) {
        hash_map<int, vector<int> *>::iterator itDegree;
        for(itDegree = mapDegreeToListIdNode.begin(); itDegree != mapDegreeToListIdNode.end(); ++itDegree) {
            delete itDegree->second;
        }
    }
}

/// Z score calculation for nodes (no hash map, direct summation)
float Netzcore::calculateNodeZScore(int vId, bool flagUseEdgeReliabilityScore, bool flagScaleWithNodeDegree) {
    unsigned int i = 0, nSample;//, j=0;
    MapIntToVertex *pMapNode;
    MapIntToEdge *pMapEdge;
    MapIntToInt *pMapIdNew;
    Vertex *pVertex, *pNeighbor;
    Edge *pEdge;
    MapIntToString mapName = network.mapIdToName;
    float sumScoreNeighbor = 0.0, sumScore = 0.0, sumScoreSquare = 0.0, zScore = 0.0;
    float sumScore_square = 0.0;
    MapIntToEdge::iterator itEdge, itEndEdge;
    Graph *g;
	
    nSample = listNetworkSampled.size(); //N_SAMPLE;
    for(i=0; i<nSample; i++) {
        g = listNetworkSampled[i];
        pMapNode = g->pMapIdToNode;
        MapStringToInt mapNameToIdSampled =  g->mapNameToId; 
        if(listNetworkSampled[i]->isSampled()) {
                pMapIdNew = g->pMapIdNewToId; // previously pMapIdToIdNew
                pVertex = (*pMapNode)[(*pMapIdNew)[vId]]; 
        } else {	
                pVertex = (*pMapNode)[mapNameToIdSampled[mapName[vId].c_str()]]; // previously pVertex = (*pMapNode)[vId];
        }
        pMapEdge = pVertex->pMapIdToEdge;
        sumScoreNeighbor = 0.0;
        for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            pNeighbor = (*pMapNode)[itEdge->first];
            if(flagUseEdgeReliabilityScore) {
                sumScoreNeighbor += pEdge->data.score * pNeighbor->data.score;
            } else {
                sumScoreNeighbor += pNeighbor->data.score;
            }
        }
        sumScore += sumScoreNeighbor;
        sumScoreSquare += sq(sumScoreNeighbor);
    }
    sumScore /= nSample;
    sumScoreSquare /= nSample;
    sumScore_square = sq(sumScore);
    pMapNode = network.pMapIdToNode;
    pVertex = (*pMapNode)[vId];
    pMapEdge = pVertex->pMapIdToEdge;
    sumScoreNeighbor = 0.0;
    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
        pEdge = itEdge->second;
        pNeighbor = (*pMapNode)[itEdge->first];
        if(flagUseEdgeReliabilityScore) {
            sumScoreNeighbor += pEdge->data.score * pNeighbor->data.score;
        } else {
            sumScoreNeighbor += pNeighbor->data.score;
        }
        //cerr << mapName[pNeighbor->id] << " " << pNeighbor->data.score << endl;
    }
    if(sumScore_square > sumScoreSquare) {
        /// problem in z scoring E(x2) < E2(x)
        if(sumScore_square - sumScoreSquare >  MIN_SCORE_SHIFT) {
            cerr << "Warning: meanSquare is less than squareMean in Z score calculation ";
            cerr << "s: " << sumScoreNeighbor << " E[x^2]: " << sumScoreSquare << " E[x]^2: " << sumScore_square << " E[x]: " << sumScore << endl;
        }
        zScore = 0.0;
    } else {
        zScore = (sumScoreNeighbor - sumScore);
        if(sumScore_square == sumScoreSquare) {
            // zScore /= MIN_SCORE_SHIFT;
            if(zScore == 0) zScore = 0;
            else if (zScore > 0) zScore = INFINITY; 
            else if (zScore < 0) zScore = -INFINITY;
        } else {
            zScore /= float(sqrt(sumScoreSquare - sumScore_square));
        }
        if(flagScaleWithNodeDegree) {
            zScore *= float(1/sqrt(pVertex->degree));
        } 
        /// in order to avoid inf - inf = nan confusion
        //if(isinf(zScore) == 1 or isinf(zScore) == 1) { //if(isnan(zScore)) {
        if(zScore > MAX_Z_SCORE) {
        //    cerr << "Warning: resetting overflowing z score " << zScore << endl;
            zScore = MAX_Z_SCORE;
        } else if(zScore < -MAX_Z_SCORE) {
        //    cerr << "Warning: resetting overflowing z score " << zScore << endl;
            zScore = -MAX_Z_SCORE;
        }
        if (flagDebug) { 
            cout << mapName[vId] << ": ( " << sumScoreNeighbor << " - " << sumScore <<  " ) / " <<  float(sqrt(sumScoreSquare - sq(sumScore))) << " z: "<< zScore << endl;
        }
    }
    
    return zScore;
}

/*
/// In Z score calculation floats are stored in hash with PRECISION so all the calculations must be based on the same precision to avoid ex2-e2x < 0
float Netzcore::calculateNodeZScore(int vId, bool flagUseEdgeReliabilityScore, bool flagScaleWithNodeDegree) {
	hash_map<int, int>::iterator itHistogram, itEndHistogram;
	hash_map<int, int> mapHistogramScore;
	unsigned int i = 0;//, j=0;
	MapIntToVertex *pMapNode;
	MapIntToEdge *pMapEdge;
	MapIntToInt *pMapId;
	Vertex *pVertex, *pNeighbor;
	Edge *pEdge;
	float sumScoreNeighbor = PRECISION_INIT, sumScore = PRECISION_INIT, sumScoreSquare = PRECISION_INIT, scoreTemp=PRECISION_INIT, zScore = PRECISION_INIT;
    float sumScoreSquare_sqrt = PRECISION_INIT, sumScore_square = PRECISION_INIT;
	MapIntToEdge::iterator itEdge, itEndEdge;
	int intScore = 0;
	
    sumScore = PRECISION_INIT;
    sumScoreSquare = PRECISION_INIT;
    scoreTemp = PRECISION_INIT;
	for(i=0; i<listNetworkSampled.size(); i++) {
		pMapNode = listNetworkSampled[i]->pMapIdToNode;
		if(listNetworkSampled[i]->isSampled()) {
			pMapId = listNetworkSampled[i]->pMapIdToIdNew;
			pVertex = (*pMapNode)[(*pMapId)[vId]]; 
		} else {
			pVertex = (*pMapNode)[vId];	
		}
		pMapEdge = pVertex->pMapIdToEdge;
		sumScoreNeighbor = PRECISION_INIT;
		for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
        	pEdge = itEdge->second;
            pNeighbor = (*pMapNode)[itEdge->first];
            if(flagUseEdgeReliabilityScore) {
    	    	sumScoreNeighbor += pEdge->data.score * pNeighbor->data.score;
            } else {
    	    	sumScoreNeighbor += pNeighbor->data.score;
            }
		}
		//scoreTemp = sumScoreNeighbor;
        intScore = convertScoreFloatToInt(sumScoreNeighbor);
		//cout << "h: " << sumScoreNeighbor << " " << intScore << endl;
		if(mapHistogramScore.find(intScore) == mapHistogramScore.end()) { 
			mapHistogramScore[intScore] = 1;
		}
		else {
			mapHistogramScore[intScore] = mapHistogramScore[intScore] + 1;
		}
	}
    //sumScore = 0.000000f;
    //sumScoreSquare = 0.000000f;
    //scoreTemp = 0.000000f;
	for(i=0, itHistogram = mapHistogramScore.begin(), itEndHistogram = mapHistogramScore.end(); itHistogram != itEndHistogram and i < mapHistogramScore.size(); ++itHistogram, i++) {
        // itHistogram is corrupted - most probably due to copy / could be that when pointer not a problem
        scoreTemp = convertScoreIntToFloat(itHistogram->first);
		//cout << " 1.: " << itHistogram->first << " 2.: " << itHistogram->second << " s: " << scoreTemp << endl;
		sumScore += scoreTemp * itHistogram->second;
		//sumScore += convertScoreIntToFloat(convertScoreFloatToInt(scoreTemp * float(itHistogram->second)));
		//sumScore = convertScoreIntToFloat(convertScoreFloatToInt(sumScore));
		sumScoreSquare += sq(scoreTemp) * itHistogram->second;
		//sumScoreSquare += convertScoreIntToFloat(convertScoreFloatToInt(sq(scoreTemp))) * float(itHistogram->second);
		//sumScoreSquare = convertScoreIntToFloat(convertScoreFloatToInt(sumScoreSquare));
		//cout << i << " n " << itHistogram->second << " nT " << listNetworkSampled.size() << " s " << scoreTemp << " sS " << sumScore << " sSS " << sumScoreSquare << endl;
	}
    sumScore /= listNetworkSampled.size();
    //sumScore = convertScoreIntToFloat(convertScoreFloatToInt(sumScore));
    sumScoreSquare /= listNetworkSampled.size();
    //sumScoreSquare = convertScoreIntToFloat(convertScoreFloatToInt(sumScoreSquare));
	pMapNode = network.pMapIdToNode;
	pVertex = (*pMapNode)[vId];
	sumScoreNeighbor = PRECISION_INIT;
	for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
    	pEdge = itEdge->second;
        pNeighbor = (*pMapNode)[itEdge->first];
        if(flagUseEdgeReliabilityScore) {
    	    sumScoreNeighbor += pEdge->data.score * pNeighbor->data.score;
        } else {
        	sumScoreNeighbor += pNeighbor->data.score;
        }
	}
    /// when histogram has one element zScore = 0 since sq(sumScore) == sumScoreSquare
    if(mapHistogramScore.size() == 1) {
	   	//cerr << "histogram contains 1 element - s: " << sumScoreNeighbor << " E[x^2]: " << sumScoreSquare << " E[x]^2: " << sumScore_square << " E[x]: " << sumScore << endl;
        zScore = 0.0;
    } else {
        sumScoreSquare_sqrt = sqrt(sumScoreSquare); 
        sumScore_square = sq(sumScore); //sq(convertScoreIntToFloat(convertScoreFloatToInt(sumScore)));
        /// float comparision equal if within epsilon == MIN_SCORE_SHIFT
        if(abs(sumScoreSquare_sqrt - sumScore) < MIN_SCORE_SHIFT) {
            zScore = 0.0;
        /// problem in z scoring E(x2) < E2(x)
        } else if(sumScoreSquare_sqrt < sumScore) {
                cerr << "Warning: meanSquare is less then squareMean in Z score calculation" << " - ";
                cerr << "s: " << sumScoreNeighbor << " E[x^2]: " << sumScoreSquare << " E[x]^2: " << sumScore_square << " E[x]: " << sumScore << " sqrt(E[x^2]): " << sumScoreSquare_sqrt << endl;
                for(i=0, itHistogram = mapHistogramScore.begin(), itEndHistogram = mapHistogramScore.end(); itHistogram != itEndHistogram and i < mapHistogramScore.size(); ++itHistogram, i++) {
                    cerr << i << " " << mapHistogramScore.size() << " / " << itHistogram->first << " : " << itHistogram->second << endl;
                }
                zScore = 0.0;
	        } else {
	            cout << vId << ": " << sumScoreNeighbor << " - " << sumScore <<  " / " <<  float(sqrt(sumScoreSquare - sq(sumScore))) << endl;
	            //cout << float(sqrt(pVertex->degree)) * (sumScoreNeighbor - sumScore) / float(sqrt(sumScoreSquare - sq(sumScore)));
                //cout << " " << pVertex->degree << " " << sumScoreNeighbor << " " << sumScore << " " << sumScoreSquare << endl;
                if(flagScaleWithNodeDegree) {
    	            zScore = float(1/sqrt(pVertex->degree)) * (sumScoreNeighbor - sumScore) / float(sqrt(sumScoreSquare - sumScore_square));
                } else {
                    zScore = (sumScoreNeighbor - sumScore) / float(sqrt(sumScoreSquare - sumScore_square));
                }
                //cerr << "z score: " << zScore << endl;
                /// in order to avoid inf - inf = nan confusion
                //if(isinf(zScore) == 1 or isinf(zScore) == 1) { //if(isnan(zScore)) {
                if(zScore > MAX_Z_SCORE) {
                    cerr << "Warning: resetting overflowing z score " << zScore << endl;
                    zScore = MAX_Z_SCORE;
                } else if(zScore < -MAX_Z_SCORE) {
                    cerr << "Warning: resetting overflowing z score " << zScore << endl;
                    zScore = -MAX_Z_SCORE;
                }
            }
	}
    
    //if(isnan(sumScoreNeighbor)) {
    //    cerr << "Warning: sum of neighbor scores are nan for node " << pVertex->id;
    //    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
    //        pNeighbor = (*pMapNode)[itEdge->first];
    //        cerr << " " << pNeighbor->data.score;
    //    }
    //    cerr << endl;
    //    sumScoreNeighbor = 0.0;
    //}

	return zScore;
}
*/

// must be incorporated in hash_map<float> in hash.h
int Netzcore::convertScoreFloatToInt(float scoreFloat, unsigned int precisionCoefficient) {
    //for(unsigned int j=0; j<precision; j++) {
    //    scoreFloat *= 10;
    //}
    return (int)((scoreFloat * precisionCoefficient) + 0.5);
}

float Netzcore::convertScoreIntToFloat(int scoreInt, unsigned int precisionCoefficient) {
    //float scoreFloat = (float)scoreInt;
    //for(unsigned int j=0; j<precision; j++) {
    //    scoreFloat /= 10;
    //}
    //return scoreFloat;
    return float(scoreInt) / precisionCoefficient;
}

void Netzcore::deleteSampledNetworks() {
    for(unsigned int i = 0; i<listNetworkSampled.size(); i++) {
        delete listNetworkSampled[i];
    }
    listNetworkSampled.clear();
}

void Netzcore::updateNetwork(int typeTransfer, int nIteration, bool flagUseEdgeReliabilityScore, bool flagNoZScoring, bool flagAccumulate, bool flagScaleZScoreWithNodeDegree) { 
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToEdge::iterator itEdge, itEndEdge;
    MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
    MapIntToEdge *pMapEdge;
    Edge *pEdge;
    hash_map<int,void*> mapIdProcessed;
    pair<float, int> tempPair(1.0, 0);
	
    if(nIteration == 0 and flagNoZScoring) {
        /// Initialize message arrays
        for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
            pVertex = itNode->second;
            // now pair stores (edge_multiplier, n_iteration) with above initial values
            //tempPair.first = pVertex->data.score; 
            (*pVertex->pMapIdToScore)[pVertex->id] = tempPair;
        }
        return; //! for -N netscore mode just initialize and do nothing else
    }
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
        updateNodeScore(typeTransfer, itNode->first, nIteration, flagUseEdgeReliabilityScore, flagNoZScoring, flagAccumulate, flagScaleZScoreWithNodeDegree); 
    }
    /// Adjust z-scores: shift them so that they are bigger than 0 and seed proteins have at least highest possible z-score
    if(!flagNoZScoring) {
        //cout << "zMin: " << minZScore << " zMax: " << maxZScore << endl;
        for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
                pVertex = itNode->second;
                pVertex->data.scoreUpdated += (float)abs(minZScore);
                //if(nIteration == 0 and mapSource.find(pVertex->id) != mapSource.end()) {
                //    pVertex->data.scoreUpdated += (maxZScore-minZScore);
                //}
                pVertex->data.scoreUpdated += pVertex->data.score;
                if(nIteration == 0) {
                    // required for later - during score calculation of a node, allways initial score of self is taken into consideration
                    pVertex->data.scoreInitial = pVertex->data.scoreUpdated;
                }
        }
        /// Initialize message arrays
        if(nIteration == 0) {
            for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
                pVertex = itNode->second;
                //tempPair.first = pVertex->data.scoreUpdated;
                (*pVertex->pMapIdToScore)[pVertex->id] = tempPair;
            }
        }
    }
    if(flagDebug)	
            printNetwork("After updateNodeScore");
    
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

float Netzcore::transferScore(int typeTransfer, float score, float a, float b) {
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

void Netzcore::scaleNodeScores() {
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex;
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        pVertex->data.score /= maxScore; 
    }
    return;
}


float Netzcore::calculateErrorAndUpdateScores(bool flagResetSeedScoresToInitial) {
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

void Netzcore::updateNodeScore(int typeTransfer, int vId, int nIteration, bool flagUseEdgeReliabilityScore, bool flagNoZScoring, bool flagAccumulate, bool flagScaleZScoreWithNodeDegree) {
    MapIntToEdge::iterator itEdge, itEndEdge;
    MapIntToFloatInt::iterator itScore, itEndScore;
    MapIntToVertex *pMapNode = network.pMapIdToNode;
    Vertex *pVertex = (*pMapNode)[vId], *pNeighbor;
    MapIntToEdge *pMapEdge = pVertex->pMapIdToEdge;
    Edge *pEdge; //, *pEdgeNeighbor;
    MapIntToFloatInt *pMapNeighborScore = pVertex->pMapIdToScore, *pMapNeighborScoreNeighbor;
    //float sumScoreNeighbor = PRECISION_INIT, scoreCalculated = PRECISION_INIT, tempScore = PRECISION_INIT;
    float sumScoreNeighbor = 0.0, scoreCalculated = 0.0, tempScore = 0.0;
    float i=0.0;
    if(flagVerbose)
        cerr << "-Checking node: " << vId << endl;
    if(flagNoZScoring) {
        for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pNeighbor = (*pMapNode)[itEdge->first];
            pEdge = itEdge->second;
            pMapNeighborScoreNeighbor = pNeighbor->pMapIdToScore;
            /*
            if(flagUseEdgeReliabilityScore) {
                ///sumScoreNeighbor += pEdge->data.score * pNeighbor->data.score; // / pVertex->degree;
                //tempScore = pNeighbor->data.score - pEdge->data.scoreAccumulatedLast;
                //tempScore *= pEdge->data.score;
                tempScore = pNeighbor->data.score * pEdge->data.score;
            } else {
                //tempScore = pNeighbor->data.score - pEdge->data.scoreAccumulatedLast;
                tempScore = pNeighbor->data.score;
                //cout << "_" << itEdge->first << " t: " << tempScore << endl; //<< " last: " << pEdge->data.scoreAccumulatedLast << " lastUp: " << pEdge->data.scoreAccumulatedLastUpdated << endl;
            }
            if(nIteration != 0) {
                tempScore /= sq(nIteration);
            }
            if(pMapNeighborScore->find(pNeighbor->id) == pMapNeighborScore->end()) { // if you would use it in for iterators would be different since pMapNeighborScore->begin() gives an iterator of the object but consequent it++ does not change iterator state
                (*pMapNeighborScore)[pNeighbor->id] = pair<float, int>(tempScore, nIteration); // tempScore
            } else {
                (*pMapNeighborScore)[pNeighbor->id].first += tempScore; 
                (*pMapNeighborScore)[pNeighbor->id].second = nIteration;
            }
            for(itScore = pMapNeighborScoreNeighbor->begin(), itEndScore = pMapNeighborScoreNeighbor->end(); itScore != itEndScore; ++itScore) {
                if(itScore->first != vId) {
                    if(nIteration != 0) {
                        (*pMapNeighborScore)[itScore->first].first += itScore->second.first/sq(nIteration);
                        (*pMapNeighborScore)[itScore->first].second = nIteration;
                    } else {
                        (*pMapNeighborScore)[itScore->first].first += itScore->second.first;
                        (*pMapNeighborScore)[itScore->first].second = nIteration;
                    }
                }
            }
            */
            if(flagVerbose)
                cerr << "--Messages from neighbor: " << itEdge->first << endl; 
            for(itScore = pMapNeighborScoreNeighbor->begin(), itEndScore = pMapNeighborScoreNeighbor->end(); itScore != itEndScore; ++itScore) {
                if(itScore->second.second < nIteration) {
                    if(pMapNeighborScore->find(itScore->first) == pMapNeighborScore->end()) {
                        // now tempScore corresponds to edge_multiplier rather than score
                        tempScore = (itScore->second).first;
                        if(flagUseEdgeReliabilityScore) {
                            tempScore *= pEdge->data.score;
                        }
                        //if(nIteration != 0) {
                        //    tempScore /= sq(nIteration);
                        //}
                        (*pMapNeighborScore)[itScore->first].first = tempScore;
                        (*pMapNeighborScore)[itScore->first].second = nIteration;
                        if(flagVerbose)
                            cerr << "---message of " << itScore->first << " score: " << tempScore << endl; 
                    }
                }
            }

            //sumScoreNeighbor += tempScore;
            //pEdgeNeighbor = (*(pNeighbor->pMapIdToEdge))[vId];
            //pEdgeNeighbor->data.scoreAccumulatedLastUpdated = tempScore;
        }
        i=0;
        for(itScore = pMapNeighborScore->begin(), itEndScore = pMapNeighborScore->end(); itScore != itEndScore; ++itScore) {
            // now using up-to-date scores of nodes at each step
            pNeighbor = (*pMapNode)[itScore->first];
            if(itScore->second.second != 0) { //nIteration
                ////sumScoreNeighbor += (itScore->second).first / itScore->second.second; //sq(itScore->second.second); //nIteration
                //sumScoreNeighbor +=  pNeighbor->data.score * (itScore->second).first  / sq(itScore->second.second); 
                sumScoreNeighbor +=  pNeighbor->data.scoreInitial * (itScore->second).first  / sq(itScore->second.second); 
                i += 1.0 / sq(itScore->second.second);
            } else {
                // consider initial self score (rather than updated score) while calculating score for that node
                //sumScoreNeighbor += pNeighbor->data.score * (itScore->second).first;
                sumScoreNeighbor += pNeighbor->data.scoreInitial * (itScore->second).first;
                i += 1.0;
            }
        } 
        sumScoreNeighbor /= i;
        scoreCalculated = transferScore(typeTransfer, sumScoreNeighbor);
        pVertex->data.scoreUpdated = scoreCalculated;
        if(flagAccumulate) {
            ///pVertex->data.scoreUpdated = pVertex->data.score + scoreCalculated;
            ///pVertex->data.scoreUpdated = pVertex->data.scoreInitial + scoreCalculated;
            ///pVertex->data.scoreUpdated += pVertex->data.scoreInitial; 
            ///pVertex->data.scoreUpdated += pVertex->data.score - pVertex->data.scoreAccumulatedLast; 
            //pVertex->data.scoreUpdated += pVertex->data.score - pVertex->data.scoreAccumulatedLast/2; 
            //pVertex->data.scoreAccumulatedLastUpdated = scoreCalculated;
            ///cout << "acc: " << scoreCalculated << endl;
            pVertex->data.scoreUpdated += pVertex->data.score - pVertex->data.scoreAccumulatedLast; 
            pVertex->data.scoreAccumulatedLast = scoreCalculated;
        } 
        if(pVertex->data.scoreUpdated < MIN_SCORE) {
            pVertex->data.scoreUpdated = MIN_SCORE;
        } else if(pVertex->data.scoreUpdated > MAX_SCORE) {
            pVertex->data.scoreUpdated = MAX_SCORE;
        }
    } else {
        scoreCalculated = transferScore(typeTransfer, calculateNodeZScore(vId, flagUseEdgeReliabilityScore, flagScaleZScoreWithNodeDegree));
        pVertex->data.scoreUpdated = scoreCalculated;
        if(scoreCalculated<minZScore) {
            minZScore = scoreCalculated;
        }
        if(scoreCalculated>maxZScore) {
            maxZScore = scoreCalculated;
        }
    }
    //cout << vId << ": " << scoreCalculated << endl;
    if(flagVerbose)
        cerr << "-Total score of node: "<< pVertex->id << " at the end of this iteration: " << pVertex->data.scoreUpdated << endl;
}

void Netzcore::updateEdgeScore(int typeTransfer, int idSource, int idTarget) {
	/*
	MapIntToVertex *pMapNode = network.pMapIdToNode;
	Vertex *pSource = (*pMapNode)[idSource], *pTarget = (*pMapNode)[idTarget];
	MapIntToEdge *pMapEdgeSource = pSource->pMapIdToEdge, *pMapEdgeTarget = pTarget->pMapIdToEdge;
	Edge *pEdgeSource = (*pMapEdgeSource)[idTarget], *pEdgeTarget = (*pMapEdgeTarget)[idSource];
	//pEdgeSource->data.scoreUpdated = transferScore(typeTransfer, );
	*/
}

void Netzcore::updateSampledNetworkNodeScores() {
    MapIntToVertex::iterator itNode, itEndNode, itNodeSampled, itEndNodeSampled;
    MapIntToVertex *pMapNode = network.pMapIdToNode, *pMapNodeSampled;
    MapIntToString mapName = network.mapIdToName;
    MapIntToInt *pMapIdNewSampled;
    Vertex *pVertex, *pVertexSampled;

    for(unsigned int i = 0; i<listNetworkSampled.size(); i++) {
        Graph *g = listNetworkSampled[i];
    	pMapNodeSampled = g->pMapIdToNode;
        pMapIdNewSampled =  g->pMapIdNewToId; 
        MapStringToInt mapNameToIdSampled =  g->mapNameToId; 
        //if(i==0) cerr << "isSampled: " << g->isSampled() << endl;
        /*
    	itNode = pMapNode->begin();
    	itEndNode = pMapNode->end();
    	itNodeSampled = pMapNodeSampled->begin();
    	itEndNodeSampled = pMapNodeSampled->end();
    	for(; itNode != itEndNode and itNodeSampled != itEndNodeSampled; ++itNode, ++itNodeSampled) {
            #ifdef CONTROL
                if(itNode->first != itNodeSampled->first) cout << "Warning: iteration inconsistency" << endl;
            #endif //CONTROL
            pVertex = itNode->second;
            pVertexSampled = itNodeSampled->second;
            pVertexSampled->data.score = pVertex->data.score;
        }
        */
        for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
            pVertex = itNode->second;
            if(g->isSampled()) {
                pVertexSampled = (*pMapNodeSampled)[(*pMapIdNewSampled)[pVertex->id]];
            } else {
                pVertexSampled = (*pMapNodeSampled)[mapNameToIdSampled[mapName[pVertex->id].c_str()]];
                //if(i == 0) {
                //    cerr << pVertex->id << ": " << pVertex->data.score << " sampled " << pVertexSampled->id << ": " << pVertexSampled->data.score << endl;
                //}
            }
            pVertexSampled->data.score = pVertex->data.score;
	}
    }
}

/*
void Netzcore::resetSampledNetworksNodeScores() {
    for(unsigned int i = 0; i<listNetworkSampled.size(); i++) {
        resetSampledNetworkNodeScores(listNetworkSampled[i]);
    }
    return;
}

void Netzcore::resetSampledNetworkNodeScores(Graph *g) {
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToVertex *pMapNode = network.pMapIdToNode, *pMapNodeSampled = g->pMapIdToNode; 
    MapIntToInt *pMapIdNewSampled = g->pMapIdNewToId; 
    Vertex *pVertex, *pVertexSampled;

    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        if(g->isSampled()) {
            pVertexSampled = (*pMapNodeSampled)[(*pMapIdNewSampled)[pVertex->id]];
            pVertexSampled->data.score = pVertex->data.score;
        } else {
            pVertexSampled = (*pMapNodeSampled)[pVertex->id];
            pVertexSampled->data.score = pVertex->data.score;
        }
    }
    return;
}
*/

void Netzcore::writeSampledNetworks() {
    ostringstream oss;
    for(unsigned int i = 0; i<listNetworkSampled.size(); i++) {
        oss.str("");
        oss << RANDOM_SAMPLING_DIRECTORY << "/" << i << ".txt";
        outputGraph(listNetworkSampled[i], oss.str());
    }
    return;
}

void Netzcore::readSampledNetworks(int nSample) {
    ostringstream oss;
    for(int i = 0; i<nSample; i++) {
        oss.str("");
        oss << RANDOM_SAMPLING_DIRECTORY << "/" << i << ".txt";
        listNetworkSampled.push_back(readSampledNetwork(oss.str()));
    }
    return;
}

Graph * Netzcore::readSampledNetwork(string fileName) {
        Graph * pNetworkSampled = new Graph();
        fstream file;
        file.open(fileName.c_str(), ios::in);
        string name1, name2;
        float score1, score2, affinity;
        if(!file) cerr << "Warning: file can not be opened " << fileName << endl;
        while(file >> name1 >> score1 >> name2 >> score2 >> affinity) {
            pNetworkSampled->addVertex(name1, score1);
            pNetworkSampled->addVertex(name2, score2);
            pNetworkSampled->addEdge(name1, name2, affinity);
        }
        file.close();
        return pNetworkSampled;
}

//void Netzcore::outputGraph(Graph const & g, string fileName) {
void Netzcore::outputGraph(Graph *g, string fileName) {
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

void Netzcore::printNetwork(string explanation) {
    printGraph(network, explanation);
}

void Netzcore::printGraph(Graph const & g, string explanation) {
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

void Netzcore::printSampledNetworks() {
	ostringstream oss;	
	for(unsigned int i = 0; i<listNetworkSampled.size(); i++) {
		oss.str("");
		oss << i;
		printGraph(*(listNetworkSampled[i]), oss.str());
	}
}
