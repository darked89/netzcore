
void Netzcore::updateNodeScore(int vId, int nIteration) {
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

void Netzcore::updateEdgeScore(int idSource, int idTarget) {
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

