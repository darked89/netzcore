#include "Netween.hpp"

#include <iostream>
#include <vector>

using namespace boost;
using namespace std;

Netween::Netween() 
{
    flagAccumulateToInitialNodeScore = false;
    flagVerbose = false;
}

Netween::Netween(string fileNode, string fileEdge, string fileOutput, bool fAccumulateToInitialNodeScore, bool fVerbose)
{
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    flagVerbose = fVerbose;
    outputFile = fileOutput;
    network.loadNodes(fileNode); 
    network.loadEdges(fileEdge);//, true); 
    VertexIterator it, itEnd;
    for(tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	//! Using default non seed score assumption
	if(getNetwork().getVertexScore(*it) > 0.01)
	    setSeed.insert(*it);
    }
}


Netween::~Netween() {
}


void Netween::run() 
{ 
    VertexIterator it, itEnd;
    unordered_set<Vertex> setIncluded;
    unordered_map<Vertex, int> mapVertexToLocalCount;
    unordered_map<Vertex, int> mapVertexToGlobalCount;

    // Initialize node involvement counts
    for(tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	mapVertexToLocalCount[*it] = 0;
	mapVertexToGlobalCount[*it] = 0;
    }

    //map<Vertex, float> mapDistance;
    map<Vertex, vector<Vertex> > mapPredecessors;

    unordered_set<Vertex>::iterator seedIt, seedItEnd;
    unordered_set<Vertex>::iterator setIt, setItEnd;
    map<Vertex, vector<Vertex> >::iterator preIt, preItEnd;
    Vertex v, v_prev;

    // For each seed find shortest paths: Ps,t (s € seeds, t € nodes)
    for(seedIt = setSeed.begin(), seedItEnd = setSeed.end(); seedIt != seedItEnd; ++seedIt) 
    {
	if(flagVerbose)
	    cout << "Checking seed " << getNetwork().getVertexName(*seedIt) << endl;
	mapPredecessors = getNetwork().getAllShortestPaths(*seedIt); 
	// Record nodes (I) involved in Ps,s' (s, s' € seeds) where P denotes shortest path
	for(preIt = mapPredecessors.begin(), preItEnd = mapPredecessors.end(); preIt != preItEnd; ++preIt)
	    if(setSeed.find(preIt->first) != seedItEnd)
		for(unsigned int i=0; i < (preIt->second).size(); ++i)
		    setIncluded.insert((preIt->second)[i]);
	
	//! check if included contains all and all needed -> seems to insert everything on the path without checking s'

	// For each i € I, check how many times i is involved in all possible Ps,t (s € seeds, t € nodes) 
	//! Count only once i for all possible distinct Ps,t for a given pair s,t
	//for(setIt = setIncluded.begin(), setItEnd = setIncluded.end(); setIt != setItEnd; ++setIt) {
	for(tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) {
	    //for(preIt = mapPredecessors.begin(), preItEnd = mapPredecessors.end(); preIt != preItEnd; ++preIt) {
		//v = preIt->first;
		v = *it;
		if(setIncluded.find(v) != setIncluded.end())
		    continue;
		if(flagVerbose)
		    cout << "- Checking shortest path to " << getNetwork().getVertexName(v) << endl;
		v_prev = mapPredecessors[v][0];
		while(v_prev != v) { 
		    for(unsigned int i=0; i < mapPredecessors[v].size(); ++i) {
			if(setIncluded.find(v) != setIncluded.end()) {
			    if(flagVerbose)
				cout << "-- Counting " << getNetwork().getVertexName(v) << endl;
			    if(setSeed.find(*it) != seedItEnd) {
				mapVertexToLocalCount[v] += 1;
			    }
			    mapVertexToGlobalCount[v] += 1;
			}
			//v_prev = mapPredecessors[v][i];
		    }
		    v_prev = v;
		    v = mapPredecessors[v][0];
		}
	    //} 
	}
	setIncluded.clear();
    }

	    /*
    Vertex ut, ut_prev;
    for(tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	//getNetwork().calculateShortestPath(*it, mapDistance, mapPredecessor); 
	mapPredecessors = getNetwork().getAllShortestPaths(*it); 
	//cout << score << endl;
	cout << getNetwork().getVertexName(*it) << endl;
	//for(vt = mapDistance.begin(), vtEnd = mapDistance.end(); vt != vtEnd; ++vt) 
	for(vt = mapPredecessors.begin(), vtEnd = mapPredecessors.end(); vt != vtEnd; ++vt) 
	{
	    //cout << getNetwork().getVertexName(vt->first) << " " << getNetwork().getVertexName(vt->second) << endl;
	    cout << getNetwork().getVertexName(vt->first) << ": ";
	    for(unsigned int i=0; i<((vector<Vertex>)(vt->second)).size(); ++i)
	       cout << getNetwork().getVertexName(vt->second[i]) << " ";
	    cout << endl;
	    //cout << vt->second << " ";
	    if(setSeed.find(vt->first) != setSeed.end()) {
		ut = vt->first;
		ut_prev = mapPredecessor[ut];
		while(ut_prev != ut) { 
		    mapVertexToLocalCount[ut_prev] += 1;
		    ut = ut_prev;
		    ut_prev = mapPredecessor[ut];
		}
	    } else {
		ut = vt->first;
		ut_prev = mapPredecessor[ut];
		while(ut_prev != ut) { 
		    mapVertexToGlobalCount[ut_prev] += 1;
		    ut = ut_prev;
		    ut_prev = mapPredecessor[ut];
		}
	    }
	}
    }
	    */

    float score = 0;
    for(tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	cout << getNetwork().getVertexName(*it)  << ": " << mapVertexToLocalCount[*it] << " " << mapVertexToGlobalCount[*it] << endl;
	if(mapVertexToLocalCount[*it] == 0) 
	    score = 0;
	else
	    score = float(mapVertexToLocalCount[*it]) / mapVertexToGlobalCount[*it]; 
	if(flagAccumulateToInitialNodeScore) 
	{
	    score += getNetwork().getVertexScore(*it);
	}
	getNetwork().setVertexScore(*it, score);
    }
    getNetwork().scaleVertexScores(SCALE_BETWEEN_ZERO_AND_ONE);
    getNetwork().outputScores(outputFile);
    return;
}



