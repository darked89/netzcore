#include "Netween.hpp"

#include <iostream>

Netween::Netween() 
{
    flagAccumulateToInitialNodeScore = false;
}

Netween::Netween(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fAccumulateToInitialNodeScore)
{
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    outputFile = fileOutput;
    network.loadNodes(fileNode); 
    network.loadEdges(fileEdge);//, true); 
    VertexIterator it, itEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	setSeed.insert(*it);
    }
}

Netween::~Netween() {
}

void Netween::run() 
{ 
    VertexIterator it, itEnd;
    boost::unordered_map<Vertex, std::pair<int, int> > mapVertexToLocalAndGlobalCounts;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	mapVertexToLocalAndGlobalCounts[*it] = 0;
    }

    std::map<Vertex, float>::iterator vt, vtEnd;
    std::map<Vertex, Vertex> mapPredecessor;
    std::map<Vertex, float> mapDistance;
    float score = 0;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	getNetwork().calculateShortestPath(*it, mapDistance, mapPredecessor); 
	score = 0;
	//std::cout << score << std::endl;
	for(vt = mapDistance.begin(), vtEnd = mapDistance.end(); vt != vtEnd; ++vt) 
	{
	    //std::cout << getNetwork().getVertexName(vt->first) << " " << vt->second << std::endl;
	    //std::cout << vt->second << " ";
	    score += vt->second;
	}
	//getNetwork().setVertexScoreUpdated(*it, score);
	score = 1 / score; // # 1000/score #! depends on the initial weights
	if(flagAccumulateToInitialNodeScore) 
	{
	    score += getNetwork().getVertexScore(*it);
	}
	getNetwork().setVertexScore(*it, score);
	//std::cout << std::endl << getNetwork().getVertexIndex(*it)  << " " << getNetwork().getVertexName(*it) << " " << 10000/score << std::endl;
    }
    getNetwork().scaleVertexScores(SCALE_BETWEEN_ZERO_AND_ONE);
    getNetwork().outputScores(outputFile);
    return;
}



