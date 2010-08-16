#include "Netshort.hpp"

#include <iostream>

Netshort::Netshort() 
{
    flagAccumulateToInitialNodeScore = false;
}

Netshort::Netshort(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fAccumulateToInitialNodeScore)
{
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    network.loadNodes(fileNode); 
    network.loadEdges(fileEdge, true); // flag causes edge scores to be inverted (1/weight)
    outputFile = fileOutput;
}

Netshort::~Netshort() {
}

void Netshort::run() 
{ 
    VertexIterator it, itEnd;
    std::map<Vertex, float>::iterator vt, vtEnd;
    std::map<Vertex, Vertex> mapPredecessor;
    std::map<Vertex, float> mapDistance;
    float score = 0;
    //std::cout << flagAccumulateToInitialNodeScore << std::endl;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	getNetwork().calculateShortestPath(*it, mapDistance, mapPredecessor); 
	score = 0;
	//std::cout << score << std::endl;
	for(vt = mapDistance.begin(), vtEnd = mapDistance.end(); vt != vtEnd; ++vt) 
	{
	    //std::cout << vt->second << " ";
	    //score += vt->second;
	    if(vt->second != 0)
		score += 1/vt->second;
	}
	//getNetwork().setVertexScoreUpdated(*it, score);
	//score = 1 / score; // # 1000/score #! depends on the initial weights
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



/////
float Netshort::testRun() {
    std::map<Vertex, Vertex> mapPredecessor;
    std::map<Vertex, float> mapDistance;
    std::map<Vertex, float>::iterator vt, vtEnd;
    getNetwork().calculateShortestPath(getNetwork().getVertex("25604"), mapDistance, mapPredecessor); 
    //mapDistance = getNetwork().calculateShortestPath(getNetwork().getVertex("25604")); 
    float score = 0;
    for(vt = mapDistance.begin(), vtEnd = mapDistance.end(); vt != vtEnd; ++vt) 
    {
	//std::cout << getNetwork().getVertexIndex(vt->first) << std::endl;
	//std::cout << getNetwork().getVertexName(vt->first) << " " << vt->second << std::endl;
	score += vt->second;
    }
    return 1000/score;
}


