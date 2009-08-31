#include "Netrank.hpp"

#include <iostream>

Netrank::Netrank() 
{
    flagAccumulateToInitialNodeScore = false;
}

Netrank::Netrank(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fAccumulateToInitialNodeScore)
{
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    network.loadNodes(fileNode); 
    network.loadEdges(fileEdge); 
    outputFile = fileOutput;
}

Netrank::~Netrank() {
}

void Netrank::run(unsigned int nIteration) 
{ 
    std::map<Vertex, float>::iterator vt, vtEnd;
    std::map<Vertex, float> mapRank;
    getNetwork().calculatePageRank(mapRank, nIteration);
    float score = 0.0;
    for(vt = mapRank.begin(), vtEnd = mapRank.end(); vt != vtEnd; ++vt) 
    {
	score = 0.0;
	if(flagAccumulateToInitialNodeScore) {
	    score = getNetwork().getVertexScore(vt->first);
	}
	score += vt->second;
	getNetwork().setVertexScore(vt->first, score);
    }
    getNetwork().scaleVertexScores(SCALE_BETWEEN_ZERO_AND_ONE);
    getNetwork().outputScores(outputFile);
    return;
}


