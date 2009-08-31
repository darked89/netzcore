#include "Netrandom.hpp"

#include <iostream>

#include <ctime>
#include <cstdlib>

Netrandom::Netrandom() 
{
    flagAccumulateToInitialNodeScore = false;
}

Netrandom::Netrandom(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fAccumulateToInitialNodeScore)
{
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    network.loadNodes(fileNode); 
    network.loadEdges(fileEdge); 
    outputFile = fileOutput;
}

Netrandom::~Netrandom() {
}

void Netrandom::run() 
{ 
    VertexIterator it, itEnd;
    srand ( time(NULL) );
    float score = 0.0;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	score = 0.0;
	if(flagAccumulateToInitialNodeScore) {
	    score = getNetwork().getVertexScore(*it);
	}
	score += float(rand())/RAND_MAX;
	getNetwork().setVertexScore(*it, score);
    }
    getNetwork().outputScores(outputFile);
    return;
}


