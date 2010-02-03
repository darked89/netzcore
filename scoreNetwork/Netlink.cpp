#include "Netlink.hpp"

#include <iostream>
#include <cstdlib>

Netlink::Netlink() : ScoreNetwork()
{
    //flagAccumulateToInitialNodeScore = false;
    //flagVerbose = false;
    threshold = 0.0;
}

Netlink::Netlink(std::string fileNode, std::string fileEdge, std::string fileOutput, float t, bool fAccumulateToInitialNodeScore, bool fVerbose) : ScoreNetwork(fileNode, fileEdge, fileOutput, false, fAccumulateToInitialNodeScore, false, fVerbose)
{
    //flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    //flagVerbose = fVerbose;
    //getNetwork().loadNodes(fileNode); 
    //getNetwork().loadEdges(fileEdge); 
    //outputFile = fileOutput;
    threshold = t;
}

Netlink::~Netlink() {
}

void Netlink::updateNodeScore(Vertex v)
{ 
    AdjVertexIterator ut, utEnd;
    float score = 0.0;
    unsigned int i = 0;
    if(flagVerbose) {
	std::cout << "Checking node " << getNetwork().getVertexName(v) << std::endl;
    }
    for(boost::tie(ut, utEnd) = getNetwork().getAdjacentVertexIteratorOfVertex(v); ut != utEnd; ++ut) {
	if(v != *ut) { // Skip self edges
	    if(flagVerbose) {
		std::cout << "- Score from neighbor " << getNetwork().getVertexName(*ut) << ": " << getNetwork().getVertexScore(*ut) << std::endl;
	    }
	    score += getNetwork().getVertexScore(*ut);
	    i += 1;
	}
    }
    // Scale by number of neighbors
    //if(i!=0) score /= i;
    if(score >= threshold) {
	score = 1.0;
    } else {
	score = 0.0;
    }
    if(flagAccumulateToInitialNodeScore) {
	score += getNetwork().getVertexScore(v);
    }
    setVertexScoreUpdated(v, score);
    return;
}

/*
void Netlink::run(float threshold)
{ 
    VertexIterator it, itEnd;
    AdjVertexIterator ut, utEnd;
    float score = 0.0;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	if(flagVerbose) {
	    std::cout << "Checking node " << getNetwork().getVertexName(*it) << std::endl;
	}
	score = 0.0;
	for(boost::tie(ut, utEnd) = getNetwork().getAdjacentVertexIteratorOfVertex(*it); ut != utEnd; ++ut) {
	    if(*it != *ut) { // Skip self edges
		if(flagVerbose) {
		    std::cout << "- Score from neighbor " << getNetwork().getVertexName(*ut) << ": " << getNetwork().getVertexScore(*ut) << std::endl;
		}
		score += getNetwork().getVertexScore(*ut);
	    }
	}
	if(score >= threshold) {
	    score = 1.0;
	}
	if(flagAccumulateToInitialNodeScore) {
	    score += getNetwork().getVertexScore(*it);
	}
	getNetwork().setVertexScore(*it, score);
    }
    getNetwork().outputScores(outputFile);
    return;
}
*/

