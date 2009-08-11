
#include "Netzcore.hpp"

#include <iostream>

Netzcore::Netzcore() : ScoreNetwork()
{
}

Netzcore::Netzcore(std::string fileNode, std::string fileEdge, std::string fileOutput, std::string prefixSampled, unsigned int nSampled, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore, bool fResetSeedScoresToInitial, bool fVerbose) : ScoreNetwork(fileNode, fileEdge, fileOutput, fUseEdgeScore, fAccumulateToInitialNodeScore, fResetSeedScoresToInitial, fVerbose)
{ 
    prefixSampledGraphs = prefixSampled;
    nodeFileName = fileNode;
    nSampledGraphs = nSampled;
}

Netzcore::~Netzcore()
{
}

void Netzcore::loadSampledGraphs()
{
    std::ostringstream oss;
    for(unsigned int i = 1; i <= nSampledGraphs; i++)
    {	
	oss.str("");
	oss << prefixSampledGraphs << i;
	//std::cout << oss.str() << std::endl;
	Graph* pG = new Graph();
	pG->loadNodes(nodeFileName);
	pG->loadEdges(oss.str());
	sampledGraphs.push_back(pG);
    }
}

void Netzcore::printSampledGraphs()
{
    std::list<Graph*>::iterator it, itEnd;
    Graph* pG;
    for(it=sampledGraphs.begin(), itEnd=sampledGraphs.end(); it != itEnd; ++it)
    {	
	pG = *it;
	pG->print(true);
    }
}

void Netzcore::updateSampledGraphScores()
{
    std::list<Graph*>::iterator it, itEnd;
    Graph* pG;
    VertexIterator it, itEnd;
    for(it=sampledGraphs.begin(), itEnd=sampledGraphs.end(); it != itEnd; ++it)
    {	
	pG = *it;
	for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) {
	    pG->setVertexScore(pG->getVertex(getNetwork().getVertexName(*it)), getNetwork().getVertexScore(*it));
	}
    }
}

void Netzcore::initializeScoring() {
    VertexIterator it, itEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	//createVertexMessageMap(*it);
	//setVertexScoreInitial(*it, getNetwork().getVertexScore(*it));
	setVertexScoreUpdated(*it, 0.0);
    }
    loadSampledGraphs();
    //printSampledGraphs();
    //float vertex_scores[] = { 0.0, 1.0, 0.0, 1.0 };
    //std::pair<float, float> tempPair = calculateMeanAndSigma(vertex_scores, vertex_scores+sizeof(vertex_scores)/sizeof(float));
    //std::cout << tempPair.first << " " << tempPair.second << std::endl;
}

void Netzcore::finalizeIteration()
{
    scaleNodeScores(true);
    updateSampledGraphScores();
}

void Netzcore::updateNodeScore(Vertex v) 
{
    Vertex u;
    AdjVertexIterator vt, vtEnd;
    std::list<Graph*>::iterator it, itEnd;
    float tempScore = 0.0;
    std::pair<float, float> tempPair(0.0,0.0);
    unsigned int i = 0;
    Graph* pG;
    if(flagVerbose)
        std::cout << "-Checking node: " << getNetwork().getVertexName(v) << std::endl;

    std::list<float> scores;
    for(it=sampledGraphs.begin(), itEnd=sampledGraphs.end(); it != itEnd; ++it)
    {	
	pG = *it;
	u = pG->getVertex(getNetwork().getVertexName(v));
	//pG->print(false);
	for(boost::tie(vt, vtEnd) = pG->getAdjacentVertexIteratorOfVertex(u); vt != vtEnd; ++vt) 
	{
	    tempScore = pG->getVertexScore(*vt);
            if(flagUseEdgeScore) 
	    {
                tempScore *= pG->getEdgeScore(u, *vt);
            }
	    scores.push_back(tempScore);
	    if(flagVerbose)
		std::cout << "--Score of neighbor " << pG->getVertexName(*vt) << ": " << tempScore << std::endl; 
	}
    }
    tempScore = 0.0;
    i = 0;
    for(boost::tie(vt, vtEnd) = getNetwork().getAdjacentVertexIteratorOfVertex(v); vt != vtEnd; ++vt) 
    {
	if(flagUseEdgeScore) 
	{
	    tempScore += getNetwork().getEdgeScore(v, *vt) * getNetwork().getVertexScore(*vt);
	} 
	else 
	{
	    tempScore += getNetwork().getVertexScore(*vt);
	}
	i += 1;
    }
    tempScore /= i;
    tempPair = calculateMeanAndSigma(scores.begin(), scores.end());
    //std::cout << "s: " << tempScore << " m: " << tempPair.first << " sig: " << tempPair.second << std::endl;
    if(tempPair.second == 0) //(isnan(tempScore))
    {
	std::cout << "Zero variance!" << std::endl;	
	tempScore = 0.0;
    } 
    else 
    {
	tempScore = tempScore-tempPair.first;
	if(tempScore != 0)
	    tempScore /= tempPair.second;
    }
    // Update scoreUpdated (for error calculation)
    setVertexScoreUpdated(v, tempScore);
    if(flagVerbose)
	std::cout << "Score calculated for " << getNetwork().getVertexName(v) << ": " << tempScore << " after iteration " << iterationCounter << std::endl;
}

template <class InputIterator>
std::pair<float, float> calculateMeanAndSigma(InputIterator first, InputIterator beyond)
{
    float sum = 0.0, squareSum = 0.0, mean=0.0;
    unsigned int count = 0;
    while(first != beyond)
    {
	//std::cout << *first << std::endl;
	sum += *first;
	squareSum += sq(*first);
	count++;
	++first;
    }
    mean = sum/count;
    return std::make_pair(mean, squareSum/count-sq(mean));
}


