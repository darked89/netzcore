
#include "Netscore.hpp"

#include <iostream>

Netscore::Netscore() : ScoreNetwork()
{
}

Netscore::Netscore(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore, bool fResetSeedScoresToInitial, bool fVerbose) : ScoreNetwork(fileNode, fileEdge, fileOutput, fUseEdgeScore, fAccumulateToInitialNodeScore, fResetSeedScoresToInitial, fVerbose)
{
    //std::cout << "Nscore fAccumulate: " << flagAccumulateToInitialNodeScore << std::endl;
    //std::cout << "Nscore file: " << fileOutput << std::endl;
    //if base constructor is not called, need to assign all like 
    //outputFile = fileOutput;
    //flagVerbose = flagVerbose;
    //this->flagAccumulateToInitialNodeScore = flagAccumulateToInitialNodeScore; 
    // Problem caused by pythonic assignment of instance variable with the same name as argument
}

Netscore::~Netscore()
{
}

void Netscore::initializeScoring() {
    VertexIterator it, itEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	getNetwork().createVertexMessageMap(*it);
    }
}

void Netscore::initializeRepeatition()
{
    VertexIterator it, itEnd;
    UIntToMessage * pMapMessage;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	//std::cout << getNetwork().getVertexName(*it) << std::endl;
	//mapMessage = getNetwork().getVertexMessageMap(*it); 
	//mapMessage[getNetwork().getVertexIndex(*it)] = std::pair<float, int>(1.0, 0);
	pMapMessage = getNetwork().getVertexMessageMap(*it); 
	(*pMapMessage)[getNetwork().getVertexIndex(*it)] = boost::make_tuple(1.0, 0, 1); //std::pair<float, int>(1.0, 0);
	//UIntToMessage::iterator mt, mtEnd;
	//for(mt=mapMessage.begin(), mtEnd=mapMessage.end(); mt != mtEnd; ++mt)
	//    std::cout << mt->first << " " << mt->second.second << std::endl;
	//mapMessage = getNetwork().getVertexMessageMap(*it); 
	//for(mt=mapMessage.begin(), mtEnd=mapMessage.end(); mt != mtEnd; ++mt)
	//    std::cout << mt->first << " " << mt->second.second << std::endl;
    }
}

void Netscore::updateNodeScore(Vertex v) 
{
    AdjVertexIterator vt, vtEnd;
    UIntToMessage::iterator it, itEnd, itSearch;
    UIntToMessage * pMapMessage = getNetwork().getVertexMessageMap(v), mapMessageNeighbor;
    float tempScore = 0.0;
    //boost::unordered_map<unsigned int, unsigned int> mapIterationToCount;
    unsigned int nMessagesRecieved = 0;
    //std::cout << "Nscore fverbose : " << flagVerbose << " fAccumulate: " << flagAccumulateToInitialNodeScore << std::endl; //<< " T F: " << true << " " << false <<std::endl;
    if(flagVerbose)
        std::cout << "-Checking node: " << getNetwork().getVertexName(v) << std::endl;
    for(boost::tie(vt, vtEnd) = getNetwork().getAdjacentVertexIteratorOfVertex(v); vt != vtEnd; ++vt) 
    {
	if(flagVerbose)
	    std::cout << "--Messages from neighbor: " << getNetwork().getVertexName(*vt) << std::endl; 
	mapMessageNeighbor = *(getNetwork().getVertexMessageMap(*vt)); 
	for(it=mapMessageNeighbor.begin(), itEnd=mapMessageNeighbor.end(); it != itEnd; ++it)
	{
	    //std::cout << boost::get<1>(it->second) << " " << iterationCounter << std::endl;
	    // Consider the message if only coming from a previous iteration
	    if(boost::get<1>(it->second) < iterationCounter)
	    {
		//tempScore = it->second.first;
		tempScore = boost::get<0>(it->second);
		if(flagUseEdgeScore) 
		{
		    tempScore *= getNetwork().getEdgeScore(v, *vt);
		}
		if(flagVerbose)
		    std::cout << "---message of " << it->first << " " << getNetwork().getVertexName(it->first) << " score: " << tempScore << std::endl; 
		// Sum up messages from the same node arriving at the same time
		itSearch = pMapMessage->find(it->first);
		if(itSearch == pMapMessage->end())
		{
		    (*pMapMessage)[it->first] = boost::make_tuple(tempScore, iterationCounter, 1); //std::pair<float, int>(tempScore, iterationCounter); // iterationCounter == it->second.second + 1
		} 
		else 
		{
		    if(boost::get<1>(itSearch->second) == iterationCounter) 
		    {
			boost::get<0>(itSearch->second) += tempScore;
			boost::get<2>(itSearch->second) += 1;
		    }
		}
	    }
	}
    }
    tempScore = 0.0;
    for(it=pMapMessage->begin(), itEnd=pMapMessage->end(); it != itEnd; ++it)
    {
	// Vertex score contains accumulated score over the iterations so only consider messages from that iteration
	if(boost::get<1>(it->second) == iterationCounter) 
	{
	    tempScore += boost::get<0>(it->second) * getNetwork().getVertexScoreInitial(getNetwork().getVertex(it->first));
	    //nMessagesRecieved += 1;
	    nMessagesRecieved += boost::get<2>(it->second);
	    if(flagVerbose)
		std::cout << "Evaluating message of " << it->first << " " << getNetwork().getVertexName(it->first) << " score: " << tempScore << "(+ " << boost::get<0>(it->second) << " * " << getNetwork().getVertexScoreInitial(getNetwork().getVertex(it->first)) <<" )" << std::endl; 
	}
    }
    if(nMessagesRecieved != 0) 
    {
        tempScore /= nMessagesRecieved;
    }
    tempScore += getNetwork().getVertexScore(v);
    // Remove initial score if accumulation to initial score is not desired
    if(flagAccumulateToInitialNodeScore == false && iterationCounter == 1) 
    {
    	tempScore -= getNetwork().getVertexScoreInitial(v);
    } 
    // Update scoreUpdated (for error calculation)
    getNetwork().setVertexScoreUpdated(v, tempScore);
    if(flagVerbose)
	std::cout << "Score calculated for " << getNetwork().getVertexName(v) << ": " << tempScore << " after iteration " << iterationCounter << std::endl;
}

