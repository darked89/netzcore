
#include "Netscore2.hpp"

#include <iostream>

Netscore2::Netscore2()
{
    typeTransfer = IDENTITY;
    minScore = INFINITY;
    maxScore = -INFINITY;
    nIteration = 0;
    iterationCounter = 0;
    flagVerbose = false;
    flagUseEdgeScore = true;
    flagAccumulateToInitialNodeScore = true;
    flagResetSeedScoresToInitial = false;

}

Netscore2::Netscore2(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore, bool fResetSeedScoresToInitial, bool fVerbose)
{
    typeTransfer = IDENTITY;
    minScore = INFINITY;
    maxScore = -INFINITY;
    nIteration = 0;
    iterationCounter = 0;
    flagVerbose = fVerbose;
    flagUseEdgeScore = fUseEdgeScore; // should? be in updateEdgeScore
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore; // should? be in updateNodeScore
    flagResetSeedScoresToInitial = fResetSeedScoresToInitial; 
    outputFile = fileOutput;
    network.loadNodes(fileNode); 
    network.loadEdges(fileEdge); 
}

Netscore2::~Netscore2()
{
}

void Netscore2::run(int nRepetition, int nIteration, float tError) 
{ 
    nIteration = nIteration;
    float error = INFINITY;
    initializeScoring();
    //cout << "Snet fAccumulate in run: " << flagAccumulateToInitialNodeScore << endl;
    for(int repeatCounter = 1; repeatCounter<=nRepetition and error > tError; ++repeatCounter) {
	initializeRepetition();
	for(iterationCounter = 1; iterationCounter<=nIteration and error > tError; ++iterationCounter) {
	    initializeIteration();
	    updateNetwork();
	    error = calculateErrorAndUpdateScores();
	    finalizeIteration();
	} 
	finalizeRepetition();
    }
    finalizeScoring();
    //cout << "Snet run file: " << this->outputFile << endl;
    getNetwork().outputScores(outputFile);
}


void Netscore2::initializeScoring() {
    VertexIterator it, itEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	getNetwork().createVertexMessageMap(*it);
    }
}

void Netscore2::initializeRepetition()
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

void Netscore2::updateNodeScore(Vertex v) 
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

void Netscore2::updateNetwork() 
{ 
    VertexIterator it, itEnd;
    //OutEdgeIterator et, etEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
        updateNodeScore(*it); 
	/*
	for(boost::tie(et, etEnd) = getNetwork().getEdgeIteratorOfVertex(*it); et != etEnd; ++et) 
	{
	    // Make sure that we update edges only once since it is a an undirected graph
	    if(getNetwork().getVertexIndex(*it) >= getNetwork().getVertexIndex(getNetwork().getTargetVertex(*et))) {
		updateEdgeScore(*et);
	    }
	}
	*/
    }
    //if(flagVerbose)	
    //        printNetwork("After updateNetwork");
}

float Netscore2::calculateErrorAndUpdateScores() {
    float error = 0.0, sumErrorNode = 0.0; //, sumErrorEdge = 0.0;
    VertexIterator it, itEnd;
    OutEdgeIterator et, etEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
        if(isnan(getNetwork().getVertexScore(*it)) or isinf(getNetwork().getVertexScore(*it))) {
	    std::cerr << "NAN or INF score" << std::endl;
	}
        //! Reseting to initial scores? 
	//if(flagResetSeedScoresToInitial) {
            //if(setSource.find(pVertex->id) != setSource.end()) { // for fFlow compatibility
                //!pVertex->data.scoreUpdated = pVertex->data.scoreInitial; 
            //}
        //}
        error = getNetwork().getVertexScoreUpdated(*it) - getNetwork().getVertexScore(*it); 
        sumErrorNode += sq(error);
	getNetwork().setVertexScore(*it, getNetwork().getVertexScoreUpdated(*it)); 
        //!pVertex->data.score = pVertex->data.scoreUpdated;
	// For the time being, not considering updation of edge scores
	/*
	for(boost::tie(et, etEnd) = getNetwork().getEdgeIteratorOfVertex(*it); et != etEnd; ++et) 
	{
	    // Make sure that we update edges only once since it is a an undirected graph
	    if(getNetwork().getVertexIndex(*it) >= getNetwork().getVertexIndex(getNetwork().getTargetVertex(*et))) {
                error = getNetwork().getEdgeScore(*et);
		sumErrorEdge += sq(error);
		//!pEdge->data.score = pEdge->data.scoreUpdated;	            
	    }
	}
	*/
    }
    //cerr << "minScore: " << minScore << " maxScore: " << maxScore << endl;
    //return max(float(sqrt( (sumErrorEdge/nEdge)*(sumErrorNode/nNode) )), float( (sqrt(sumErrorNode/nNode)+ sqrt(sumErrorEdge/nEdge)) /2)); //float(sqrt(sumErrorEdge/nEdge)); //float(sqrt(sumErrorNode/nNode)+sqrt(sumErrorEdge/nEdge))/2;
    return float(sqrt( sumErrorNode/getNetwork().getSize() ));
}

float Netscore2::transferScore(float score, float a, float b) {
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
        //cout << "Warning: unidentified type" << endl;
	throw new TypeException(); //typeEx;
        return 0.0;
    }
}

void Netscore2::scaleNodeScores() {
    VertexIterator it, itEnd;
    float value = 0;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	value = getNetwork().getVertexScore(*it);
	maxScore = (value > maxScore)?value:maxScore;
	minScore = (value < minScore)?value:minScore;
    }
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	getNetwork().setVertexScore(*it, getNetwork().getVertexScore(*it)/maxScore);
    }
    return;
}


