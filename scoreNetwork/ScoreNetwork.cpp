
#include "ScoreNetwork.hpp"
#include <iostream>
//#include <cstdlib> 
//#include <ctime>
//#include <sstream>

using namespace std;

ScoreNetwork::ScoreNetwork() 
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

ScoreNetwork::ScoreNetwork(string fileNode, string fileEdge, string fileOutput, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore, bool fResetSeedScoresToInitial, bool fVerbose) 
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
    //cout << "Snet file:" << fileOutput << " " << outputFile << endl;
    //cout << "Snet fverbose:" << fVerbose << " " << flagVerbose << endl;
    //cout << "Snet fAccumulate: " << this->flagAccumulateToInitialNodeScore << endl;
}

ScoreNetwork::~ScoreNetwork() {
}

void ScoreNetwork::run(int nRepeatition, int nIteration, float tError) 
{ 
    nIteration = nIteration;
    float error = INFINITY;
    initializeScoring();
    //cout << "Snet fAccumulate in run: " << flagAccumulateToInitialNodeScore << endl;
    for(int repeatCounter = 1; repeatCounter<=nRepeatition and error > tError; ++repeatCounter) {
	initializeRepeatition();
	for(iterationCounter = 1; iterationCounter<=nIteration and error > tError; ++iterationCounter) {
	    initializeIteration();
	    updateNetwork();
	    error = calculateErrorAndUpdateScores();
	    finalizeIteration();
	} 
	finalizeRepeatition();
    }
    finalizeScoring();
    //cout << "Snet run file: " << this->outputFile << endl;
    getNetwork().outputScores(outputFile);
}

/*
// check which is called
void ScoreNetwork::initializeIteration() {
    cout << endl << "initialize of base" << endl;  
}

//virtual void ScoreNetwork::finalizeIteration() { }

// should it be in destructor
void ScoreNetwork::finalizeScoring() {
    scaleNodeScores();
}
*/
/*
void ScoreNetwork::updateNodeScore(Vertex v) {}
void ScoreNetwork::updateEdgeScore(Edge e) {}
*/

void ScoreNetwork::updateNetwork() 
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

float ScoreNetwork::calculateErrorAndUpdateScores() {
    float error = 0.0, sumErrorNode = 0.0; //, sumErrorEdge = 0.0;
    VertexIterator it, itEnd;
    OutEdgeIterator et, etEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
        if(isnan(getNetwork().getVertexScore(*it)) or isinf(getNetwork().getVertexScore(*it))) {
	   cerr << "NAN or INF score" << endl;
	}
        //! Reseting to initial scores? 
	//if(flagResetSeedScoresToInitial) {
            //if(setSource.find(pVertex->id) != setSource.end()) { // for fFlow compatibility
                //!pVertex->data.scoreUpdated = pVertex->data.scoreInitial; 
            //}
        //}
        error = getVertexScoreUpdated(*it) - getNetwork().getVertexScore(*it); 
        sumErrorNode += sq(error);
	getNetwork().setVertexScore(*it, getVertexScoreUpdated(*it)); 
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

float ScoreNetwork::transferScore(float score, float a, float b) {
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

void ScoreNetwork::scaleNodeScores() {
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


