#ifndef NETSCORE2_HPP_9054862562754458964905368042156
#define NETSCORE2_HPP_9054862562754458964905368042156

#include "graph/Graph.hpp"

#include "Exceptions.hpp"

//#define CONTROL 0
#define sq(x) (x*x)


class Netscore2 {
public:
    Netscore2(); 
    Netscore2(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeReliabilityScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    ~Netscore2();
    void run(int nRepetition, int nIteration, float tError = PRECISION); 
    Graph & getNetwork() { return network; };

private:
    // MEMBERS
    Graph network;
    int nIteration;
    int iterationCounter;
    float minScore;
    float maxScore;
    bool flagUseEdgeScore;
    bool flagAccumulateToInitialNodeScore;
    bool flagResetSeedScoresToInitial;
    bool flagVerbose;
    std::string outputFile;

    void updateEdgeScore(Edge e) {}; 
    void finalizeScoring() {};
    void finalizeRepetition() {};
    void initializeIteration() {};
    void finalizeIteration() {}; 

    float calculateErrorAndUpdateScores();
    void updateNetwork();
    void scaleNodeScores();
    float transferScore(float score, float a = 0.5, float b = 0.1); 

    static const float PRECISION = 0.00001;
    enum TransferType { IDENTITY, POLYNOMIAL, LOGARITHMIC, EXPONENTIAL };
    int typeTransfer;
    
 
    void initializeScoring();
    void initializeRepetition();
    void updateNodeScore(Vertex v);
};

#endif // NETSCORE2_HPP_

