#ifndef NETSCORE_HPP_9054862562754458964905368042156
#define NETSCORE_HPP_9054862562754458964905368042156

#include "ScoreNetwork.hpp"

class Netscore: public ScoreNetwork {
public:
    Netscore(); 
    Netscore(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeReliabilityScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    ~Netscore();

private:
    void initializeScoring();
    void initializeRepeatition();
    void updateNodeScore(Vertex v);

    // MEMBERS
    //VertexToMessageMap 
};

#endif // NETSCORE_HPP_

