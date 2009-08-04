#ifndef NETSHORT_HPP_9054862562754458964905368042156
#define NETSHORT_HPP_9054862562754458964905368042156

#include "graph/Graph.hpp"

class Netshort {
private:
    // MEMBERS
    Graph network;
    
    std::string outputFile;
    bool flagAccumulateToInitialNodeScore; // if true, accumulated shortest path score is added to node score
    //void accumulateAllShortestPathScores(); 

public:
    Netshort(); 
    Netshort(std::string fileNode, std::string fileEdge, std::string fileOutput, bool flagAccumulateToInitialNodeScore = false);
    ~Netshort();
    Graph & getNetwork() { return network; };
    void run(); 
    float testRun(); 
};

#endif // NETSHORT_HPP_

