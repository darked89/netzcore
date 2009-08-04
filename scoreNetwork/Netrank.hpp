#ifndef NETRANK_HPP_9054862562754458964905368042156
#define NETRANK_HPP_9054862562754458964905368042156

#include "graph/Graph.hpp"

class Netrank {
private:
    // MEMBERS
    Graph network;
    
    std::string outputFile;
    bool flagAccumulateToInitialNodeScore; // if true, accumulated shortest path score is added to node score

public:
    Netrank(); 
    Netrank(std::string fileNode, std::string fileEdge, std::string fileOutput, bool flagAccumulateToInitialNodeScore = false);
    ~Netrank();
    Graph & getNetwork() { return network; };
    void run(unsigned int nIteration = 20); 
};

#endif // NETRANK_HPP_

