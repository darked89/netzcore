#ifndef NETRANDOM_HPP_9054862562754458964905368042156
#define NETRANDOM_HPP_9054862562754458964905368042156

#include "graph/Graph.hpp"

class Netrandom {
private:
    // MEMBERS
    Graph network;
    
    std::string outputFile;
    bool flagAccumulateToInitialNodeScore; // if true, accumulated shortest path score is added to node score

public:
    Netrandom(); 
    Netrandom(std::string fileNode, std::string fileEdge, std::string fileOutput, bool flagAccumulateToInitialNodeScore = false);
    ~Netrandom();
    Graph & getNetwork() { return network; };
    void run(); 
};

#endif // NETRANDOM_HPP_

