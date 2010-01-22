#ifndef NETLINK_HPP_9054862562754458964905368042156
#define NETLINK_HPP_9054862562754458964905368042156

#include "graph/Graph.hpp"
#include "ScoreNetwork.hpp"

class Netlink: public ScoreNetwork {
private:
    // MEMBERS
    float threshold;
    // inherited from ScoreNetwork
    //Graph network;
    //std::string outputFile;
    //bool flagAccumulateToInitialNodeScore; // if true, linker degree score is added to node score
    //bool flagVerbose;

public:
    Netlink(); 
    Netlink(std::string fileNode, std::string fileEdge, std::string fileOutput, float t, bool flagAccumulateToInitialNodeScore = false, bool flagVerbose = false);
    ~Netlink();
    //Graph & getNetwork() { return network; };
    //void run(int threshold); 
    void updateNodeScore(Vertex v);
};

#endif // NETLINK_HPP_

