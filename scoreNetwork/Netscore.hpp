#ifndef NETSCORE_HPP_9054862562754458964905368042156
#define NETSCORE_HPP_9054862562754458964905368042156

#include "ScoreNetwork.hpp"

#include <boost/unordered_map.hpp>

// multiplication of weight of edges on the path, iteration of arrival, number of equidistant arrivals
typedef boost::tuple<float, int, int> Message;
typedef boost::unordered_map<unsigned int,  Message> UIntToMessage;
//typedef std::map<unsigned int,  Message> UIntToMessage;
typedef boost::unordered_map<Vertex, UIntToMessage* > VertexToMessageMap;
//typedef std::map<Vertex, UIntToMessage* > VertexToMessageMap;

class Netscore: public ScoreNetwork {
public:
    Netscore(); 
    Netscore(std::string fileOutput, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    Netscore(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    ~Netscore();
    //setNetwork(const Graph * g);

    // methods below were protected before, made public to be able to use netzcore & netscore in combination
    void initializeScoring();
    void finalizeScoring();
    void initializeRepetition();
    void finalizeIteration();
    void updateNodeScore(Vertex v);

private:
    UIntToMessage * getVertexMessageMap(Vertex const v) { return messageMaps[v]; };
    void createVertexMessageMap(Vertex const v) { messageMaps[v] = new UIntToMessage(); };
    float getVertexScoreInitial(Vertex const v) const { return mapScoreInitial.at(v); };
    void setVertexScoreInitial(Vertex const v, float vData) { mapScoreInitial[v] = vData; };

    // MEMBERS
    VertexToMessageMap messageMaps;
    VertexToFloat mapScoreInitial;
};

#endif // NETSCORE_HPP_

