#ifndef NETZCORE_HPP_9054862562754458964905368042156
#define NETZCORE_HPP_9054862562754458964905368042156

#include "ScoreNetwork.hpp"

//#include <boost/unordered_map.hpp>

#include <list>

//typedef boost::unordered_map<unsigned int,  Message> UIntToMessage;
///typedef boost::unordered_map<Vertex, UIntToMessage* > VertexToMessageMap;

template <class InputIterator>
std::pair<float, float> calculateMeanAndSigma(InputIterator first, InputIterator beyond);

class Netzcore: public ScoreNetwork {
public:
    Netzcore(); 
    Netzcore(std::string fileNode, std::string fileOutput, std::string prefixSampled, unsigned int nSample, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    Netzcore(std::string fileNode, std::string fileEdge, std::string fileOutput, std::string prefixSampled, unsigned int nSample, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    ~Netzcore();

    // methods below were protected before, made public to be able to use netzcore & netscore in combination
    void initializeScoring();
    void finalizeScoring();
    void initializeIteration();
    void finalizeIteration();
    void updateNodeScore(Vertex v);

private:
    //UIntToMessage * getVertexMessageMap(Vertex const v) { return messageMaps[v]; };
    //void createVertexMessageMap(Vertex const v) { messageMaps[v] = new UIntToMessage(); };

    //float getVertexScoreInitial(Vertex const v) const { return mapScoreInitial.at(v); };
    //void setVertexScoreInitial(Vertex const v, float vData) { mapScoreInitial[v] = vData; };

    void loadSampledGraphs();
    void printSampledGraphs();
    void updateSampledGraphScores();

    // MEMBERS
    //VertexToMessageMap messageMaps;
    //VertexToFloat mapScoreInitial;
    std::list<Graph*> sampledGraphs;
    std::string nodeFileName;
    std::string prefixSampledGraphs;
    unsigned int nSampledGraphs;
    float mean;
    float sigma;
};

#endif // NETZCORE_HPP_

