#ifndef SCORENETWORK_HPP_29045684568920568920456590268
#define SCORENETWORK_HPP_29045684568920568920456590268

#include <cmath> // for INFINITY

#include "graph/Graph.hpp"

#include "Exceptions.hpp"

//#define CONTROL 0
#define sq(x) (x*x)

//typedef boost::unordered_map<Vertex, float> VertexToFloat;
//typedef std::map<Vertex, float> VertexToFloat;
typedef boost::unordered_map<Vertex, float> VertexToFloat;

class ScoreNetwork {
public:
    ScoreNetwork(); 
    ScoreNetwork(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false);
    virtual ~ScoreNetwork();
    Graph & getNetwork() { return network; };
    void run(int nRepetition, int nIteration, float tError = PRECISION); 
    void scaleNodeScoresWrapper() { scaleNodeScores(); };
    /*
    //void outputGraph(Graph const & g, std::string fileName);
    void outputGraph(Graph *g, std::string fileName);
    void printNetwork(std::string explanation = "");
    void printGraph(Graph const & g, std::string explanation = "");
    */

protected:
    int nIteration;
    int iterationCounter;
    float minScore;
    float maxScore;
    bool flagUseEdgeScore;
    bool flagAccumulateToInitialNodeScore;
    bool flagResetSeedScoresToInitial;
    bool flagVerbose;
    std::string outputFile;

    //void filterNetwork(int degreeNodeMin = -INT_INFINITY, int degreeNodeMax = INT_INFINITY, float scoreNodeMin = -INFINITY, float scoreNodeMax = INFINITY, float scoreEdgeMin = -INFINITY, float scoreEdgeMax = INFINITY, bool flagRemoveSelfEdges = true);
    //void scaleScores(float scoreAllowedMaxNode = INFINITY, float scoreAllowedEdgeMax = INFINITY, bool flagUseZScoring = false);
    //void resetScores();

    virtual void updateNodeScore(Vertex v) {}; // = 0; // not making it abstract (interface) for testing purposes
    virtual void updateEdgeScore(Edge e) {}; 
    virtual void initializeScoring() {}; //= 0;
    virtual void finalizeScoring() {};
    virtual void initializeRepeatition() {}; // = 0; 
    virtual void finalizeRepeatition() {};
    virtual void initializeIteration() {};
    virtual void finalizeIteration() {}; 

    float calculateErrorAndUpdateScores();
    void updateNetwork();
    void scaleNodeScores();
    float transferScore(float score, float a = 0.5, float b = 0.1); 

    float getVertexScoreUpdated(Vertex const v) const { return mapScoreUpdated.at(v); };
    void setVertexScoreUpdated(Vertex const v, float vData) { mapScoreUpdated[v] = vData; };

private:
    // CONSTANTS (MACROS)
    //static const float MIN_SCORE = 0.0;
    //static const float MAX_SCORE = 1.0;
    static const float PRECISION = 0.00001;
    //static const enum TransferType { IDENTITY, POLYNOMIAL, LOGARITHMIC, EXPONENTIAL };
    enum TransferType { IDENTITY, POLYNOMIAL, LOGARITHMIC, EXPONENTIAL };
    int typeTransfer;
    
    // MEMBERS
    Graph network;
    //unordered_set<string> setSource;
    VertexToFloat mapScoreUpdated;
};

#endif // SCORENETWORK_HPP_

