#ifndef SCORENETWORK_H_
#define SCORENETWORK_H_

#include "graph/globals.h"
#include "graph/Graph.h"
#include <cmath> // for INFINITY

#define CONTROL 0

#define N_PERTURBATION_MIN 100 // 5
#define PRECISION_COEFFICIENT 1000000
#define PRECISION_INIT 0.000000f
#define MIN_SCORE_SHIFT (10.0/PRECISION_COEFFICIENT) //0.00001
#define N_SAMPLE 20 //75 //100
#define N_TRIAL_MAX 10
#define MAX_Z_SCORE 5.0

#define MIN_SCORE 0.0
#define MAX_SCORE INFINITY

#define IDENTITY 1
#define POLYNOMIAL 2
#define LOGARITHMIC 3
#define EXPONENTIAL 4

#define TOPOLOGY_PRESERVING 100
#define TOPOLOGY_PRESERVING_SAME_DEGREE_INTERCHANGE 101
#define RANDOM 102
#define RANDOM_PRESERVING_ONE_DEGREED_NODES 103
#define DEGREE_DISTRIBUTION_PRESERVING 104
#define DEGREE_DISTRIBUTION_PRESERVING_EDGE_INTERCHANGE 105

#define RANDOM_SAMPLING_DIRECTORY "randomNetworks"

//#define SIF_FORMAT 0
//#define COMPACT_FORMAT 1
//#define DETAILED_FORMAT 2

class ScoreNetwork {
private:
    Graph network;
    vector<Graph*> listNetworkSampled;
    hash_map<int, void*> mapSource;
    int nIteration;
    int iterationCounter;

    float minZScore;
    float maxZScore;
    float minScore;
    float maxScore;
    bool flagDebug;
    bool flagVerbose;
    int typeTransfer;
    bool flagUseEdgeReliabilityScore;
    bool flagAccumulate;
    bool flagResetSeedScoresToInitial;

    void loadProteins(string const &fileName);
    void loadInteractions(string const &fileName);    
    void filterNetwork(int degreeNodeMin = -INT_INFINITY, int degreeNodeMax = INT_INFINITY, float scoreNodeMin = -INFINITY, float scoreNodeMax = INFINITY, float scoreEdgeMin = -INFINITY, float scoreEdgeMax = INFINITY, bool flagRemoveSelfEdges = true);
    //void scaleScores(float scoreAllowedMaxNode = INFINITY, float scoreAllowedEdgeMax = INFINITY, bool flagUseZScoring = false);
    void scaleNodeScores();
    //void outputGraph(Graph const & g, string fileName);
    void resetScores();
    float transferScore(float score, float a = 0.5, float b = 0.1); 
    float calculateErrorAndUpdateScores();
    void updateNetwork();
    void outputGraph(Graph *g, string fileName);
    virtual void updateEdgeScore(int vIdSource, int vIdTarget);
    virtual void updateNodeScore(int vId);
    virtual void initializeIteration();
    virtual void finalizeIteration() {}; // = 0; // can not create an instance
    virtual void finalizeScoring();

public:
    ScoreNetwork(); 
    ScoreNetwork(string fileProtein, string fileInteraction, int typeTransfer = IDENTITY, bool flagUseReliabilityScore = true, bool flagAccumulate = true, bool flagResetSeedScoresToInitial = false, bool flagDebug = false, bool flagVerbose = false);
    virtual ~ScoreNetwork();
    void run(int nIteration, float tError = MIN_SCORE_SHIFT); 
    void printNetwork(string explanation = "");
    void printGraph(Graph const & g, string explanation = "");
};

#endif // SCORENETWORK_H_

