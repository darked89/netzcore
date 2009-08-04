#ifndef NETZCORE_H_
#define NETZCORE_H_

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

class Netzcore {
private:
    Graph network;
    vector<Graph*> listNetworkSampled;
    hash_map<int, void*> mapSource;
    float minZScore;
    float maxZScore;
    float minScore;
    float maxScore;
    bool flagDebug;
    bool flagVerbose;
    void loadProteins(string const &fileName);
    void loadInteractions(string const &fileName);    
    void filterNetwork(int degreeNodeMin = -INT_INFINITY, int degreeNodeMax = INT_INFINITY, float scoreNodeMin = -INFINITY, float scoreNodeMax = INFINITY, float scoreEdgeMin = -INFINITY, float scoreEdgeMax = INFINITY, bool flagRemoveSelfEdges = true);
    //void scaleScores(float scoreAllowedMaxNode = INFINITY, float scoreAllowedEdgeMax = INFINITY, bool flagUseZScoring = false);
    void resetScores();
    float transferScore(int type, float score, float a = 0.5, float b = 0.1); 
    float calculateErrorAndUpdateScores(bool flagResetSeedScoresToInitial);
    void scaleNodeScores();
    void updateNodeScore(int typeTransfer, int vId, int nIteration, bool flagUseReliabilityScore, bool flagNoZScoring, bool flagAccumulate, bool flagScaleZScoreWithNodeDegree);
    void updateEdgeScore(int typeTransfer, int vIdSource, int vIdTarget);
    void sampleNetwork(int nSample, int typeRandomization);
    void writeSampledNetworks();
    void readSampledNetworks(int nSample);
    //void outputGraph(Graph const & g, string fileName);
    void outputGraph(Graph *g, string fileName);
    Graph * readSampledNetwork(string fileName);
    float calculateNodeZScore(int vId, bool flagUseReliabilityScore, bool flagScaleZScoreWithNodeDegree);
    void updateSampledNetworkNodeScores();
    void deleteSampledNetworks();
    void updateNetwork(int typeTransfer, int nIteration, bool flagUseReliabilityScore, bool flagNoZScoring, bool flagAccumulate, bool flagScaleZScoreWithNodeDegree);
    
    int convertScoreFloatToInt(float scoreFloat, unsigned int precision=PRECISION_COEFFICIENT);
    float convertScoreIntToFloat(int scoreInt, unsigned int precision=PRECISION_COEFFICIENT);

public:
    Netzcore(); 
    Netzcore(string fileProtein, string fileInteraction, bool flagDebug = false, bool flagVerbose = false);
    virtual ~Netzcore();
    void run(int nIteration, float tError = MIN_SCORE_SHIFT, int typeRandomization = TOPOLOGY_PRESERVING, int typeTransfer = IDENTITY, int nSample = N_SAMPLE, bool flagSampleOnce = true, bool flagUseReliabilityScore = true, bool flagNormalizeOnce = false, bool flagAccumulate = true, bool flagResetSeedScoresToInitial = false, bool flagUpdateSampledOnce = true, bool flagScaleZScoreWithNodeDegree = false, bool flagNetscore = false, bool flagSaveRandomNetworks = false, bool flagLoadRandomNetworks = false); 
    void printNetwork(string explanation = "");
    void printGraph(Graph const & g, string explanation = "");
    void printSampledNetworks();
};

#endif // NETZCORE_H_

