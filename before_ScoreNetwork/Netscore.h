#ifndef NETSCORE_H_
#define NETSCORE_H_

#include "graph/globals.h"
#include "graph/Graph.h"
#include <cmath> // for INFINITY

#define CONTROL 0

#define IDENTITY 1
#define POLYNOMIAL 2
#define LOGARITHMIC 3
#define EXPONENTIAL 4

#define STEPWISE_ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION 10
#define ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_DIVISION 20
#define ACCUMULATING_CONTRIBUTION_LINKER_PARAMETER_CHECK 30
#define STEPWISE_CONTRIBUTION_LINKER_PARAMETER_CHECK 40
#define ACCUMULATING_CONTRIBUTION_LINKER_SCORE_CHECK_DEGREE_DIVISION 50
#define STEPWISE_CONTRIBUTION_LINKER_SCORE_CHECK_DEGREE_DIVISION 60
#define ACCUMULATING_CONTRIBUTION_DEGREE_DIVISION 70

#define MIN_SCORE_AFTER_SHIFT 0.0001

//class Netscore: public Graph {
class Netscore {
private:
    Graph network;
	bool flagDebug;
	hash_map<int, float> mapIdToCCoefficient;
	hash_map<int, int> mapIdToClusterId;
	hash_map<int, vector<int> > mapClusterIdToId;
    float transferScore(int type, float score, float a = 0.5, float b = 0.1); 
    float calculateError(int type);
    //float decideErrorScoreNode(int type, Vertex *pVertex);
    //float decideErrorScoreEdge(int type, Edge *pEdge);
    void updateNodeScore(int type, int typeTransfer, int vId, int tDegreeLinkerNodeMin, int tDegreeLinkerNodeMinIncrease, float scoreAllowedMaxNode, float tScoreNeighborMin, float tScoreNeighborMax);
    //float decideInitialScoreNode(int type, Vertex *pVertex);
    float decideStepScoreNode(int type, Edge *pEdge, Vertex *pNeighbor, int tDegreeLinkerNodeMin, float tScoreNeighborMin, float tScoreNeighborMax, int *nCount);
    float decideScoreNode(int type, Vertex *pVertex, float sumScore, int nCount, int tDegreeLinkerNodeMin);
    void updateEdgeScore(int type, int typeTransfer, int idSource, int idTarget, int tDegreeLinkerEdgeMin, int tDegreeLinkerEdgeMinIncrease, float scoreAllowedMaxEdge, float tScoreNeighborMin, float tScoreNeighborMax);
    void updateDualScore(int vId, int tDegreeLinkerNodeMinIncrease, float scoreAllowedMaxNode);
    //float decideInitialScoreEdge(int type, Edge *pEdge);
    float decideStepScoreEdge(int type, Vertex *pNeighbor, Edge *pEdge, int tDegreeLinkerNodeMin, float tScoreNeighborMin, float tScoreNeighborMax, int *eCount);
    float decideScoreEdge(int type, Edge *pEdge, float sumScore, int eCount, int tDegreeLinkerEdgeMin);
    int combination2(int x);
    void calculateClusteringCoefficient();
    void clusterNetwork();
    void expandCluster(int vId, float cCoefficient, int idCluster, hash_map<int, void*> *pMapIdProcessed);
    void calculateClusterBasedNodeZScores();
    void calculateNodeZScores(int nSample);
    void sampleNetwork(int nSample, hash_map<int, vector<int> > *pMapIdToListDegree);
    float calculateScoreDifference();
public:
    Netscore(); 
    Netscore(string fileProtein, string fileInteraction, bool flagDebug = false);
    virtual ~Netscore();
    void loadProteins(string const &fileName);
    void loadInteractions(string const &fileName);
    //void filterNetwork(int degreeNodeMin, int degreeNodeMax, int scoreNodeMin, int scoreNodeMax, int scoreEdgeMin, int scoreEdgeMax);
    void filterNetwork(int degreeNodeMin = -INT_INFINITY, int degreeNodeMax = INT_INFINITY, float scoreNodeMin = -INFINITY, float scoreNodeMax = INFINITY, float scoreEdgeMin = -INFINITY, float scoreEdgeMax = INFINITY);
    void scaleScores(bool flagUseZScoring = false); //(float scoreAllowedMaxNode = INFINITY, float scoreAllowedEdgeMax = INFINITY, bool flagUseZScoring = false);
    void resetScores();
    void run(int type, int typeTransfer, int nIteration, float tError, int tDegreeLinkerNodeMin = 1, int tDegreeLinkerNodeMinIncrease = 1, int tDegreeLinkerEdgeMin = 1, int tDegreeLinkerEdgeMinIncrease = 1, float scoreMaxAllowedNode = 1.0, float scoreMaxAllowedEdge = 1.0, float tScoreNeighborMin = -INFINITY, float tScoreNeighborMax = INFINITY); 
    //void updateNetwork(int nIteration, float tError);
    void updateNetwork(int type, int typeTransfer, int tDegreeLinkerNodeMin, int tDegreeLinkerNodeMinIncrease, int tDegreeLinkerEdgeMin, int tDegreeLinkerEdgeMinIncrease, float scoreMaxAllowedNode, float scoreMaxAllowedEdge, float tScoreNeighborMin, float tScoreNeighborMax);
    void printNetwork(string explanation = "");
    void printNodes(string explanationPrefix = "");
    void printEdges();
    float getEdgeScore(string vNameSource, string vNameTarget);
};

#endif // NETSCORE_H_

