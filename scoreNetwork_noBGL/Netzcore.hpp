#ifndef NETZCORE_H_
#define NETZCORE_H_

class Netzcore: public ScoreNetwork {
private:
    vector<Graph*> listNetworkSampled;
    void updateNodeScore(int vId, int nIteration);
    void updateEdgeScore(int vIdSource, int vIdTarget);
    float calculateNodeZScore(int vId, bool flagUseReliabilityScore, bool flagScaleZScoreWithNodeDegree);
    void updateSampledNetworkNodeScores();
    int convertScoreFloatToInt(float scoreFloat, unsigned int precision=PRECISION_COEFFICIENT);
    float convertScoreIntToFloat(int scoreInt, unsigned int precision=PRECISION_COEFFICIENT);

public:
    Netzcore(); 
    Netzcore(string fileProtein, string fileInteraction, int typeRandomization = TOPOLOGY_PRESERVING, int typeTransfer = IDENTITY, int nSample = N_SAMPLE, bool flagSampleOnce = true, bool flagUseReliabilityScore = true, bool flagNormalizeOnce = false, bool flagAccumulate = true, bool flagResetSeedScoresToInitial = false, bool flagUpdateSampledOnce = true, bool flagScaleZScoreWithNodeDegree = false, bool flagNetscore = false, bool flagSaveRandomNetworks = false, bool flagLoadRandomNetworks = false, bool fDebug = false, bool fVerbose = false);
    virtual ~Netzcore();
    void run(int nIteration, float tError = MIN_SCORE_SHIFT, int typeRandomization = TOPOLOGY_PRESERVING, int typeTransfer = IDENTITY, int nSample = N_SAMPLE, bool flagSampleOnce = true, bool flagUseReliabilityScore = true, bool flagNormalizeOnce = false, bool flagAccumulate = true, bool flagResetSeedScoresToInitial = false, bool flagUpdateSampledOnce = true, bool flagScaleZScoreWithNodeDegree = false, bool flagNetscore = false, bool flagSaveRandomNetworks = false, bool flagLoadRandomNetworks = false); 
    void printNetwork(string explanation = "");
    void printGraph(Graph const & g, string explanation = "");
};

#endif // NETZCORE_H_

