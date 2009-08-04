#ifndef SCORENETWORK_H_
#define SCORENETWORK_H_

#include "Graph.hpp"
#include <hash_set>

#include <cmath> // for INFINITY

//#define CONTROL 0


/** BEGIN STRING HASHING METHODS **/
struct isEqualString
{
  bool operator()(std::string s1, std::string s2) const
  {
    return s1.compare(s2) == 0;
  }
};

namespace std
{                                                           
    template<> struct hash< std::string >    
    {
	size_t operator()( const std::string& x ) const
	{
	    return hash< const char* >()( x.c_str() ); 
	}
    };                      
}          

/** END STRING HASHING METHODS **/

using namespace std;

class ScoreNetwork {
private:
    // CONSTANTS (MACROS)
    //static const int MIN_SCORE = 0.0;
    //static const int MAX_SCORE = 1.0;
    static const enum { IDENTITY, POLYNOMIAL, LOGARITHMIC, EXPONENTIAL };
    
    // MEMBERS
    Graph network;
    hash_set<string> mapSource;
    int nIteration;
    int iterationCounter;

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
    //void filterNetwork(int degreeNodeMin = -INT_INFINITY, int degreeNodeMax = INT_INFINITY, float scoreNodeMin = -INFINITY, float scoreNodeMax = INFINITY, float scoreEdgeMin = -INFINITY, float scoreEdgeMax = INFINITY, bool flagRemoveSelfEdges = true);
    //void scaleScores(float scoreAllowedMaxNode = INFINITY, float scoreAllowedEdgeMax = INFINITY, bool flagUseZScoring = false);
    void scaleNodeScores();
    //void outputGraph(Graph const & g, string fileName);
    //void resetScores();
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

