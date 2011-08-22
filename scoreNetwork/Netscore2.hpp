/*
    GUILD (Genes Underlying Inheritance Linked Disorders) implements several 
    graph based algorithms for scoring relevance of a node in the network in 
    terms of a phenotype using known associations in the node's neighborhood 
    for that phenotype. GUILD has been applied to the prioritization of genes 
    for several human disorders. 2011 - Emre Guney (Unviersitat Pompeu Fabra)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NETSCORE2_HPP_9054862562754458964905368042156
#define NETSCORE2_HPP_9054862562754458964905368042156

#include "graph/Graph.hpp"

#include "Exceptions.hpp"

//#define CONTROL 0
#define sq(x) (x*x)


class Netscore2 {
public:
    Netscore2(); 
    Netscore2(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeReliabilityScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    ~Netscore2();
    void run(int nRepetition, int nIteration, float tError = PRECISION); 
    Graph & getNetwork() { return network; };

private:
    // MEMBERS
    Graph network;
    int nIteration;
    int iterationCounter;
    float minScore;
    float maxScore;
    bool flagUseEdgeScore;
    bool flagAccumulateToInitialNodeScore;
    bool flagResetSeedScoresToInitial;
    bool flagVerbose;
    std::string outputFile;

    void updateEdgeScore(Edge e) {}; 
    void finalizeScoring() {};
    void finalizeRepetition() {};
    void initializeIteration() {};
    void finalizeIteration() {}; 

    float calculateErrorAndUpdateScores();
    void updateNetwork();
    void scaleNodeScores();
    float transferScore(float score, float a = 0.5, float b = 0.1); 

    static const float PRECISION = 0.00001;
    enum TransferType { IDENTITY, POLYNOMIAL, LOGARITHMIC, EXPONENTIAL };
    int typeTransfer;
    
 
    void initializeScoring();
    void initializeRepetition();
    void updateNodeScore(Vertex v);
};

#endif // NETSCORE2_HPP_

