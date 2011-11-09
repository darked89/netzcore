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

#include "FunctionalFlow.h"
#include <cstdio>

void functionalFlow(string fileProtein, string fileInteraction, string fileOutput, int nIteration, float seedScoreThreshold) {
	// flagRemoveSeedEffect, flagDebug
	FunctionalFlow fFlow = FunctionalFlow(fileProtein, fileInteraction, fileOutput, seedScoreThreshold, false, false);
	fFlow.run(nIteration);
	//fFlow.printNetwork();
	//fFlow.printNodesAndScores();
	fFlow.outputScores();
}

int main(int argc, char **argv)
{
	int nIteration; 
	float threshold;
	sscanf(argv[4], "%d", &nIteration);
	sscanf(argv[5], "%f", &threshold);
	//cout << "Proteins: " << argv[1] << " Interactions: " << argv[2] << " i: " << nIteration << endl;
	functionalFlow(argv[1], argv[2], argv[3], nIteration, threshold);
	return 0;
}

