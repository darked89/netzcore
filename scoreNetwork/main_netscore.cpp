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

//#include "ScoreNetwork.hpp"
#include "Netscore.hpp"

#include <iostream>
#include <ctime>

using namespace std;

void runNetscore();

int main() 
{
    clock_t t1 = clock();
    runNetscore();
    clock_t t2 = clock();
    cout << "Time: " << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;
    return 0;
}

void runNetscore() 
{
    clock_t t1 = clock();
    //Netscore sN("../../data/input/node_scores.txt", "../../data/input/edge_weights.txt");
    Netscore sN("../../data/toy_data/test_proteins_small.txt", "../../data/toy_data/test_interactions_small.txt", "../../data/output/test.txt", true, true, false, true);
    clock_t t2 = clock();
    cout << "Time to load graph: " << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;
    sN.run(1,3);
    //cout << sN.run(1,1) << endl;
}


