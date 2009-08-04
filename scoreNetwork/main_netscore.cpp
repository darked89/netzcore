
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


