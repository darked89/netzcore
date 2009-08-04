
#include "graph/globals.h"
#include "graph/Graph.h"
#include <ctime>

void graph_demo() {
    const int V = 1000;
    Graph network;
    MapIntToVertex *pMapNode = network.pMapIdToNode;
    MapIntToEdge *pMapEdge;
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToEdge::iterator itEdge, itEndEdge;
    hash_map<int,void*> mapIdProcessed;
    Vertex *pVertex;
    Edge *pEdge;

    clock_t t1 = clock();

    for(int i=0; i<V; i++) {
	network.addVertex();
    }
    for(int i=1; i<V+1; i++) {
        for(int j=i; j<V+1; j++) {
	    network.addEdge(i, j); 
	}
    }

    clock_t t2 = clock();
    std::cout << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << std::endl;

    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        pMapEdge = pVertex->pMapIdToEdge;
	for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
	    pEdge = itEdge->second;
	    if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
		    //std::cout << pEdge->idSource << pEdge->idTarget;
		    ;
	    }
	}
	mapIdProcessed[itNode->first] = NULL;
    }  

    clock_t t3 = clock();
    std::cout << (t3-t2) << " (" << (t3-t2)/(double)CLOCKS_PER_SEC << "s)" << std::endl;
}

int
main()
{
  graph_demo();
  return 0;
}

