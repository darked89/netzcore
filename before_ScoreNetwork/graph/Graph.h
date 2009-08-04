#ifndef GRAPH_H_
#define GRAPH_H_

#include "globals.h"
#include "hash.h"
#include "Vertex.h"
#include "Edge.h"

typedef hash_map<int, Vertex*> MapIntToVertex; // hash_map<char*, int, hash<char*>, eqstr> Map;
//typedef hash_map<string, int, hash<string>, isEqualString> MapStringToInt; 
typedef hash_map<const char*, int, hash<const char*>, eqstr> MapStringToInt;
//typedef MapIntToVertex::value_type value_type;
//typedef MapIntToVertex::iterator iterator;
typedef hash_map<int, string> MapIntToString;
typedef hash_map<int, int> MapIntToInt;

class Graph {
friend ostream& operator<<(ostream& output, const Graph& g);
private:
void removeEdgeDirected(int vIdSource, int vIdTarget);
public:
	MapIntToVertex *pMapIdToNode; // for mapping of node ids to Vertices
	MapStringToInt mapNameToId; // for mapping of node names to node ids
	MapIntToString mapIdToName; // for mapping of node ids to node names
	MapIntToInt * pMapIdNewToId; // for mapping ids in sampled graph
	MapIntToInt * pMapIdToIdNew; // for mapping ids in sampled graph
	int idLast; // id assigned to last added node
    bool flagDirected;
    unsigned int nNode;
    unsigned int nEdge;    
    bool flagSampled;
//public:
	Graph();
    Graph(bool fDirected);
    Graph(Graph const & g, bool flagNoEdges = false);
	virtual ~Graph();
    vector<int> getNodeIdList();
    hash_map<int, vector<int> *> getMapDegreeToNodeIdList();
    hash_map<int, hash_map<int, void*> *> getMapDegreeToNodeIdMap();
    hash_map<int, int> getDegreeHistogram();
    vector<Data> getEdgeDataList();
	void initializeMapNewId();
	bool isSampled() const;
	void swapInMapNewId(int idNode, int idNodeNew);
	void convertMapNewIdToMapId();
	bool containsVertex(string name);
	bool containsVertex(int id);
	int addVertex();
	int addVertex(string name);
	void addVertex(int id, Data vData);
	int addVertex(string name, Data vData);
	bool addEdge(int idSource, int idTarget);
	bool addEdge(int idSource, int idTarget, Data eData);
	bool addEdge(string nameSource, string nameTarget);
	bool addEdge(string nameSource, string nameTarget, Data eData);
	void removeVertex(int vId);
	void removeEdge(int vIdSource, int vIdTarget); 
};

#endif /*GRAPH_H_*/
