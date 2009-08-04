#ifndef VERTEX_H_
#define VERTEX_H_

#include "hash.h"
#include "globals.h"
#include "Edge.h"

typedef hash_map<int, Edge*> MapIntToEdge;
//typedef MapIntToEdge::value_type value_type;
//typedef MapIntToEdge::iterator iterator;
typedef hash_map<int, pair<float, int> > MapIntToFloatInt;

class Vertex
{
friend ostream& operator<<(ostream& output, const Vertex& v);
public:
	int id;
	Data data;
	int degree;
	//float cClustering;
	MapIntToEdge *pMapIdToEdge; // neighbor id to edge
        MapIntToFloatInt *pMapIdToScore;
//public:
	Vertex();
	Vertex(int id);
	Vertex(int vId, Data vData);
	Vertex(Vertex const & v);
	virtual ~Vertex();
	bool containsNeighbor(int id);
	bool addNeighbor(int vId);
	bool addNeighbor(int vId, Data eData);
	void removeNeighbor(int vId);
};

#endif /*VERTEX_H_*/
