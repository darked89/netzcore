#ifndef EDGE_H_
#define EDGE_H_

#include "globals.h"

class Edge {
friend ostream& operator<<(ostream& output, const Edge& e);
public:
	int idSource;
	int idTarget;
	Data data;
	Edge();
	Edge(int sId, int eId);
	Edge(int sId, int eId, Data eData);
	Edge(Edge const & e);
	virtual ~Edge();
};

#endif /*EDGE_H_*/
