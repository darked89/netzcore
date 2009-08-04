#include "Vertex.h"

Vertex::Vertex() {
	id = 0;
	degree = 0;	
	pMapIdToEdge = new MapIntToEdge;
	pMapIdToScore = new MapIntToFloatInt;
}

Vertex::Vertex(int vId) {
	id = vId;
	degree = 0;
	pMapIdToEdge = new MapIntToEdge;
	pMapIdToScore = new MapIntToFloatInt;
}

Vertex::Vertex(int vId, Data vData) {
	id = vId;
	data = vData;
	degree = 0;
	pMapIdToEdge = new MapIntToEdge;
	pMapIdToScore = new MapIntToFloatInt;
}

Vertex::Vertex(Vertex const & v) {
	id = v.id;
	data = v.data;
	degree = v.degree;
	pMapIdToEdge = new MapIntToEdge;
	MapIntToEdge *pMapEdge = v.pMapIdToEdge;
	MapIntToEdge::iterator it = pMapEdge->begin(), itEnd = pMapEdge->end();
	for(; it != itEnd; ++it) {
		(*pMapIdToEdge)[it->first] = new Edge(*(it->second));
	}	 
}

Vertex::~Vertex() {
	// iterate over Edges and remove them
	MapIntToEdge::iterator it = pMapIdToEdge->begin(), itEnd = pMapIdToEdge->end();
	for(; it != itEnd; ++it)
		delete it->second;
	// removeMap
	delete pMapIdToEdge;
}

bool Vertex::containsNeighbor(int id) {
	MapIntToEdge::iterator it;
	it = pMapIdToEdge->find(id);
	if(it != pMapIdToEdge->end()) {
	    return true;
	}
	return false;
}

bool Vertex::addNeighbor(int vId) {
	if(!containsNeighbor(vId)) {
		(*pMapIdToEdge)[vId] = new Edge(id, vId);
		degree++;
		return true;
	} else {
		cerr << "Warning: Neighbor already contained (ignoring) " << endl;
		return false;
	}
}

bool Vertex::addNeighbor(int vId, Data eData) {
	if(!containsNeighbor(vId)) {
		(*pMapIdToEdge)[vId] = new Edge(id, vId, eData);
		degree++;
		return true;
	} else {
		cerr << "Warning: Neighbor already contained (ignoring) " << endl;
		return false;
	}
}

void Vertex::removeNeighbor(int vId){
	MapIntToEdge::iterator it;
	it = pMapIdToEdge->find(vId);
    if(it != pMapIdToEdge->end()) {
    	delete it->second;
	    pMapIdToEdge->erase(it);
	    degree--;
    } else {
        cerr << "Warning: trying to remove a non-existant edge" << endl;
    }
}

ostream& operator<<(ostream& output, const Vertex& v) {
    output << "(" <<  v.id << "{" << v.data.score << "}: ";
	MapIntToEdge::iterator it, itEnd;
	MapIntToEdge *map = v.pMapIdToEdge;
	for(it = map->begin(), itEnd = map->end(); it != itEnd; ++it)
		output << it->first << ", ";
	output << ")";
    return output; 
}

