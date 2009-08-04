#include "Graph.h"

Graph::Graph() {
	pMapIdToNode = new MapIntToVertex;
	pMapIdNewToId = NULL;
	pMapIdToIdNew = NULL;
	idLast = 0;
    flagDirected = false;
    nNode = 0;
    nEdge = 0;
    //flagSampled = false; // check nullity of pMapIdNewToId instead 
}

Graph::Graph(bool fDirected) {
	pMapIdToNode = new MapIntToVertex;
	pMapIdNewToId = NULL;
	pMapIdToIdNew = NULL;
	idLast = 0;
    flagDirected = fDirected;
    nNode = 0;
    nEdge = 0;
    //flagSampled = false;
}

Graph::Graph(Graph const & g, bool flagNoEdges) {
	mapNameToId = g.mapNameToId;
	mapIdToName = g.mapIdToName;
	idLast = g.idLast;
    flagDirected = g.flagDirected;
    nNode = g.nNode;
    if(flagNoEdges)
        nEdge = 0;
    else
        nEdge = g.nEdge;
    pMapIdToNode = new MapIntToVertex;
    MapIntToVertex::iterator it = (g.pMapIdToNode)->begin(), itEnd = (g.pMapIdToNode)->end();
	for(; it != itEnd; ++it) {
        if(flagNoEdges)
		    (*pMapIdToNode)[it->first] = new Vertex((it->second)->id, (it->second)->data);
        else
		    (*pMapIdToNode)[it->first] = new Vertex(*(it->second));
	}
	MapIntToInt *pMapId;
	MapIntToInt::iterator itId, itIdEnd;
	pMapId = g.pMapIdNewToId;
	if(pMapId == NULL) {
		pMapIdNewToId = NULL;
	} else {
		pMapIdNewToId = new MapIntToInt;
		for(itId=pMapId->begin(), itIdEnd=pMapId->end(); itId!=itIdEnd; ++itId) {
			(*pMapIdNewToId)[itId->first] = itId->second;
		}	
	}
	pMapId = g.pMapIdToIdNew;
	if(pMapId == NULL) {
		pMapIdToIdNew = NULL;
	} else {
		pMapIdToIdNew = new MapIntToInt;
		for(itId=pMapId->begin(), itIdEnd=pMapId->end(); itId!=itIdEnd; ++itId) {
			(*pMapIdToIdNew)[itId->first] = itId->second;
		}
	}
	//flagSampled = false;
}

Graph::~Graph() {
	MapIntToVertex::iterator it = pMapIdToNode->begin(), itEnd = pMapIdToNode->end();
	for(; it != itEnd; ++it)
		delete it->second;
	mapIdToName.clear();
	mapNameToId.clear();
	///*
	if(pMapIdNewToId != NULL) {
		delete pMapIdNewToId;
		pMapIdNewToId = NULL;
		//cout << "deleted" << endl;
	}
	if(pMapIdToIdNew != NULL) {
		delete pMapIdToIdNew;
		pMapIdToIdNew = NULL;
		//cout << "deleted" << endl;
	}
	//*/
}

vector<int> Graph::getNodeIdList() {
	MapIntToVertex::iterator itNode, itEndNode;
	vector<int> listIdNode;
	for(itNode = pMapIdToNode->begin(), itEndNode = pMapIdToNode->end(); itNode != itEndNode; ++itNode) {
		listIdNode.push_back(itNode->first);
    }
    return listIdNode;
}

hash_map<int, vector<int> *> Graph::getMapDegreeToNodeIdList() {
	MapIntToVertex::iterator itNode, itEndNode;
    Vertex *pVertex;
    hash_map<int, vector<int> *> mapDegreeToListIdNode;
    hash_map<int, vector<int> *>::iterator itDegree;

    /// for 0 degreed nodes that could be added
    mapDegreeToListIdNode[0] = new vector<int>;
	for(itNode = pMapIdToNode->begin(), itEndNode = pMapIdToNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        if(mapDegreeToListIdNode.find(pVertex->degree) == mapDegreeToListIdNode.end()) {
            mapDegreeToListIdNode[pVertex->degree] = new vector<int>(1, pVertex->id);
        } else {
            (mapDegreeToListIdNode[pVertex->degree])->push_back(pVertex->id);
        }
    }
    return mapDegreeToListIdNode;
}

hash_map<int, hash_map<int, void*> *> Graph::getMapDegreeToNodeIdMap() {
	MapIntToVertex::iterator itNode, itEndNode;
    Vertex *pVertex;
    hash_map<int, hash_map<int, void*> *> mapDegreeToMapIdNode;
    hash_map<int, void*> *pMapIdNodeTemp;
    hash_map<int, hash_map<int, void*> *>::iterator itDegree;
    /// for handling nodes with a degree of 1
    //pMapIdNodeTemp = new hash_map<int,void*>; instead add when required
    //mapDegreeToMapIdNode[0] = pMapIdNodeTemp; 
	for(itNode = pMapIdToNode->begin(), itEndNode = pMapIdToNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        if(mapDegreeToMapIdNode.find(pVertex->degree) == mapDegreeToMapIdNode.end()) {
            pMapIdNodeTemp = new hash_map<int,void*>;
            (*pMapIdNodeTemp)[pVertex->id] = NULL;
            mapDegreeToMapIdNode[pVertex->degree] = pMapIdNodeTemp; 
        } else {
            (*(mapDegreeToMapIdNode[pVertex->degree]))[pVertex->id] = NULL;
        }
    }
    return mapDegreeToMapIdNode;
}

hash_map<int, int> Graph::getDegreeHistogram() {
	MapIntToVertex::iterator itNode, itEndNode;
    Vertex *pVertex;
    hash_map<int, int> mapDegreeToCount;
	for(itNode = pMapIdToNode->begin(), itEndNode = pMapIdToNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        if(mapDegreeToCount.find(pVertex->degree) == mapDegreeToCount.end()) {
            mapDegreeToCount[pVertex->degree] = 1;
        } else {
            mapDegreeToCount[pVertex->degree] += 1;
        }
    }
    return mapDegreeToCount;
}

vector<Data> Graph::getEdgeDataList() {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge; 
	MapIntToEdge *pMapEdge;
    Vertex *pVertex;
	Edge *pEdge;
    vector<Data> listDataEdge;
    hash_map<int,void*> mapIdProcessed;
	for(itNode = pMapIdToNode->begin(), itEndNode = pMapIdToNode->end(); itNode != itEndNode; ++itNode) {
        pVertex = itNode->second;
        pMapEdge = pVertex->pMapIdToEdge;
        for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
	    	if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
                listDataEdge.push_back(pEdge->data);
            }
        }
	    mapIdProcessed[itNode->first] = NULL;
    }
    return listDataEdge;
}

void Graph::initializeMapNewId() {
	MapIntToVertex::iterator it, itEnd;
	pMapIdNewToId = new MapIntToInt;
	for(it=pMapIdToNode->begin(), itEnd=pMapIdToNode->end(); it!=itEnd; ++it) {
		(*pMapIdNewToId)[it->first] = it->first;
	}
	//flagSampled = true;
}

bool Graph::isSampled() const {
	if(pMapIdNewToId == NULL) {
		return false;
	} else {
		return true;
	}
}

void Graph::swapInMapNewId(int idNode, int idNodeNew){
	MapIntToInt::iterator itIdNode, itIdNodeNew;
	MapIntToInt *pMapId;
	int posNode, posNodeNew;
	pMapId = pMapIdNewToId;
	itIdNode = pMapId->find(idNode);
	posNode = itIdNode->second;
	itIdNodeNew = pMapId->find(idNodeNew);
	posNodeNew = itIdNodeNew->second;
	pMapId->erase(itIdNode);
	(*pMapId)[idNode] = posNodeNew;
	pMapId->erase(itIdNodeNew);
	(*pMapId)[idNodeNew] = posNode;
	//for(itIdNode=pMapIdNewToId->begin(); itIdNode != pMapIdNewToId->end(); ++itIdNode) {
	//	cout << itIdNode->second << "/" << itIdNode->first << endl;
	//}
}

void Graph::convertMapNewIdToMapId() {
	MapIntToInt::iterator itIdNode, itIdNodeEnd;
	pMapIdToIdNew = new MapIntToInt; 
	for(itIdNode=pMapIdNewToId->begin(), itIdNodeEnd=pMapIdNewToId->end(); itIdNode != itIdNodeEnd; ++itIdNode) {
		(*pMapIdToIdNew)[itIdNode->second] = itIdNode->first;
		//cout << itIdNode->second << "/" << itIdNode->first << endl;
	}
}

bool Graph::containsVertex(string name) {
	MapStringToInt::iterator iterator;
	iterator = mapNameToId.find(name.c_str());
	if(iterator != mapNameToId.end()) {
	    return true;
	}
	return false;
}

bool Graph::containsVertex(int id) {
	MapIntToString::iterator iterator;
	iterator = mapIdToName.find(id);
	if(iterator != mapIdToName.end()) {
	    return true;
	}
	return false;
}

int Graph::addVertex() {
	idLast++;
	(*pMapIdToNode)[idLast] = new Vertex(idLast);
	nNode++;
    return idLast;
}

int Graph::addVertex(string name) {
	if(!containsVertex(name)) {
		idLast++;
		(*pMapIdToNode)[idLast] = new Vertex(idLast);
		mapNameToId[name.c_str()] = idLast; 
		mapIdToName[idLast] = name;
		nNode++;
		return idLast;
	}
	cerr << "Warning: Vertex already contained (ignoring): " << name << endl;
	return mapNameToId[name.c_str()];
}

void Graph::addVertex(int id, Data data) {
	if(!containsVertex(id)) { 
		(*pMapIdToNode)[id] = new Vertex(id, data);
		nNode++;
	} else {
		cerr << "Warning: Vertex already contained (ignoring): " << id << endl;
	}
}

int Graph::addVertex(string name, Data data) {
	if(!containsVertex(name)) {
		idLast++;
	    (*pMapIdToNode)[idLast] = new Vertex(idLast, data); 
		mapNameToId[name.c_str()] = idLast; 
		mapIdToName[idLast] = name;
		nNode++;
		return idLast;
	}
	cerr << "Warning: Vertex already contained (ignoring): " << name << endl;
	return mapNameToId[name.c_str()];
}

bool Graph::addEdge(int idSource, int idTarget) {
	bool valReturn = false;
	if(containsVertex(idSource) && containsVertex(idTarget)) {
		valReturn = (*pMapIdToNode)[idSource]->addNeighbor(idTarget);
	    if(!flagDirected and idSource != idTarget)
		    (*pMapIdToNode)[idTarget]->addNeighbor(idSource);
	    nEdge++;
	} else {
		cerr << "Warning: At least one of the vertex is missing (ignoring edge): " << idSource << ", " << idTarget << endl;
		if(!containsVertex(idSource))
			cerr << "Warning: vertex is missing: " << idSource << endl;
 		if(!containsVertex(idTarget))
			cerr << "Warning: vertex is missing: " << idTarget << endl;
	}
	return valReturn;
}

bool Graph::addEdge(int idSource, int idTarget, Data eData) {
	bool valReturn = false;
	if(containsVertex(idSource) && containsVertex(idTarget)) {
		valReturn = (*pMapIdToNode)[idSource]->addNeighbor(idTarget, eData);
	    if(!flagDirected and idSource != idTarget)
		    (*pMapIdToNode)[idTarget]->addNeighbor(idSource, eData);
	    nEdge++;
	} else {
		cerr << "Warning: At least one of the vertex is missing (ignoring edge): " << idSource << ", " << idTarget << endl;
		if(!containsVertex(idSource))
			cerr << "Warning: vertex is missing: " << idSource << endl;
 		if(!containsVertex(idTarget))
			cerr << "Warning: vertex is missing: " << idTarget << endl;
	}
	return valReturn;
}

bool Graph::addEdge(string nameSource, string nameTarget) {
	bool valReturn = false;
	if(containsVertex(nameSource) && containsVertex(nameTarget)) {
	    int idSource, idTarget;
	    idSource = mapNameToId[nameSource.c_str()];
        idTarget = mapNameToId[nameTarget.c_str()];
        valReturn = (*pMapIdToNode)[idSource]->addNeighbor(idTarget);
	    if(!flagDirected and idSource != idTarget)
		    (*pMapIdToNode)[idTarget]->addNeighbor(idSource);
	    nEdge++;
	} else {
		cerr << "Warning: At least one of the vertex is missing (ignoring edge): " << nameSource << ", " << nameTarget << endl;
		if(!containsVertex(nameSource))
			cerr << "Warning: vertex is missing: " << nameSource << endl;
 		if(!containsVertex(nameTarget))
			cerr << "Warning: vertex is missing: " << nameTarget << endl;
	}
	return valReturn;
}

bool Graph::addEdge(string nameSource, string nameTarget, Data eData) {
	bool valReturn = false;
	if(containsVertex(nameSource) && containsVertex(nameTarget)) {
	    int idSource, idTarget;
	    idSource = mapNameToId[nameSource.c_str()];
        idTarget = mapNameToId[nameTarget.c_str()];
        valReturn = (*pMapIdToNode)[idSource]->addNeighbor(idTarget, eData);
	    if(!flagDirected and idSource != idTarget)
		    (*pMapIdToNode)[idTarget]->addNeighbor(idSource, eData);
	    nEdge++;
	} else {
		cerr << "Warning: At least one of the vertices is missing " << nameSource << " " << nameTarget << endl;
	}
	return valReturn;
}

void Graph::removeVertex(int vId) {
    Vertex *pVertex = (*pMapIdToNode)[vId];
    if(!flagDirected) {
        // remove neighbors edges to this edge
    	MapIntToEdge::iterator it = pVertex->pMapIdToEdge->begin(), itEnd = pVertex->pMapIdToEdge->end();
    	for(; it != itEnd; ++it) {
            removeEdgeDirected(it->first, vId); // (*pMapIdToNode)[it->first]->removeNeighbor(vId); // removeEdge(vId, it->first); 
            nEdge--; //!! should be erased if removeEdge is used. 
        }
    }
    // remove node from the graph
    delete pVertex; // does erase call destructor somehow maybe with delete??
    pMapIdToNode->erase(vId);
    string name = mapIdToName[vId];
	mapNameToId.erase(name.c_str());
	mapIdToName.erase(vId);
	nNode--;
	//cout << " size: " << pMapIdToNode->size() << endl;
	//cout<< "removeVertex " << *this << endl;
}

void Graph::removeEdge(int vIdSource, int vIdTarget) {
	removeEdgeDirected(vIdSource, vIdTarget);
	if(!flagDirected and vIdSource != vIdTarget) {
		removeEdgeDirected(vIdTarget, vIdSource);
	}
	nEdge--;
	// Vertex *pVertex = (*pMapIdToNode)[vId];
	// if(pVertex->pMapIdToEdge->size() == 0) removeVertex(vId) 
}

// nEdge is not modified since this method supposed to be called only from removeEdge or removeVertex 
void Graph::removeEdgeDirected(int vIdSource, int vIdTarget) {
	Vertex *pVertex = (*pMapIdToNode)[vIdSource];
	pVertex->removeNeighbor(vIdTarget);
	// cout << *pVertex << " removed edge with: " << vIdTarget << endl;
}

ostream& operator<<(ostream& output, const Graph& g) {
	MapIntToVertex *pMap = g.pMapIdToNode;
	MapIntToVertex::iterator it, itEnd;
    MapIntToString map2 = g.mapIdToName;
	itEnd = pMap->end(); 
	for(it = pMap->begin(); it != itEnd; ++it)
		output << "[" << map2[it->first] << "]" << *(it->second) << " - ";
    output << endl;
    return output;  
}

