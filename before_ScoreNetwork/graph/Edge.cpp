#include "Edge.h"
#include "globals.h"

Edge::Edge() {
	idSource = 0;
	idTarget = 0;
}

Edge::Edge(int sId, int eId) {
	idSource = sId;
	idTarget = eId;
}

Edge::Edge(int sId, int eId, Data eData) {
	idSource = sId;
	idTarget = eId;
	data = eData;
}

Edge::Edge(Edge const & e) {
	idSource = e.idSource;
	idTarget = e.idTarget;
	data = e.data;
}

Edge::~Edge(){
}

ostream& operator<<(ostream& output, const Edge& e) {
    output << "(" <<  e.idSource << " - " << e.idTarget << ", " << e.data.score << ")";
    return output;  // for multiple << operators.
}
