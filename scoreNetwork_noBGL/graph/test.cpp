#include "Edge.h"
#include "Vertex.h"
#include "Graph.h"

void test() {
	
	Vertex *pV1 = new Vertex(1, Data(2.3));
	Vertex *pV2 = new Vertex(2, Data(5.3));
	Vertex v1 = Vertex(3, Data(5.6));
	
	Edge *pE1 = new Edge(1, 2, Data(0.3));
	Edge e2 = Edge(3, 2, Data(3.6));
	
	cout << *pV1 << " " << *pV2 << " " << v1 << " " << endl; 
	cout << *pE1 << " " << e2 << " " << endl;
	
	pV1->addNeighbor(2, Data(0.3));
	pV2->addNeighbor(1, Data(0.3));
	v1.addNeighbor(2, Data(3.6));
	pV2->addNeighbor(3, Data(3.6));

	cout << *pV1 << " " << *pV2 << " " << v1 << " " << endl;
	
	pV2->removeNeighbor(1);
	pV2->removeNeighbor(3);
	cout << *pV1 << " " << *pV2 << " " << v1 << " " << endl;

    Graph graph = Graph();
    graph.addVertex("v1", Data(1.2));
    graph.addVertex("v1", Data(1.2));
    graph.addVertex("v2", Data(2.3));
    graph.addVertex("v3", Data(3.4));
    graph.addVertex("v4", Data(4.5));
    graph.addEdge("v1", "v2", Data(1.2));
    graph.addEdge("v1", "v2", Data(1.2));
    graph.addEdge("v1", "v3", Data(2.2));
    graph.addEdge("v2", "v3", Data(3.2));
    graph.addEdge("v2", "v4", Data(4.2));
    graph.addEdge("v3", "v4", Data(5.2));

    cout << graph << endl;    

    graph.removeEdge(1, 2);
    cout << graph << endl;    
    graph.removeVertex(1);
    cout << graph << endl;    

	/*
	//Graph *pG = new Graph();
	//g.printNodes();
	//string str = g.toString();
	//cout << "G2:" << str << endl;
	//cout << g.toString() << endl;
	//Vertex v = Vertex();
	//Vertex v2 = Vertex(1);
	//Vertex v3 = Vertex(2, "abcd"); 
	
	Graph g = Graph();
	g.addVertex();
	g.addVertex("");
	g.addVertex("abc");
	g.addVertex("emr");
	g.addEdge("abc", "emr");
	cout << "G: " << g << endl;
	
	const int x = 5;
	vector<int> list;
	vector<int>::iterator Iter;
	list.push_back(2);
	list.push_back(5);
	
	cout << "myVector =";
	for (Iter = list.begin(); Iter != list.end(); Iter++ )
	{
	     cout << " " << *Iter;
	}
	cout << endl;

	typedef hash_map<uInt,Vertex> Map; // hash_map<int,double> Map; // hash_map<char*, int, hash<char*>, eqstr> Map;
	typedef Map::value_type value_type;
	typedef Map::iterator iterator;	
	Map map;  
	//map.insert(value_type(8, 8.888));
	//map.insert(value_type(1, 1.111));
	//map.insert(value_type(12, 12.12));
	//map.insert(value_type(3, 3.3));
	//map.insert(value_type(122, 122.122));
	map.insert(value_type(5, Vertex(5, "emr")));
	map.insert(value_type(v3.getId(), v3));
  
	std::cout << "\nSize is: " << map.size();
	std::cout << "\nElements are:";
	for (iterator it = map.begin(); it != map.end(); ++it)
	{
		std::cout << "\n\tKey = " << it->first << "   Value = " << it->second;
	}
	*/
	return;
}
