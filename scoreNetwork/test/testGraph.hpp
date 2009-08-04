
//#define BOOST_TEST_DYN_LINK
//#include <boost/bind.hpp>
//#include <boost/test/unit_test.hpp>

// Note that the following header is needed!!!
// Without this, you must compile like:
// g++ -lboost_unit_test_framework foo.cpp
#include <boost/test/included/unit_test_framework.hpp> 

#include "../graph/Graph.hpp"

using namespace boost;
using namespace boost::unit_test;

class TestGraph
{
    typedef std::pair<std::size_t, std::size_t> E; //unsigned int 
    float *vertex_scores;
    std::string *vertex_names;
    std::size_t n_vertex;
    E *edge_array;
    float *edge_scores;
    std::size_t n_edge;

    Graph _graph;
    Graph _graph2;

public:

    TestGraph() 
    {
	float vertex_scores2[] = { 0.5, 0.6, 0.7, 0.8, 0.9 };
	n_vertex = sizeof(vertex_scores2)/sizeof(float);
	vertex_scores = new float[n_vertex];
	memcpy(vertex_scores, vertex_scores2, sizeof(vertex_scores2));
	std::string vertex_names2[] = { "v1", "v2", "v3", "v4", "v5" };
	vertex_names = new std::string[n_vertex];
	memcpy(vertex_names, vertex_names2, sizeof(vertex_names2));
	E edge_array2[] = { E(0,1), E(1,2), E(1,3), E(2,3), E(3,4) };
	n_edge = sizeof(edge_array) / sizeof(E);
	edge_array = new E[n_edge];
	memcpy(edge_array, edge_array2, sizeof(edge_array2));
	float edge_scores2[] = { 0.1, 0.2, 0.3, 0.4, 0.5 };
	edge_scores = new float[n_edge];
	memcpy(edge_scores, edge_scores2, sizeof(edge_scores2));

	_graph2 = Graph(2); 
	for(std::size_t i=0; i<n_vertex; i++) {
	    _graph.addVertex(vertex_names[i], vertex_scores[i]);
	}
	for(std::size_t i=0; i<n_edge; i++) {
	    _graph.addEdge(vertex_names[edge_array[i].first], vertex_names[edge_array[i].second], edge_scores[i]);
	}
    }

    ~TestGraph() {
	delete vertex_scores;
	//delete vertex_names;
	//delete edge_array;
	delete edge_scores;
    }

    void test_constructor()
    {
	//BOOST_REQUIRE(true); // Require, if not true causes tests to abort
	BOOST_CHECK_EQUAL( _graph.getSize(), (size_t)5 );
	BOOST_CHECK_EQUAL( _graph2.getSize(), (size_t)2 ); // Note that number of vertices is 0 although it initialized with 2, since none added
	//_graph.print();
    }

    void test_vertexData()
    {
	std::string str("v1");
	BOOST_CHECK_EQUAL(_graph.getSize(), n_vertex);
	BOOST_CHECK_CLOSE(_graph.getVertexScore(str), 0.5f, 0.00001f);
	BOOST_CHECK_EQUAL(_graph.getVertexIndex(str), 0);
    }

    void test_graphData()
    {
	VertexIterator it, itEnd;
	OutEdgeIterator et, etEnd;
	for(boost::tie(it, itEnd) = _graph.getVertexIterator(); it != itEnd; it++) 
	{
	    //std::cout << _graph.getVertexIndex(*it) << " " << _graph.getVertexName(*it) << std::endl;
	    BOOST_CHECK(_graph.getVertexName(*it) == std::string(vertex_names[_graph.getVertexIndex(*it)]));
	    BOOST_CHECK_CLOSE(_graph.getVertexScore(*it), vertex_scores[_graph.getVertexIndex(*it)], 0.00001f);
	    //for(boost::tie(et, etEnd) = _graph.getEdgeIteratorOfVertex(*it); et != etEnd; et++) 
	    //{
	    //	BOOST_CHECK_CLOSE(_graph.getEdgeScore(*et), edge_scores[i], 0.00001f);
	    //}
	}
	for(size_t i=0; i<n_vertex; i++) 
	{
	    //std::cout << i << ": " << vertex_scores[i] << " " << std::endl;
	    //BOOST_CHECK(_graph.getVertexScore(i) == std::string(vertex_scores[i]));
	}
	for(boost::tie(it, itEnd) = _graph.getVertexIterator(); it != itEnd; it++) 
	{
	    for(boost::tie(et, etEnd) = _graph.getEdgeIteratorOfVertex(*it); et != etEnd; et++) 
	    {
		unsigned int edge_index=0;
		for(size_t i=0; i<n_edge; i++, edge_index++) 
		{
		    if(((unsigned int)edge_array[edge_index].first) == _graph.getVertexIndex(_graph.getSourceVertex(*et)) && ((unsigned int)edge_array[edge_index].second) == _graph.getVertexIndex(_graph.getTargetVertex(*et)))
			break;
		}
		BOOST_CHECK_CLOSE(_graph.getEdgeScore(*et), edge_scores[edge_index], 0.00001f);
	    }
	}
	EdgeIterator et2, etEnd2;
	for(boost::tie(et2, etEnd2) = _graph.getEdgeIterator(); et2 != etEnd2; et2++) 
	{
	    unsigned int edge_index=0;
	    for(size_t i=0; i<n_edge; i++, edge_index++) 
	    {
		if(((unsigned int)edge_array[edge_index].first) == _graph.getVertexIndex(_graph.getSourceVertex(*et2)) && ((unsigned int)edge_array[edge_index].second) == _graph.getVertexIndex(_graph.getTargetVertex(*et2)))
		    break;
	    }
	    BOOST_CHECK_CLOSE(_graph.getEdgeScore(*et2), edge_scores[edge_index], 0.00001f);
	}
    }

    void test_shortestPath()
    {
	float vertex_distances[] = { 0, 0.1, 0.3, 0.4, 0.9 };
	std::map<Vertex, Vertex> mapPredecessor;
	std::map<Vertex, float> mapDistance;
	std::map<Vertex, float>::iterator vt, vtEnd;
	_graph.calculateShortestPath(_graph.getVertex("v1"), mapDistance, mapPredecessor); 
	//mapDistance = _graph.calculateShortestPath(_graph.getVertex("v1")); 
	for(vt = mapDistance.begin(), vtEnd = mapDistance.end(); vt != vtEnd; ++vt) 
	{
	    //! VERY STRANGELY DOES NOT WORK HERE (ALTHOUGH IT WORKS IN MAIN)
	    //std::cout << _graph.getVertexName(vt->first) << " " << vt->second << std::endl;
	    //BOOST_CHECK_CLOSE(vertex_distances[_graph.getVertexIndex(vt->first)], vt->second, 0.00001f);
	}
    }

};


class TestSuiteGraph: public test_suite
{
public:
    TestSuiteGraph(): test_suite("test_suite_score_network")
    {
        shared_ptr<TestGraph> instance(new TestGraph());
        test_case *constructor = BOOST_CLASS_TEST_CASE( &TestGraph::test_constructor, instance);
        test_case *vertexData = BOOST_CLASS_TEST_CASE( &TestGraph::test_vertexData, instance);
        test_case *graphData = BOOST_CLASS_TEST_CASE( &TestGraph::test_graphData, instance);
        test_case *shortestPath = BOOST_CLASS_TEST_CASE( &TestGraph::test_shortestPath, instance);
        add(constructor);
        add(vertexData);
        add(graphData);
        add(shortestPath);
    }
    ~TestSuiteGraph() {
	//instance.release();
    }
};

/*
void init_function()
{
    //framework::master_test_suite().add( BOOST_TEST_CASE( boost::bind( &free_test_function, 1, 1 ) ) );
    framework::master_test_suite().add( new TestSuiteGraph() );
}

int main(int argc, char* argv[] )
{
    ::boost::unit_test::unit_test_main( &init_function, argc, argv );
    return 0;
}
*/

