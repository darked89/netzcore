#define BOOST_TEST_MODULE Graph test // defines main
#include <boost/test/included/unit_test.hpp>
#include "Graph.hpp"  

//typedef boost::unordered_map<std::string, float> StrToFloat;

typedef std::pair<std::size_t, std::size_t> E; //unsigned int 
float vertex_scores[] = { 0.5, 0.6, 0.7, 0.8, 0.9 };
std::string vertex_names[] = { "v1", "v2", "v3", "v4", "v5" };
std::size_t n_vertex = sizeof(vertex_scores)/sizeof(float);
E edge_array[] = { E(0,1), E(1,2), E(1,3), E(2,3), E(3,4) };
float edge_scores[] = { 0.1, 0.2, 0.3, 0.4, 0.5 };
std::size_t n_edge = sizeof(edge_array) / sizeof(E);

BOOST_AUTO_TEST_CASE( constructors_test )
{
    Graph myObject;
    BOOST_CHECK_EQUAL( myObject.getSize(), (size_t)0 );
    Graph myObject2(2);
    BOOST_CHECK_EQUAL( myObject.getSize(), (size_t)0 ); // Note that number of vertices is 0 although it initialized with 2, since none added
    //BOOST_CHECK( true );
    //BOOST_CHECK_EQUAL(sizeof(size_t), sizeof(unsigned int));
    //BOOST_CHECK_EQUAL(sizeof(size_t), sizeof(vertices_size_type));
}

struct GraphTestFixture
{
    Graph myObject;
    //std::string str1, str2;

    GraphTestFixture()
    {
	for(std::size_t i=0; i<n_vertex; i++) {
	    myObject.addVertex(vertex_names[i], vertex_scores[i]);
	}
	for(std::size_t i=0; i<n_edge; i++) {
	    myObject.addEdge(vertex_names[edge_array[i].first], vertex_names[edge_array[i].second], edge_scores[i]);
	}
	//myObject.print();
    }
};

BOOST_FIXTURE_TEST_CASE( TestCase1, GraphTestFixture )
{
    std::string str("v1");
    BOOST_CHECK_EQUAL(myObject.getSize(), n_vertex);
    BOOST_CHECK_CLOSE(myObject.getVertexScore(str), 0.5f, 0.00001f);
    BOOST_CHECK_EQUAL(myObject.getVertexIndex(str), 0);
}  

BOOST_FIXTURE_TEST_CASE( TestCase2, GraphTestFixture )
{
    VertexIterator it, itEnd;
    for(boost::tie(it, itEnd) = myObject.getVertexIterator(); it != itEnd; it++) 
    {
	for(size_t i=0; i<n_vertex; i++) 
	{
	    if(vertex_scores[i]==myObject.getVertexScore(vertex_names[i]))
		BOOST_CHECK(myObject.getVertexName(*it) == std::string(vertex_names[myObject.getVertexIndex(*it)]));
	}
    	BOOST_CHECK_CLOSE(myObject.getVertexScore(*it), vertex_scores[myObject.getVertexIndex(*it)], 0.00001f);
    }
}

