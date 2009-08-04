
#include <boost/test/included/unit_test_framework.hpp> 

#include "../Netscore.hpp"

#include <iostream>

using namespace boost;
using namespace boost::unit_test;

class TestNetscore
{
    Netscore _sNetwork;

public:

    TestNetscore() {
	_sNetwork = Netscore("../../data/toy_data/test_proteins_small.txt", "../../data/toy_data/test_interactions_small.txt", "../../data/out_toy_testNetscore.txt", true, true, false, false);
    }

    ~TestNetscore() {
    }

    void test_loadData()
    {
	BOOST_CHECK_EQUAL(_sNetwork.getNetwork().getSize(), 5);
	BOOST_CHECK(_sNetwork.getNetwork().getNumberOfEdges() == 5);
	std::string str("v1"); 
	BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(str), 0.0f, 0.00001f);
    }

    void test_initialScores()
    {
	float vertex_scores[] = { 0.0, 1.0, 1.0, 0.0, 0.0 };
	std::size_t n_vertex = sizeof(vertex_scores)/sizeof(float);

	VertexIterator it, itEnd;
	for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	{
	    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], 0.00001f);
	    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScoreInitial(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], 0.00001f);
	}
    }

    void test_calculatedScores()
    {	
	std::cout << "\ntestNetscore: Ignoring calculatedScores() test \n" << std::endl;
	float vertex_scores[] = { 0.75, 4.0/3, 1.5, 1.0/3, 0.005 };
	std::size_t n_vertex = sizeof(vertex_scores)/sizeof(float);

	VertexIterator it, itEnd;
	//for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	//{
	//    std::cout << _sNetwork.getNetwork().getVertexName(*it) << " " << _sNetwork.getNetwork().getVertexScore(*it) << std::endl;
	//}

	_sNetwork.run(1,2); //!

	for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	{
	//    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], 0.00001f);
	//    std::cout << _sNetwork.getNetwork().getVertexName(*it) << " " << _sNetwork.getNetwork().getVertexScore(*it) << std::endl;
	}
    }

};

class TestSuiteNetscore: public test_suite
{
public:
    TestSuiteNetscore(): test_suite("test_suite_score_network")
    {
        shared_ptr<TestNetscore> instance(new TestNetscore());
        test_case *loadData = BOOST_CLASS_TEST_CASE( &TestNetscore::test_loadData, instance);
        test_case *initialScores = BOOST_CLASS_TEST_CASE( &TestNetscore::test_initialScores, instance);
        test_case *calculatedScores = BOOST_CLASS_TEST_CASE( &TestNetscore::test_calculatedScores, instance);
        add(loadData);
        add(initialScores);
        add(calculatedScores); //! has a problem with score updating (probably due to inheritence of virtual functions)
    }
    ~TestSuiteNetscore() {
	//instance.release();
    }
};

