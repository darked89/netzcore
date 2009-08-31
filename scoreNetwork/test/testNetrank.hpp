
#include <boost/test/included/unit_test_framework.hpp> 

#include "../Netrank.hpp"

using namespace boost;
using namespace boost::unit_test;

class TestNetrank
{
    Netrank _sNetwork;

public:

    TestNetrank() {
	_sNetwork = Netrank("../../data/toy_data/test_proteins_small.sif", "../../data/toy_data/test_interactions_small.sif", "../../test/data/out_toy_testNetrank.sif");
    }

    ~TestNetrank() {
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

	VertexIterator it, itEnd;
	for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	{
	    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], 0.00001f);
	    //BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScoreInitial(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], 0.00001f);
	}
    }

    void test_calculatedScores()
    {	
	float vertex_scores[] = { 0.03+0.85*0.5/3, 0.03+0.85*0.5/2, 0.03+0.85*0.5/3, 0.03+0.85*(0.5/3+0.5/2), 0.03 };

	_sNetwork.run(1);
	VertexIterator it, itEnd;
	for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	{
	    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], 0.00001f);
	}
    }

};

class TestSuiteNetrank: public test_suite
{
public:
    TestSuiteNetrank(): test_suite("test_suite_score_network")
    {
        shared_ptr<TestNetrank> instance(new TestNetrank());
        test_case *loadData = BOOST_CLASS_TEST_CASE( &TestNetrank::test_loadData, instance);
        test_case *initialScores = BOOST_CLASS_TEST_CASE( &TestNetrank::test_initialScores, instance);
        test_case *calculatedScores = BOOST_CLASS_TEST_CASE( &TestNetrank::test_calculatedScores, instance);
        add(loadData);
        add(initialScores);
        add(calculatedScores);
    }
    ~TestSuiteNetrank() {
	//instance.release();
    }
};

