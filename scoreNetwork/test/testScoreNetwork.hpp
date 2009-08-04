
#include <boost/test/included/unit_test_framework.hpp> 

#include "../ScoreNetwork.hpp"

using namespace boost;
using namespace boost::unit_test;

class TestScoreNetwork
{
    ScoreNetwork _sNetwork;

public:

    TestScoreNetwork() {
	//_sNetwork.loadNodes("../../data/toy_data/test_proteins_small.txt");
	//_sNetwork.loadEdges("../../data/toy_data/test_interactions_small.txt");
	_sNetwork = ScoreNetwork("../../data/toy_data/test_proteins_small.txt", "../../data/toy_data/test_interactions_small.txt", "../../data/out_toy_testScoreNetwork.txt");
    }

    ~TestScoreNetwork() {
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
	}
    }

    void test_scaledScores()
    {	
	float vertex_scores[] = { 0.0, 1.0, 1.0, 0.0, 0.0 };
	std::size_t n_vertex = sizeof(vertex_scores)/sizeof(float);

	_sNetwork.scaleNodeScoresWrapper();
	VertexIterator it, itEnd;
	for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	{
	    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], 0.00001f);
	}
    }

};

class TestSuiteScoreNetwork: public test_suite
{
public:
    TestSuiteScoreNetwork(): test_suite("test_suite_score_network")
    {
        shared_ptr<TestScoreNetwork> instance(new TestScoreNetwork());
        test_case *loadData = BOOST_CLASS_TEST_CASE( &TestScoreNetwork::test_loadData, instance);
        test_case *initialScores = BOOST_CLASS_TEST_CASE( &TestScoreNetwork::test_initialScores, instance);
        test_case *scaledScores = BOOST_CLASS_TEST_CASE( &TestScoreNetwork::test_scaledScores, instance);
        add(loadData);
        add(initialScores);
        add(scaledScores);
    }
    ~TestSuiteScoreNetwork() {
	//instance.release();
    }
};

