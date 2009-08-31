
#include <boost/test/included/unit_test_framework.hpp> 

#include "../ScoreNetwork.hpp"

using namespace boost;
using namespace boost::unit_test;

class TestScoreNetwork
{
    ScoreNetwork _sNetwork;

public:

    TestScoreNetwork() {
	_sNetwork = ScoreNetwork("../../data/toy_data/test_proteins_small_2.sif", "../../data/toy_data/test_interactions_small.sif", "../../data/test/out_toy_testScoreNetwork.sif");
    }

    ~TestScoreNetwork() {
    }

    void test_loadData()
    {
	BOOST_CHECK_EQUAL(_sNetwork.getNetwork().getSize(), 5);
	BOOST_CHECK(_sNetwork.getNetwork().getNumberOfEdges() == 5);
	std::string str("v1"); 
	BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(str), 0.0f, ScoreNetwork::getPrecision());
    }

    void test_initialScores()
    {
	float vertex_scores[] = { 0.1, 0.3, 0.5, 0.7, 0.9 };
	//std::size_t n_vertex = sizeof(vertex_scores)/sizeof(float);

	VertexIterator it, itEnd;
	for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	{
	    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], ScoreNetwork::getPrecision());
	}
    }

    void test_scaledScores()
    {	
	float vertex_scores[] = { 0.0, 1.0, 1.0, 0.0, 0.0 };

	_sNetwork.scaleNodeScoresWrapper(SCALE_BY_MAX_SCORE);

	VertexIterator it, itEnd;
	for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	{
	    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], ScoreNetwork::getPrecision());
	}
    }

    void test_scaledScores_by_max()
    {	
	float vertex_scores[] = { 0.111111, 0.333333, 0.555555, 0.777777, 1.0 };

	_sNetwork.scaleNodeScoresWrapper(SCALE_BY_MAX_SCORE);

	VertexIterator it, itEnd;
	for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	{
	    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], ScoreNetwork::getPrecision());
	}
    }

    void test_scaledScores_between_zero_one()
    {	
	float vertex_scores[] = { 0.0, 0.25, 0.5, 0.75, 1.0 };

	_sNetwork.scaleNodeScoresWrapper(SCALE_BETWEEN_ZERO_AND_ONE);

	VertexIterator it, itEnd;
	for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	{
	    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], ScoreNetwork::getPrecision());
	}
    }

    void test_scaledScores_between_initial_min_max()
    {	
	float vertex_scores[] = { 0.1, 0.3, 0.5, 0.7, 0.9 };

	_sNetwork.scaleNodeScoresWrapper(SCALE_BETWEEN_INITIAL_MIN_AND_MAX_SCORE);

	VertexIterator it, itEnd;
	for(boost::tie(it, itEnd) = _sNetwork.getNetwork().getVertexIterator(); it != itEnd; it++) 
	{
	    BOOST_CHECK_CLOSE(_sNetwork.getNetwork().getVertexScore(*it), vertex_scores[_sNetwork.getNetwork().getVertexIndex(*it)], ScoreNetwork::getPrecision());
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
        //test_case *scaledScores = BOOST_CLASS_TEST_CASE( &TestScoreNetwork::test_scaledScores, instance);
        test_case *scaledScores_by_max = BOOST_CLASS_TEST_CASE( &TestScoreNetwork::test_scaledScores_by_max, instance);
        test_case *scaledScores_between_zero_one = BOOST_CLASS_TEST_CASE( &TestScoreNetwork::test_scaledScores_between_zero_one, instance);
        test_case *scaledScores_between_initial_min_max = BOOST_CLASS_TEST_CASE( &TestScoreNetwork::test_scaledScores_between_initial_min_max, instance);
        add(loadData);
        add(initialScores);
        //add(scaledScores);
        add(scaledScores_by_max);
        add(scaledScores_between_zero_one);
        add(scaledScores_between_initial_min_max);
    }
    ~TestSuiteScoreNetwork() {
	//instance.release();
    }
};

