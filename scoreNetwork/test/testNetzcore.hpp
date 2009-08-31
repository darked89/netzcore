
#include <boost/test/included/unit_test_framework.hpp> 

#include "../Netzcore.hpp"

#include <iostream>

#include <list>

using namespace boost;
using namespace boost::unit_test;


/*
template <class InputIterator>
std::pair<float, float> calculateMeanAndSigma(InputIterator first, InputIterator beyond)
{
    float sum = 0.0, squareSum = 0.0, mean=0.0;
    unsigned int count = 0;
    while(first != beyond)
    {
	sum += *first;
	squareSum += sq(*first);
	count++;
	++first;
    }
    mean = sum/count;
    return std::make_pair(mean, squareSum/count-sq(mean));
}
*/

class TestNetzcore
{
    Netzcore _sNetwork;

public:

    TestNetzcore() {
	_sNetwork = Netzcore("../../data/toy_data/test_proteins_small.sif", "../../data/toy_data/test_interactions_small.sif", "../../data/out_toy_testNetzcore.sif", "../../test/data/sampled_graphs/sampled_graph.sif.", 4);
    }

    ~TestNetzcore() {
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

    void test_calculateMeanAndSigma() 
    {
	float vertex_scores[] = { 0.0, 1.0, 0.0, 1.0 };
	std::list<float> scores(vertex_scores, vertex_scores+sizeof(vertex_scores)/sizeof(float));
	std::pair<float, float> resultPair(0.5, 0.5);
	//std::pair<float, float> tempPair = calculateMeanAndSigma(vertex_scores, vertex_scores+sizeof(vertex_scores)/sizeof(float));
	//std::pair<float, float> tempPair = calculateMeanAndSigma(scores.begin(), scores.end());
	//BOOST_CHECK_CLOSE(resultPair.first, tempPair.first, 0.00001f);
	//BOOST_CHECK_CLOSE(resultPair.second, tempPair.second, 0.00001f);
    }

    void test_calculatedScores()
    {	
	std::cout << "\ntestNetzcore: Ignoring calculatedScores() test \n" << std::endl;
	//float vertex_scores[] = { 0.75, 4.0/3, 1.5, 1.0/3, 0.005 };

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

class TestSuiteNetzcore: public test_suite
{
public:
    TestSuiteNetzcore(): test_suite("test_suite_score_network")
    {
        shared_ptr<TestNetzcore> instance(new TestNetzcore());
        test_case *loadData = BOOST_CLASS_TEST_CASE( &TestNetzcore::test_loadData, instance);
        test_case *initialScores = BOOST_CLASS_TEST_CASE( &TestNetzcore::test_initialScores, instance);
        test_case *calculateMeanAndSigma = BOOST_CLASS_TEST_CASE( &TestNetzcore::test_calculateMeanAndSigma, instance);
        test_case *calculatedScores = BOOST_CLASS_TEST_CASE( &TestNetzcore::test_calculatedScores, instance);
        add(loadData);
        add(initialScores);
	add(calculateMeanAndSigma);
        add(calculatedScores); //! has a problem with score updating (probably due to inheritence of virtual functions)
    }
    ~TestSuiteNetzcore() {
	//instance.release();
    }
};

