
#include <boost/test/included/unit_test_framework.hpp> 
#include <boost/test/included/unit_test.hpp>

#include "test/testGraph.hpp"
#include "test/testScoreNetwork.hpp"
#include "test/testNetshort.hpp"
#include "test/testNetscore.hpp"
#include "test/testNetrank.hpp"
//#include "test/testNetscore2.hpp"
#include "test/testNetzcore.hpp"

using namespace boost;
using namespace boost::unit_test;

test_suite* init_unit_test_suite(int argc, char** argv)
{
    test_suite* suite_graph(BOOST_TEST_SUITE("Graph Suite"));
    suite_graph->add(new TestSuiteGraph());
    test_suite* suite_score_network(BOOST_TEST_SUITE("Score Network Suite"));
    suite_score_network->add(new TestSuiteScoreNetwork());
    test_suite* suite_netshort(BOOST_TEST_SUITE("Netshort Suite"));
    suite_netshort->add(new TestSuiteNetshort());
    test_suite* suite_netscore(BOOST_TEST_SUITE("Netscore Suite"));
    suite_netscore->add(new TestSuiteNetscore());
    //test_suite* suite_netscore2(BOOST_TEST_SUITE("Netscore2 Suite"));
    //suite_netscore2->add(new TestSuiteNetscore2());
    test_suite* suite_netrank(BOOST_TEST_SUITE("Netrank Suite"));
    suite_netrank->add(new TestSuiteNetrank());

    test_suite* suite_netzcore(BOOST_TEST_SUITE("Netzcore Suite"));
    suite_netzcore->add(new TestSuiteNetzcore());

    framework::master_test_suite().add(suite_graph);
    framework::master_test_suite().add(suite_score_network);
    framework::master_test_suite().add(suite_netshort);
    framework::master_test_suite().add(suite_netscore);
    framework::master_test_suite().add(suite_netrank);
    //framework::master_test_suite().add(suite_netscore2);
    framework::master_test_suite().add(suite_netzcore);
    
    return 0;
}

