/*
    GUILD (Genes Underlying Inheritance Linked Disorders) implements several 
    graph based algorithms for scoring relevance of a node in the network in 
    terms of a phenotype using known associations in the node's neighborhood 
    for that phenotype. GUILD has been applied to the prioritization of genes 
    for several human disorders. 2011 - Emre Guney (Unviersitat Pompeu Fabra)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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

