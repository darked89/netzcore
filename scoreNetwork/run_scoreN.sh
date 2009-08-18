#./scoreN -n ../../data/input/node_scores.txt -e ../../data/input/edge_relevance_scores.txt -o ../../data/output/netscore_onlyEdgeRelevance_noNodeInitialScoreAccumulation_i3.txt -s s -i 3 &> ../../data/output/netscore_onlyEdgeRelevance_noNodeInitialScoreAccumulation_i3.err
#./scoreN -n ../../data/input/node_scores.txt -e ../../data/input/edge_relevance_scores.txt -o ../../data/output/netshort_onlyEdgeRelevance_noNodeInitialScoreAccumulation.txt -s d &> ../../data/output/netshort_onlyEdgeRelevance_noNodeInitialScoreAccumulation.err
#./scoreN -n ../../data/input/node_scores.txt -e ../../data/input/edge_relevance_scores.txt -o ../../data/output/netrank_onlyEdgeRelevance_noNodeInitialScoreAccumulation.txt -s r &> ../../data/output/netrank_onlyEdgeRelevance_noNodeInitialScoreAccumulation.err
#./scoreN -n ../../data/input/node_scores.txt -e ../../data/input/edge_relevance_scores.txt -o ../../data/output/netzcore_onlyEdgeRelevance_noNodeInitialScoreAccumulation.txt -s z &> ../../data/output/netzcore_onlyEdgeRelevance_noNodeInitialScoreAccumulation.err
#./scoreN -n ../../data/input/node_scores.txt -e ../../data/input/edge_relevance_scores.txt -o ../../data/output/netzcore_noEdgeRelevance_noNodeInitialScoreAccumulation.txt -s z &> ../../data/output/netzcore_noEdgeRelevance_noNodeInitialScoreAccumulation.err

# Runs on test files
#./scoreN -n ../../data/toy_data/test_proteins_small.txt -e ../../data/toy_data/test_interactions_small.txt -o ../../data/output/test.txt -s d 
#./scoreN -n ../../data/toy_data/test_proteins_small.txt -e ../../data/toy_data/test_interactions_small.txt -o ../../data/output/test.txt -s s -i 3
#./scoreN -n ../../data/toy_data/test_proteins_small.txt -e ../../data/toy_data/test_interactions_small.txt -o ../../data/output/test.txt -s r 
#./scoreN -n ../../data/toy_data/test_proteins_small.txt -e ../../data/toy_data/test_interactions_small.txt -o ../../data/output/test.txt -s z 

