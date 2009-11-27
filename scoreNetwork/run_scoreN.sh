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

#./scoreN -n ../../data/toy_data/test_proteins_small.sif -e ../../data/toy_data/test_interactions_small_with_scores.sif -o test.txt -d ../../data/toy_data/sampled_graphs_toy/ -x 4 -s z -i 1
#./scoreN -n ../../data/toy_data/test_proteins_middle.sif -e ../../data/toy_data/test_interactions_middle.sif -o test.txt -d ../../data/toy_data/sampled_graphs_toy_middle/ -x 4 -s z -i 1

./scoreN -n ../../data/toy_data/test_proteins_small.sif -e ../../data/toy_data/test_interactions_small.sif -o test_netscore.txt -d ../../data/toy_data/sampled_graphs_toy/ -x 4 -s z -r 1 -i 1
./scoreN -n ../../data/toy_data/test_proteins_small.sif -e ../../data/toy_data/test_interactions_small.sif -o test_netzcore.txt -d ../../data/toy_data/sampled_graphs_toy/ -x 4 -s z -r 1 -i 1
./scoreN -n ../../data/toy_data/test_proteins_small.sif -e ../../data/toy_data/test_interactions_small.sif -o test_netzscore.txt -d ../../data/toy_data/sampled_graphs_toy/ -x 4 -s h -r 1 -i 1
./scoreN -n ../../data/toy_data/test_proteins_small.sif -e ../../data/toy_data/test_interactions_small.sif -o test_netz1score.txt -d ../../data/toy_data/sampled_graphs_toy/ -x 4 -s 1 -r 1 -i 1

