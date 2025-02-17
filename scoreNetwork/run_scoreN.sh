#!/bin/bash

data_dir="../../data/toy_data/"
node_file=$data_dir"test_proteins_small.sif"
network_file=$data_dir"test_interactions_small.sif"
#node_file=$data_dir"test_proteins_middle.sif"
#network_file=$data_dir"test_interactions_middle.sif"
#node_file=$data_dir"test_proteins_larger.sif"
#network_file=$data_dir"test_interactions_larger.sif"
sampling_prefix=$data_dir"sampled_graphs_toy/sampled_graph.sif."

#./scoreN -n $node_file -e ../../data/toy_data/test_interactions_small_with_scores.sif -o test.txt -d $sampling_prefix -x 4 -s z -i 1
#./scoreN -n ../../data/toy_data/test_proteins_middle.sif -e ../../data/toy_data/test_interactions_middle.sif -o test.txt -d ../../toy_data/sampling_graphs_toy_middle/ -x 4 -s z -i 1
#./scoreN -n ../../data/toy_data/test_proteins_small.sif -e ../../data/toy_data/test_interactions_small.sif  -o test/test_netween.txt -s w

# Runs on test files

# Netshort
#./scoreN -n $node_file -e $network_file -o test/test_netshort.txt -s d 

# Netrank
./scoreN -n $node_file -e $network_file -o test/test_netrank.txt -s r 

# Netween
#./scoreN -n $node_file -e $network_file -o test/test_netween.txt -s w 

# Netlink
#./scoreN -n $node_file -e $network_file -o test/test_netlink.txt -t 2 -s l

# Netscore
#./scoreN -n $node_file -e $network_file -o test/test_netscore.txt -d -s s -r 1 -i 1

# Netzcore
#./scoreN -n $node_file -e $network_file -o test/test_netzcore.txt -d $sampling_prefix -x 4 -s z -r 1 -i 1

# Netzscore
#./scoreN -n $node_file -e $network_file -o test/test_netzscore.txt -d $sampling_prefix -x 4 -s h -r 1 -i 1

# Netz1score
#./scoreN -n $node_file -e $network_file -o test/test_netz1score.txt -d $sampling_prefix -x 4 -s 1 -r 1 -i 1

