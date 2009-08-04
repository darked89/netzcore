#./main -p ../data/toy_data/test_proteins_small.txt -i ../data/toy_data/test_interactions_small.txt -n 3 -e 0.01 -N # ~netscore 
#./main -p ../data/toy_data/test_proteins_small.txt -i ../data/toy_data/test_interactions_small.txt -n 3 -e 0.01 -X -z -g -S 3 # default top preserving -R 101
#./main -p ../data/toy_data/test_proteins_small.txt -i ../data/toy_data/test_interactions_small.txt -n 3 -e 0.01 -X -z -g -S 3 -R 102 # random
./main -p ../data/toy_data/test_proteins_small.txt -i ../data/toy_data/test_interactions_small.txt -n 3 -e 0.01 -X -z -g -S 3 -R 101
