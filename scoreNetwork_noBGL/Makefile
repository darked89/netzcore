
CC=g++
CFLAGS=-O2
LDFLAGS=
GRAPH_SRC=./graph/Graph.cpp ./graph/Vertex.cpp ./graph/Edge.cpp 
GRAPH_OBJ=./Graph.o ./Vertex.o ./Edge.o 
SRC=ScoreNetwork.cpp
OBJ=ScoreNetwork.o

#all: graph scoreNetwork 
all: $(GRAPH_OBJ) $(OBJ)
	$(CC) $(CFLAGS) $(GRAPH_OBJ) $(OBJ) main.cpp -o test

#scoreNetwork: ScoreNetwork.cpp
#	$(CC) -c ScoreNetwork.cpp -o $@

#scoreNetwork: $(SRC)
#	$(CC) $(CFLAGS) -c $(SRC) -o $@

$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c $(SRC)

#graph: $(GRAPH_SRC)
#    	$(CC) -c $(GRAPH_SRC) $(CFLAGS) $(LDFLAGS)

#graph: $(GRAPH_SRC)
#	#$(CC) -c graph/*.cpp 
#	$(CC) -c $(GRAPH_SRC)

$(GRAPH_OBJ): $(GRAPH_SRC)
	$(CC) $(CFLAGS) -c $(GRAPH_SRC)

clean:
	rm -f *.o

