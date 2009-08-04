#
#CC=g++
#CFLAGS=-c -g -Wall
#LDFLAGS=
#SOURCES=main.cpp netzcore.cpp graph/Graph.cpp graph/Vertex.cpp graph/Edge.cpp
#OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=main

#all: $(SOURCES) $(EXECUTABLE)

#$(EXECUTABLE): $(OBJECTS) 
#	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

#.cpp.o:
#		$(CC) $(CFLAGS) $< -o $@

#object:
#	g++ -c Edge.cpp Vertex.cpp Graph.cpp

all: netzcore
	#g++ graph/*.cpp *.cpp -o main -Wno-deprecated -g -Wall
	#g++ *.cpp test/*.cpp -o main -g

netzcore: graph.o edge.o vertex.o
	g++ *.o Netzcore.cpp main.cpp -o main -Wno-deprecated -g -Wall -O2

graph.o:
	g++ -c graph/Graph.cpp -Wno-deprecated -g -Wall -O2

vertex.o:
	g++ -c graph/Vertex.cpp -Wno-deprecated -g -Wall -O2

edge.o:
	g++ -c graph/Edge.cpp -Wno-deprecated -g -Wall -O2

clear:
	rm *.o

