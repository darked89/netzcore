
CC=g++
CFLAGS=-I $(HOME)/lib/boost_1_38_0/ -O3 -g #-Wall
LDFLAGS=-c

SRCGRAPH=graph/Graph.cpp 
SRCSNETWORK=ScoreNetwork.cpp
SRCNETSCORE=Netscore.cpp
#SRCNETSCORE2=Netscore2.cpp
SRCNETSHORT=Netshort.cpp
SRCNETRANK=Netrank.cpp
SRCNETZCORE=Netzcore.cpp

HEADERGRAPHEXTRA=graph/page_rank_with_initial_ranks.hpp graph/clustering_coefficient.hpp
OBJGRAPH=$(patsubst %.cpp, %.o, $(SRCGRAPH))
OBJSNETWORK=$(patsubst %.cpp, %.o, $(SRCSNETWORK))
OBJNETSCORE=$(patsubst %.cpp, %.o, $(SRCNETSCORE))
OBJNETSHORT=$(patsubst %.cpp, %.o, $(SRCNETSHORT))
OBJNETRANK=$(patsubst %.cpp, %.o, $(SRCNETRANK))
#OBJNETSCORE2=$(patsubst %.cpp, %.o, $(SRCNETSCORE2))
OBJNETZCORE=$(patsubst %.cpp, %.o, $(SRCNETZCORE))

OBJ=$(OBJGRAPH) $(OBJSNETWORK) $(OBJNETSCORE) $(OBJNETSHORT) $(OBJNETRANK) $(OBJNETZCORE) #$(OBJNETSCORE2)
SRCMAIN=main.cpp
SRCMAINNETSCORE=main_netscore.cpp

HEADERTEST=$(wildcard test/test*.hpp)
SRCTEST=main_test.cpp
#SRCTEST=$(wildcard test*.cpp) 
#OBJTEST=$(patsubst %.cpp, %.o, $(SRCTEST))

all: scoreN testN netscore 

scoreN: $(OBJGRAPH) $(OBJSNETWORK) $(OBJNETSCORE) $(OBJNETSHORT) $(OBJNETRANK) $(OBJNETZCORE) $(SRCMAIN)
	$(CC) $(CFLAGS) $(OBJGRAPH) $(OBJSNETWORK) $(OBJNETSCORE) $(OBJNETSHORT) $(OBJNETRANK) $(OBJNETZCORE) $(SRCMAIN) -o scoreN

netscore: $(OBJNETSCORE) $(SRCMAINNETSCORE)
	$(CC) $(CFLAGS) $(OBJGRAPH) $(OBJSNETWORK) $(OBJNETSCORE) $(SRCMAINNETSCORE) -o netscore

testN: $(OBJ) $(HEADERTEST)
	$(CC) $(CFLAGS) $(OBJ) $(SRCTEST) -o testN

$(OBJNETZCORE): $(SRCNETZCORE) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRCNETZCORE)

$(OBJNETRANK): $(SRCNETRANK) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRCNETRANK)

$(OBJNETSHORT): $(SRCNETSHORT) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRCNETSHORT)

$(OBJNETSCORE): $(SRCNETSCORE) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRCNETSCORE)

#$(OBJNETSCORE2): $(SRCNETSCORE2) 
#	$(CC) $(CFLAGS) $(LDFLAGS) $(SRCNETSCORE2)

$(OBJSNETWORK): $(SRCSNETWORK)
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRCSNETWORK)

$(OBJGRAPH): $(SRCGRAPH) $(HEADERGRAPHEXTRA)
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRCGRAPH) -o $(OBJGRAPH)

$(HEADERGRAPHHEXTRA):
	touch Graph.hpp

clean:
	rm -f *.o
	rm -f graph/*.o

