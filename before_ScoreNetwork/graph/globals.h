#ifndef GLOBALS_H_
#define GLOBALS_H_
 
#include <string>
#include <sstream>
#include <iostream>

using namespace std;

//typedef unsigned int uInt;

#define INT_INFINITY 9999999
#define sq(x) (x*x)

typedef struct data_struct { 
	float score;
	float scoreUpdated;
	//float scoreDual;
	float scoreInitial;
	float scoreAccumulatedLast; // in since initial scores is not modified after zScoring used to get initial zScore, in edges for knowing which neighbor give what amount of contribution at the last step 
	//float scoreAccumulatedLastUpdated;  
	data_struct() { score = 0.0; scoreUpdated = 0.0; scoreInitial = 0.0; scoreAccumulatedLast = 0.0; /*scoreAccumulatedLastUpdated=0.0; scoreDual = 0.0;*/ };
	data_struct(float valScore) { score = valScore; scoreUpdated = 0.0; scoreInitial = valScore; /*scoreAccumulatedLast = 0.0; scoreAccumulatedLastUpdated=0.0; scoreDual = 0.0;*/ };
} Data;

/*
template <class T>
inline std::string to_string(const T& t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}
*/

#endif /*GLOBALS_H_*/

