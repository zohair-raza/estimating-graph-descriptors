#ifndef _WORKER_HPP_
#define _WORKER_HPP_

#include "source.hpp"
#include <unordered_map>
#include <set>
#include <stdlib.h>     /* srand, rand */
#include <random>
#include <vector>
#include "time.h"

class Worker {
_PRIVATE:

	const int k; // maximum number of edges that can be stored
	std::vector<Edge> samples; // sampled edges
	long n; // maximum number of edges sent to so far with doStore=1
	std::vector<long double> globalCnt; // global sub-graphs count
	//std::unordered_map<VID, int> nodeToDeg; // degree counts
	//std::unordered_map<VID, float> nodeToCnt[8]; // local sub-graphs count
	std::unordered_map<VID, std::set<VID>> nodeToNeighbors; // sampled graph
	std::default_random_engine generator; // random real number generator
	std::uniform_real_distribution<double> distribution;
public:
	Worker(int k, unsigned int seed);
	void updateCnt(const Edge &iEdge);
	int deleteEdge();
	
	long double count[4]; // count[i] is the count for a graphlet with i+3 edges
	
	void processEdge(const Edge &iEdge);
	std::vector<long double> getGlobalCnt();
	//std::unordered_map<VID, float> & getLocalCnt(int i);
};



#endif // #ifndef _WORKER_HPP_
