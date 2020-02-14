#include "worker.hpp"
#include <fstream>
#include <sstream>

Worker::Worker(int k, unsigned int seed): k(k), n(0), generator(seed), distribution(0.0, 1.0) {
	srand(seed);
	samples.reserve(k);
};


inline void Worker::incCnt(VID node, double count, int type){
	if (nodeToCnt[type].find(node) == nodeToCnt[type].end()) {
		nodeToCnt[type][node] = count;
		return;
	}
	nodeToCnt[type][node] += count;
}


void Worker::updateCnt(const Edge &iEdge){

	VID src = iEdge.src;
	VID dst = iEdge.dst;

	//cout << " on arrival of " << src << "-" << dst << endl;
	
	// can't find any sub-graphs using no neighbors
	if(nodeToNeighbors.find(src) == nodeToNeighbors.end() && nodeToNeighbors.find(dst) == nodeToNeighbors.end()) {
		return;
	}
	
	// swap so we iterate through lower number of neighbors
    if(nodeToNeighbors[src].size() > nodeToNeighbors[dst].size()) {
        VID temp = dst;
        dst = src;
        src = temp;
    }
	
	long double curSampleNum = k >= n ? n : k;
	
	count[0] = n / curSampleNum;
	count[1] = count[0]*((n - 1)/(curSampleNum - 1));
	
	// u = src, v = dst
	
	std::set<VID> &Nu = nodeToNeighbors[src];
	std::set<VID> &Nv = nodeToNeighbors[dst];
	
	
	long double cur_sum_c = 0.0;
	long double cur_sum_up = 0.0;
	long double cur_sum_vp = 0.0;
	
	//cout << " start" << endl;
	
	
	std::set<VID>::iterator NuIterator = Nu.begin();
	std::set<VID>::iterator NvIterator = Nv.begin();
	
	while (NuIterator != Nu.end() && NvIterator != Nv.end()) {
		VID neighbor_of_u = *NuIterator; // get neighbor
		VID neighbor_of_v = *NvIterator; // get neighbor
		
		if (neighbor_of_u == neighbor_of_v){
			// triangle found
			
			cur_sum_c += 1;
			incCnt(neighbor_of_u, count[1], 0);
			
			NuIterator++;
			NvIterator++;
		}
		else if (neighbor_of_u > neighbor_of_v) {
			NvIterator++;
		} else {
			NuIterator++;
		}
	}
	
	cur_sum_c *= count[1];
	
	for (std::set<VID>::iterator NuIterator = Nu.begin(); NuIterator != Nu.end(); NuIterator++) {
		VID neighbor_of_u = *NuIterator; // get neighbor
		incCnt(neighbor_of_u, count[0], 1);
	}
	cur_sum_vp = Nu.size() * count[0];
	
	for (std::set<VID>::iterator NvIterator = Nv.begin(); NvIterator != Nv.end(); NvIterator++) {
		VID neighbor_of_v = *NvIterator; // get neighbor
		incCnt(neighbor_of_v, count[0], 1);
	}
	
	cur_sum_up = Nv.size() * count[0];
	
	
	if (cur_sum_up > 0) {
		if (nodeToCnt[1].find(src) == nodeToCnt[1].end()) {
			nodeToCnt[1][src] = cur_sum_up;
		} else {
			nodeToCnt[1][src] += cur_sum_up;
		}
	}
	
	if (cur_sum_vp > 0) {
		if (nodeToCnt[1].find(dst) == nodeToCnt[1].end()) {
			nodeToCnt[1][dst] = cur_sum_vp;
		} else {
			nodeToCnt[1][dst] += cur_sum_vp;
		}
	}
	
	
	if(cur_sum_c > 0) {
		if (nodeToCnt[0].find(src) == nodeToCnt[0].end()) {
			nodeToCnt[0][src] = cur_sum_c;
		} else {
			nodeToCnt[0][src] += cur_sum_c;
		}

		if (nodeToCnt[0].find(dst) == nodeToCnt[0].end()) {
			nodeToCnt[0][dst] = cur_sum_c;
		} else {
			nodeToCnt[0][dst] += cur_sum_c;
		}
	}
	
	return;
}

int Worker::deleteEdge() {
	int index = rand() % k;
	Edge removedEdge = samples[index];
	nodeToNeighbors[removedEdge.src].erase(removedEdge.dst);
	nodeToNeighbors[removedEdge.dst].erase(removedEdge.src);
	return index;
}

void Worker::processEdge(const Edge &iEdge){

	VID src = iEdge.src;
	VID dst = iEdge.dst;

	if(src == dst) { //ignore self loop
		return;
	}

	updateCnt(iEdge); //count sub-graphs involved

	bool isSampled = false;
	if(n < k) { // always sample
		isSampled = true;
	}
	else {

		if(distribution(generator) < k / (1.0+n)) {
			isSampled = true;
		}
	}

	if(isSampled) {


		if(n < k) {
			samples.push_back(Edge(iEdge));
		}

		else {
			int index = deleteEdge();
			samples[index] = iEdge;
		}

		if(nodeToNeighbors.find(src)==nodeToNeighbors.end()) {
			nodeToNeighbors[src] = std::set<VID>();
		}
		nodeToNeighbors[src].insert(dst);

		if(nodeToNeighbors.find(dst)==nodeToNeighbors.end()) {
			nodeToNeighbors[dst] = std::set<VID>();
		}
		nodeToNeighbors[dst].insert(src);
	}

	n++;

	return;

}

std::unordered_map<VID, long double> & Worker::getLocalCnt(int i)
{
	return nodeToCnt[i];
}
