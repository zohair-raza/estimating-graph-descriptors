#include "worker.hpp"
#include <fstream>
#include <sstream>

Worker::Worker(int k, unsigned int seed): k(k), n(0), generator(seed), distribution(0.0, 1.0) {
	srand(seed);
	samples.reserve(k);
	terms.resize(3,0);
};

void Worker::updateDeg(const Edge &iEdge){
	VID src = iEdge.src;
	VID dst = iEdge.dst;

	if (nodeToDeg.find(src) == nodeToDeg.end()) {
		nodeToDeg[src] = 1;
	} else {
		nodeToDeg[src] += 1;	
	}

	if (nodeToDeg.find(dst) == nodeToDeg.end()) {
		nodeToDeg[dst] = 1;
	} else {
		nodeToDeg[dst] += 1;	
	}
}

void Worker::updateCnt(const Edge &iEdge){

	VID src = iEdge.src;
	VID dst = iEdge.dst;

	//cout << " on arrival of " << src << "-" << dst << endl;
	
	// can't find any sub-graphs using no other neighbors
	
	if(nodeToNeighbors.find(src) == nodeToNeighbors.end() && nodeToNeighbors.find(dst) == nodeToNeighbors.end()) {
		return;
	}
	
	// swap so we iterate through lower number of neighbors
    if(nodeToNeighbors[src].size() > nodeToNeighbors[dst].size()) {
        VID temp = dst;
        dst = src;
        src = temp;
    }
	
	long double  curSampleNum = k >= n ? n : k;
	long double  prob = ((curSampleNum) / (n)); //(curSampleNum / n * (curSampleNum - 1) / (n - 1));
	count[0] = 1.0/prob;
	for (int i = 1; i < 3; i++){
		prob *= ((curSampleNum - i) / (n - i));
		count[i] = 1.0/prob;
	}
	

	// u = src, v = dst
	
	std::set<VID> &Nu = nodeToNeighbors[src];
	std::set<VID> &Nv = nodeToNeighbors[dst];
	
	long double du = nodeToDeg[src];
	long double dv = nodeToDeg[dst];

	long double termCounts[3] = {0,0,0};

	// count triangles and get set of triangle nodes
	std::set<VID>::iterator NuIterator = Nu.begin();
	std::set<VID>::iterator NvIterator = Nv.begin();
	
	
	while (NuIterator != Nu.end() && NvIterator != Nv.end()) {
		VID neighbor_of_u = *NuIterator; // get neighbor
		VID neighbor_of_v = *NvIterator; // get neighbor
		
		if (neighbor_of_u == neighbor_of_v){
			// triangle found
			long double dw = nodeToDeg[neighbor_of_u];
			termCounts[1] += -1/(du*dv*dw);
			
			NuIterator++;
			NvIterator++;
		}
		else if (neighbor_of_u > neighbor_of_v) {
			NvIterator++;
		} else {
			NuIterator++;
		}
	}
	// count two-paths
	for (std::set<VID>::iterator NuIterator = Nu.begin(); NuIterator != Nu.end(); NuIterator++) {
		VID neighbor_of_u = *NuIterator; // get neighbor
		long double dw = nodeToDeg[neighbor_of_u];
		termCounts[0] += 1/(du*dv*dw*du);
	}
	
	for (std::set<VID>::iterator NvIterator = Nv.begin(); NvIterator != Nv.end(); NvIterator++) {
		VID neighbor_of_v = *NvIterator; // get neighbor
		long double dw = nodeToDeg[neighbor_of_v];
		termCounts[0] += 1/(du*dv*dw*dv);
	}

	// count three-paths starting at v, and get set of two-paths starting at u
	for (std::set<VID>::iterator NuIterator = Nu.begin(); NuIterator != Nu.end(); NuIterator++) {
		VID neighbor_of_u = *NuIterator; // get neighbor
		std::set<VID> &Nnou = nodeToNeighbors[neighbor_of_u]; // neigborhood of neighbor of v
		for (std::set<VID>::iterator NnouIterator = Nnou.begin(); NnouIterator != Nnou.end(); NnouIterator++) {
			VID neighbor_of_nou = *NnouIterator; // get neighbor of neighbor
			if (neighbor_of_nou != src && neighbor_of_nou != dst){
				if (Nv.find(neighbor_of_nou) != Nv.end()) {
					long double dw = nodeToDeg[neighbor_of_u];
					long double dx = nodeToDeg[neighbor_of_nou];
					termCounts[2] += 1/(du*dv*dw*dx);
				}		
			}
		}
	}
	
	for ( int i = 0 ; i < 3 ; i++ ){
		long double addTerm = termCounts[i]*count[i];
		if (!isnan(addTerm)){
			terms[i] += addTerm;
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

void Worker::processEdge1(const Edge &iEdge){

	VID src = iEdge.src;
	VID dst = iEdge.dst;

	if(src == dst) { //ignore self loop
		return;
	}

	//cout << "1: src: " << src << " dst: " << dst << endl;
	

	updateDeg(iEdge);
	return;
}


void Worker::processEdge2(const Edge &iEdge){

	VID src = iEdge.src;
	VID dst = iEdge.dst;

	if(src == dst) { //ignore self loop
		return;
	}

	//cout << "2: src: " << src << " dst: " << dst << endl;
	

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

std::vector<long double> Worker::getTerms()
{
	return terms;
}

