#include "worker.hpp"
#include <fstream>
#include <sstream>

Worker::Worker(int k, unsigned int seed): k(k), n(0), generator(seed), distribution(0.0, 1.0) {
	srand(seed);
	samples.reserve(k);
	
	globalCnt.resize(6,0);
};


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
	for (int i = 1; i < 5; i++){
		prob *= ((curSampleNum - i) / (n - i));
		count[i-1] = 1.0/prob;
	}
	
	// u = src, v = dst
	
	
	
	
	std::set<VID> &Nu = nodeToNeighbors[src];
	std::set<VID> &Nv = nodeToNeighbors[dst];
	
	std::set<VID> C = std::set<VID>(); // set of nodes that form triangles on u and v	
	
	long double countSums[6] = {0,0,0,0,0,0};
	
	// count triangles and get set of triangle nodes
	for (std::set<VID>::iterator NuIterator = Nu.begin(); NuIterator != Nu.end(); NuIterator++) {
		VID neighbor_of_u = *NuIterator; // get neighbor
		
		if (Nv.find(neighbor_of_u) != Nv.end()){
			C.insert(neighbor_of_u); // put neighbor in intersection
		}
	}
	countSums[0] = C.size() * count[0]; // update triangle count
	
	// count three-paths starting at u, and get set of two-paths starting at v  
	long double cur_sum = 0.0; // 3 path
	long double cur_sum_2 = 0.0; // tailed triangle
	long double cur_sum_3 = 0.0; // chordal 4-cycle
	long double cur_sum_4 = 0.0; // 4 cycle
	for (std::set<VID>::iterator NvIterator = Nv.begin(); NvIterator != Nv.end(); NvIterator++) {
		VID neighbor_of_v = *NvIterator; // get neighbor
		std::set<VID> &Nnov = nodeToNeighbors[neighbor_of_v]; // neigborhood of neighbor of v
		for (std::set<VID>::iterator NnovIterator = Nnov.begin(); NnovIterator != Nnov.end(); NnovIterator++) {
			VID neighbor_of_nov = *NnovIterator; // get neighbor of neighbor
			if (neighbor_of_nov != src && neighbor_of_nov != dst){
				cur_sum += 1;
				if (Nv.find(neighbor_of_nov) != Nv.end()){
					cur_sum_2 += 1;
					if (Nu.find(neighbor_of_nov) != Nu.end()){
						cur_sum_3 += 1;
					}
				}
			}
		}
	}
	
	// count three-paths starting at v, and get set of two-paths starting at u
	for (std::set<VID>::iterator NuIterator = Nu.begin(); NuIterator != Nu.end(); NuIterator++) {
		VID neighbor_of_u = *NuIterator; // get neighbor
		std::set<VID> &Nnou = nodeToNeighbors[neighbor_of_u]; // neigborhood of neighbor of v
		for (std::set<VID>::iterator NnouIterator = Nnou.begin(); NnouIterator != Nnou.end(); NnouIterator++) {
			VID neighbor_of_nou = *NnouIterator; // get neighbor of neighbor
			if (neighbor_of_nou != src && neighbor_of_nou != dst){
				cur_sum += 1;
				bool nv_check = false;
				if (Nv.find(neighbor_of_nou) != Nv.end()) {
					nv_check = true;
					cur_sum_4 += 1;
				}
				
				if (Nu.find(neighbor_of_nou) != Nu.end()) {
					cur_sum_2 += 1;
										
					if (nv_check) {
						cur_sum_3 += 1;
					}
				}				
			}
		}
	}
	
	// three-paths where (u,v) in middle = (Nu.size()-C.size())*(Nv.size()-C.size()) + (Nu.size()-C.size())*C.size() + (Nv.size()-C.size())*C.size())
	
	
	// count all three-paths including those where the edge (u,v) is in the middle
	countSums[1] = (cur_sum + Nu.size()*Nv.size() - C.size()) * count[0];
	
	
	// double counted triangles on u and v earlier loops
	cur_sum_2 /= 2;
	
	// count tailed triangles where (u,v) is part of the triangle
	for (std::set<VID>::iterator CIterator = C.begin(); CIterator != C.end(); CIterator++) {
		VID tri_on_uv = *CIterator;
		cur_sum_2 += nodeToNeighbors[tri_on_uv].size();
		
	}
	cur_sum_2 -= 2*C.size(); // accounting for u and v counted in each neighborhood of tri_on_uv
	cur_sum_2 += C.size() * (Nu.size() + Nv.size() - 2);
	
	countSums[2] = cur_sum_2*count[1]; // tailed triangles
	
	
	// count chordal 4 cycles
	
	cur_sum_3 += (C.size())*(C.size() -1)/2;
	
	countSums[3] = cur_sum_4*count[1];
	countSums[4] = cur_sum_3*count[2];
	


	// count 4-cliques
	cur_sum = 0.0;
	//for (std::set<VID>::iterator CIterator = C.begin(); CIterator != C.end(); CIterator++) {
	// std::set<VID>::iterator CIterator = C.begin();

	while ( !C.empty() ){
		VID tri_on_uv = *C.begin();

		std::set<VID> &Ntouv = nodeToNeighbors[tri_on_uv]; // neigborhood of tri_on_uv
		for (std::set<VID>::iterator NtouvIterator = Ntouv.begin(); NtouvIterator != Ntouv.end(); NtouvIterator++) {
			VID neighbor_of_touv = *NtouvIterator; 
			if (C.find(neighbor_of_touv) != C.end()){
				cur_sum += 1;
			}
		}
		C.erase( C.begin() ); // all cliques on u,v, and tri_on_uv counted, so we dont need it anymore
	}
	countSums[5] = cur_sum*count[3];
	
	// increment counts for incoming edge
	for (int i = 0;i < 6;i++){
		if(countSums[i] > 0) {
			
			//cout << globalCnt[i];
			//cout << " update " <<  i << endl;
			globalCnt[i] += countSums[i];
		}
	}
	
	
	
	//cout << " worker reached here **************************" << endl;
	
	
	
	
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

std::vector<long double> Worker::getGlobalCnt()
{
	return globalCnt;
}

