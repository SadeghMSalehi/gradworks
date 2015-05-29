//
//  glmain.cpp
//  interviews
//
//  Created by Joohwi Lee on 3/1/14.
//
//

#include "glmain.h"
#include "height_union.h"
#include "graph.h"
#include "stringhelper.h"
#include <deque>

using namespace std;

StringHelper _ss;

bool node_comp(std::vector<int> i, std::vector<int> j) {
	return i[1] < j[1];
}

void shortest_path() {
	typedef vector<int> Edges;
	typedef vector<int> Node;
	
	vector<Node> V = _ss.lstlst2int("1,0 2,0 3,0 4,0 5,0 6,0", " ", ",");
	vector<Edges> E = _ss.lstlst2int("2,3,6/1,3,4/1,2,4,6/2,3,5/4,6/1,3,5", "/", ",");

	vector<Edges> W = _ss.lstlst2int("7,9,14 7,10,15 9,10,11,2 15,11,6 6,9 14,2,9", " ", ",");
	
	deque<Node> q;
	int inf = 1e5;
	for (int j = 0; j < V.size(); j++) {
		Node n = V[j];
		n[1] = inf;
		q.push_back(n);
	}
	q[0][1] = 0;
	
	while (!q.empty()) {
		// pick a node with the minimum distance
		std::sort(q.begin(), q.end(), node_comp);
		Node v = q.front(); q.pop_front();
		V[v[0]-1][1] = v[1];

		// look up neighbors
		Edges n = E[v[0]-1];
		Edges w = W[v[0]-1];
		
		int dv = v[1];
		
		// process neighbors
		for (int j = 0; j < n.size(); j++) {
			int nj = n[j];
			int wj = w[j];
			
			for (int k = 0; k < q.size(); k++) {
				if (q[k][0] == nj) {
					Node& n = q[k];
					if (dv + wj <= n[1]) {
						n[1] = v[1] + wj;
					}
					q[k] = n;
				}
			};
			
		}
	}
	
	for (int j = 0; j < V.size(); j++) {
		cout << V[j][1] << endl;
	}
}


struct UnionFind {
	vector<int> objs;
	vector<int> objids;
	
	int root(int p) {
		p = objids[p];
		while (p != objids[p]) {
			int q = objids[p];
			objids[p] = objids[objids[p]];
			p = q;
		}
		return p;
	}
	
	bool find(int p, int q) {
		return root(p) == root(q);
	}
	
	void unite(int p, int q) {
		objids[p] = q;
	}
};

void union_find() {
	UnionFind uf;
	uf.objs = _ss.split2int("0 1 2 3 4 5 6 7 8 9", " ");
	uf.objids = _ss.split2int("0 1 2 3 4 5 6 7 8 9", " ");
	

	uf.unite(0,1);
	uf.unite(1,2);
	
	_ss.print(uf.objids);
	cout << uf.find(0,1) << " " << uf.find(1,2) << " " << uf.find(0,3) << endl;
}



int main(int argc, char* argv[]) {
	// shortest_path();
	union_find();
}