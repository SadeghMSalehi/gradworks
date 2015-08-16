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
#include "iostream"
#include "fstream"
#include "algorithm"
#include <deque>

using namespace std;

StringHelper _ss;


vector<string> readWords(string file) {
    ifstream fi(file);
    vector<string> words;
    
    string w;
    while (fi >> w) {
        words.push_back(w);
    }
    return words;
}

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



struct TrieNode {
    char ch;
    bool hasWord;
    std::vector<TrieNode*> nodes;
    
    TrieNode(char c) {
        ch = c;
        nodes.resize(26);
    }
    
    void markWord() {
        hasWord = true;
    }
    
    int char2int(char s) {
        return int(s - 'a');
    }
    
    TrieNode* getNext(char s) {
        return nodes[char2int(s)];
    }
    
    TrieNode* addNode(char s) {
        TrieNode* node = new TrieNode(s);
        nodes[char2int(s)] = node;
        return node;
    }
};


struct Trie {
    TrieNode* root;
    
    Trie() {
        root = new TrieNode('0');
    }
    
    bool insert(string s) {
        TrieNode* iter = root;
        int j = 0;
        while (iter->getNext(s[j]) != NULL && j < s.size()) {
            iter = iter->getNext(s[j]);
            j++;
        }
        if (j == s.size()) {
            iter->hasWord = true;
        }
        while (j < s.size()) {
            iter = iter->addNode(s[j]);
            j++;
        }
        iter->markWord();
        return true;
    }
    
    bool hasKey(string s) {
        if (root == NULL) return false;
        TrieNode* iter = root;
        int j = 0;
        while (j < s.size() && iter->getNext(s[j]) != NULL) {
            iter = iter->getNext(s[j]);
            j++;
        }
        return (j == s.size() && iter != NULL && iter->hasWord);
    }
    
    void print(TrieNode* iter, vector<char>& sofar) {
        sofar.push_back(iter->ch);
        if (iter->hasWord) {
            cout << string(sofar.begin(), sofar.end()) << endl;
        }
        for (int j = 0; j < iter->nodes.size(); j++) {
            if (iter->nodes[j] != NULL) {
                print(iter->nodes[j], sofar);
            }
        }
        sofar.pop_back();
    }
    
    void print() {
        vector<char> collector;
        for (int j = 0; j < root->nodes.size(); j++) {
            if (root->nodes[j] != NULL) {
                print(root->nodes[j], collector);
            }
        }
    }
    
    static void test() {
        vector<string> word2 = readWords("word2.txt");
        vector<string> word3 = readWords("word3.txt");
        Trie trie;
        for (int j = 0; j < word2.size(); j++) {
            trie.insert(word2[j]);
        }
        int cnt = 0, ncnt = 0;
        for (int j = 0; j < word3.size(); j++) {
            if (trie.hasKey(word3[j])) {
                cnt ++;
            } else {
                cout << word3[j] << endl;
                ncnt ++;
            }
        }
        cout << "count: " << cnt << endl;
        
        trie.print();
    }
};



struct Segment {
    float x;
    float y;
    
    Segment(float a, float b): x(a), y(b) {}
    
    bool contain(float p) {
        return x <= p && p <= y;
    }
    
    static bool xcomp(Segment* a, Segment* b) {
        return a->x < b->x;
    }

    static bool ycomp(Segment* a, Segment* b) {
        return a->y > b->y;
    }
};

ostream& operator<<(ostream& os, Segment* s) {
    os << s->x << "-" << s->y;
    return os;
}


struct IntervalTreeNode {
    
    float mid;
    
    vector<Segment*> xlist;
    vector<Segment*> ylist;

    IntervalTreeNode* tl;
    IntervalTreeNode* tr;
    
    IntervalTreeNode(): mid(0), tl(NULL), tr(NULL) {
    }
    
    
    void find(float x, vector<Segment*>& result) {
        if (x == mid) {
            for (int j = 0; j < xlist.size(); j++) {
                result.push_back(xlist[j]);
            }
        } else if (x > mid) {
            for (int j = 0; j < ylist.size(); j++) {
                if (ylist[j]->y >= x) {
                    result.push_back(ylist[j]);
                }
            }
            if (tr) {
                tr->find(x, result);
            }
        } else if (x < mid) {
            for (int j = 0; j < xlist.size(); j++) {
                if (xlist[j]->x <= x) {
                    result.push_back(xlist[j]);
                }
            }
            if (tl) {
                tl->find(x, result);
            }
        }
    }
    
    void buildTree(vector<Segment*>&  segs) {
        vector<float> endpoints;
        for (int j = 0; j < segs.size(); j++) {
            endpoints.push_back(segs[j]->x);
            endpoints.push_back(segs[j]->y);
        }
        
        int v = endpoints.size()/2;
        std::nth_element(endpoints.begin(), endpoints.begin() + v, endpoints.end());
        this->mid = endpoints[v];
        
        vector<Segment*> leftSeg;
        vector<Segment*> rightSeg;
        
        for (int j = 0; j < segs.size(); j++) {
            if (segs[j]->contain(mid)) {
                xlist.push_back(segs[j]);
                ylist.push_back(segs[j]);
                
                sort(xlist.begin(), xlist.end(), &Segment::xcomp);
                sort(ylist.begin(), ylist.end(), &Segment::ycomp);
            } else if (segs[j]->y < mid) {
                leftSeg.push_back(segs[j]);
            } else if (segs[j]->x > mid) {
                rightSeg.push_back(segs[j]);
            }
        }
        
        if (leftSeg.size() > 0) {
            tl = new IntervalTreeNode();
            tl->buildTree(leftSeg);
        }
        if (rightSeg.size() > 0) {
            tr = new IntervalTreeNode();
            tr->buildTree(rightSeg);
        }
    }
    
    
    static void test() {
        vector<Segment*> segs;
        for (int j = 0; j < 10; j++) {
            segs.push_back(new Segment(j, j+3));
        }
        
        IntervalTreeNode* root = new IntervalTreeNode();
        root->buildTree(segs);
        
        for (float j = -.5; j < 14; j+=.5) {
            vector<Segment*> output;
            root->find(j, output);
            _ss.print(output);
        }
    }
};







int main(int argc, char* argv[]) {
	// shortest_path();
	// union_find();
    IntervalTreeNode::test();
}