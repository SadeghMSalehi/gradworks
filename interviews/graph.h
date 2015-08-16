//
//  graph.h
//  interviews
//
//  Created by Joowhi Lee on 5/24/15.
//
//

#ifndef __interviews__graph__
#define __interviews__graph__

#include <stdio.h>
#include <vector>
#include <algorithm>
#include "stringhelper.h"

struct Node {
    int _id;
    double _val;
    
    Node(int i): _id(i), _val(0) {}
};

struct Edge {
    int _s;
    int _e;
    double _w;
    
    Edge(int s, int e, double w): _s(s), _e(e), _w(w) {}
    
    static bool weight_comp(Edge& e1, Edge& e2) {
        return e1._w < e2._w;
    }
};

typedef std::vector<Node> NodeList;
typedef std::vector<Edge> EdgeList;
typedef std::vector<EdgeList> AdjacencyList;

struct Graph {
    NodeList _nodes;
    AdjacencyList _edges;
    
    void addNode(std::string nodes) {
        StringHelper ss;
        std::vector<int> n = ss.convert2int(ss.split(nodes, " "));
        for (int j = 0; j < n.size(); j++) {
            addNode(Node(n[j]));
        }
    }
    
    void addNode(Node n) {
        _nodes.push_back(n);
    }
    
    void initializeAdjacencyList() {
        _edges.resize(_nodes.size());
    }
    
    void addEdge(Edge e) {
        _edges[e._s].push_back(e);
    }

    EdgeList& getNeighbors(Node n) {
        return _edges[n._id];
    }
    
    void kruskal_mst() {
        // sort edges in the increasing order
        std::sort(_edges.begin(), _edges.end(), Edge::weight_comp);
        
        // collect edges without cycle
        for (int j = 0; j < _edges.size(); j++) {
            
        }
        
    }
};




#endif /* defined(__interviews__graph__) */
