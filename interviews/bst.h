//
//  bst.h
//  interviews
//
//  Created by Joohwi Lee on 2/20/14.
//
//

#ifndef __interviews__bst__
#define __interviews__bst__

#include <iostream>
#include <utility>





template <class T1, class T2>
class TreeNode {
public:
    typedef std::pair<T1,T2> ContentType;

    TreeNode* left;
    TreeNode* right;

    ContentType value;

    TreeNode(): left(NULL), right(NULL) {}
};

template <class T1, class T2>
class BinarySearchTree {
public:
    typedef TreeNode<T1, T2> NodeType;
    typedef typename NodeType::ContentType ValueType;

    BinarySearchTree() : root(NULL) {}

    NodeType* search(T1 key);
    void insert(T1 key, T2 value);
    void remove(T1 key);
    void traverse();

private:
    void traverse(NodeType* iter);
    NodeType* root;
};



template <class T1, class T2>
void BinarySearchTree<T1,T2>::insert(T1 key, T2 value) {
    NodeType* currNode = root;
    NodeType* prevNode = NULL;

    bool left = false;
    while (currNode != NULL) {
        prevNode = currNode;
        if (currNode->value.first < key) {
            currNode = currNode->right;
            left = false;
        } else if (currNode->value.first > key) {
            currNode = currNode->left;
            left = true;
        }
    }

    if (prevNode == NULL) {
        // only for root
        root = new NodeType();
        root->value = std::make_pair(key, value);
    } else if (left) {
        prevNode->left = new NodeType();
        prevNode->left->value = std::make_pair(key, value);
    } else if (!left) {
        prevNode->right = new NodeType();
        prevNode->right->value = std::make_pair(key, value);
    }
}


template <class T1, class T2>
void BinarySearchTree<T1,T2>::traverse(NodeType* iter) {
    if (iter == NULL) {
        return;
    }
    traverse(iter->left);
    std::cout << " " << iter->value.first << " ";
    traverse(iter->right);
}

template <class T1, class T2>
void BinarySearchTree<T1,T2>::traverse() {
    traverse(root->left);
    std::cout << " " << root->value.first << " ";
    traverse(root->right);
}

#endif /* defined(__interviews__bst__) */