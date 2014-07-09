//
//  trie.h
//  interviews
//
//  Created by Joowhi Lee on 4/4/14.
//
//

#ifndef __interviews__trie__
#define __interviews__trie__

#include <iostream>
#include <string>


class TrieNode {
public:
    // assume that only lower cases characters are stored
    std::string prefix;
    TrieNode* alphabets[27];
    TrieNode* child;
};

template <class V>
class TrieValue {
public:
private:
    V value;
};


template <class V>
class Trie {
public:
    void add(std::string k, V v);
    bool remove(std::string k);
    bool has(std::string k);
    V operator[](std::string key);
};

class NoKeyException {
    
};


template <class V>
void Trie<V>::add(std::string k, V v) {
    
}

template <class V>
bool Trie<V>::remove(std::string k) {
    if (!has(k)) {
        return false;
    }
}

template <class V>
V Trie<V>::operator[](std::string k) {
    if (!has(k)) {
        throw NoKeyException();
    }
}



void trie_test();

#endif /* defined(__interviews__trie__) */
