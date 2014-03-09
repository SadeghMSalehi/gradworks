//
//  heap.h
//  interviews
//
//  Created by Joohwi Lee on 2/20/14.
//
//

#ifndef __interviews__heap__
#define __interviews__heap__

#include <iostream>
#include <utility>
#include <vector>

template <class TKey, class TValue>
class Heap {
public:
    typedef std::pair<TKey,TValue> Node;

private:
    std::vector<Node> heap;
};

#endif /* defined(__interviews__heap__) */
