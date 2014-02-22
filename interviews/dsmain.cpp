//
//  dsmain.cpp
//  interviews
//
//  Created by Joohwi Lee on 2/20/14.
//
//

#include "dsmain.h"
#include "bst.h"

#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]) {
    /* initialize random seed: */
    srand(time(NULL));

    BinarySearchTree<int, int> bst;
    for (int i = 0; i < 10; i++) {
        // insert random numbers into the tree
        int v = rand() % 1000 + 1;
        cout << v << ",";
        bst.insert(v, v);
    }

    cout << "Traverse: ";
    bst.traverse();
}