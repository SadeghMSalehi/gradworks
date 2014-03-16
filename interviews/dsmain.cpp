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
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

using namespace std;

void bst_test() {
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

void permute(string p, string m) {
    static int c = 0;

    // trivial base case
    if (m.length() <= 1) {
        return;
    }

    // base case
    if (m.length() == 2) {
        cout << ++c << " ";
        cout << p << m << endl;
        cout << (++c) << " ";
        cout << p << m[1] << m[0] << endl;
    } else {
        for (int i = 0; i < m.length(); i++) {
            permute(p + m[i], m.substr(0, i) + m.substr(i + 1, m.length() - 1 - i));
        }
    }
}

void string_permutation_test() {
    string a = "12345";

    permute("", a);
}

void last_klines(const char* filename) {
    ifstream in(filename);

    int n = 0;
    std::vector<string> queue;
    queue.resize(10);

    char buff[1024];
    while (in.is_open()) {
        in.getline(buff, sizeof(buff));
        if (in.good() && !in.eof()) {
            string s(buff);
            ++n;
            if (n == 10) {
                n = 0;
            }
            queue[n] = s;
        } else {
            break;
        }
    }
    in.close();

    for (int i = 0; i < 10; i++) {
        cout << queue[n] << endl;
        n++;
        if (n == 10) {
            n = 0;
        }
    }
}

void run_length_encoding() {
    const char* input = "abcccccaaa";
    const char* p = input;

    char c = 0;
    int count = 0;
    for (p = input; *p != (char) NULL; ++p) {
        if (*p == c) {
            count++;
        } else if (c > 0) {
            cout << c << (count + 1);
            count = 0;
        }
        c = *p;
    }
    cout << c << (count + 1);
}

void remove_dups() {
    list<char> ls;
    const char* str = "FOLLOW UP";
    for (const char* p = str; p != NULL; ++p) {
        ls.push_back(*p);
    }

    list<char>::iterator iter = ls.begin();
    list<char>::iterator prevIter;

    for (; iter != ls.end(); ++iter) {
        if (iter == ls.begin()) {
            prevIter = ls.begin();
            continue;
        } else {
            if (*iter == *prevIter) {
                continue;
            } else {
                // delete between (prevIter ~ iter)
            }
        }
    }
}

void big_sum() {
    int nums[2] = { 1234, 9999 };

    list<int> numList[2];
    for (int i = 0; i < 2; i++) {
        int denom = 1, remainder = 0;

        denom = nums[i];
        while (denom > 0) {
            remainder = denom % 10;
            denom = denom / 10;
            numList[i].push_back(remainder);
        }
    }

    list<int>::iterator iter[2];
    iter[0] = numList[0].begin();
    iter[1] = numList[1].begin();

    list<int> sol;

    int carry = 0;
    while (! (iter[0] == numList[0].end() && iter[1] == numList[1].end())) {
        int a = 0, b = 0;
        if (iter[0] != numList[0].end()) {
            a = *iter[0];
            ++iter[0];
        }
        if (iter[1] != numList[1].end()) {
            b = *iter[1];
            ++iter[1];
        }

        int digitSum = a + b;
        int newCarry = 0;
        if (digitSum > 10) {
            newCarry = 1;
            digitSum = (a + b) % 10;
        }

        sol.push_back(digitSum + carry);
        carry = newCarry;
    }

    std::reverse(sol.begin(), sol.end());
    list<int>::iterator solIter = sol.begin();

    for (; solIter != sol.end(); ++solIter) {
        cout << *solIter;
    }
    cout << endl;


}

void find_loop() {
    struct node {
        int val;
        node* next;

        node(int v) : val(v), next(NULL) {
        }
    };

    struct linked_list {
        node* head;

        linked_list(): head(NULL) {}

        void push_back(node* n) {
            if (head == NULL) {
                head = n;
                return;
            }
            node* iter = head;
            while (iter->next != NULL) {
                iter = iter->next;
            }
            iter->next = n;
        }
    };

    linked_list l;
    l.push_back(new node(1));
    l.push_back(new node(2));
    l.push_back(new node(3));
    node* c = new node(4);
    l.push_back(c);
    l.push_back(new node(5));
    l.push_back(new node(6));
    l.push_back(new node(7));
    l.push_back(c);

    node* iter = l.head;
    for (int i = 0; i < 20; i++) {
        if (iter == NULL) {
            break;
        }
        cout << iter->val << " ";
        iter = iter->next;
    }

    node* iter1 = l.head;
    node* iter2 = l.head->next;

    cout << endl;
    int count = 1;
    while (iter1->val != iter2->val) {
        cout << ++count << " ";
        iter1 = iter1->next;
        iter2 = iter2->next->next;
    }

    cout << endl;
    cout << iter1->val << endl;
    cout << iter2->val << endl;

}


static int __idx;

class Attr {
public:
    Attr() {
        _z = ++__idx;
    }
    virtual ~Attr() {};

    virtual int z() {
        return _z;
    }

private:
    int _z;

};

template<class S>
struct node: public S {
    int x;
    node(int n) {
        x = n;
    }
    void operator()(const node<S>& s) {
        cout << s.x << " ";
    }
};

template<class S>
struct RichVector: public vector<S> {
    ~RichVector() {}

    void printAll(const char* header = NULL) {
        if (header) {
            cout << header << ": ";
        }
        for (int i = 0; i < this->size(); i++) {
            cout << (*this)[i].x << " ";
        }
        cout << endl;
    }

    void swap(int i, int j) {
        S tmp = (*this)[i];
        (*this)[i] = (*this)[j];
        (*this)[j] = tmp;
    }

    int partition(int startPos, int endPos, int pivotIdx) {
        // swap elements with respect to the pivot k
        // 1) swap the pivot and the last element
        S pivot = (*this)[pivotIdx];
        swap(pivotIdx, endPos);
        int newPivot = startPos;
        for (int i = startPos; i < endPos; i++) {
            if ((*this)[i].x <= pivot.x) {
                swap(i, newPivot);
                newPivot ++;
            }
        }
        swap(newPivot, endPos);
        return newPivot;
    }

    void quickSort(int startPos, int endPos, int pivotIdx) {
        if (startPos == endPos) {
            return;
        }
        int newPivot = partition(startPos, endPos, pivotIdx);

        if (startPos < newPivot - 1) {
            quickSort(startPos, newPivot - 1, rand() % (newPivot - startPos) + startPos);
        }
        if (newPivot + 1 < endPos) {
            quickSort(newPivot + 1, endPos, rand() % (endPos - newPivot) + newPivot + 1);
        }
    }

    S& binarySearch(S s, int startPos, int endPos) {

        int mid = (endPos - startPos) / 2 + startPos;

        /// compare at i
        if (this->at(mid).x == s.x) {
            return this->at(mid);
        }
        if (this->at(mid).x < s.x) {
            return binarySearch(s, mid + 1, endPos);
        } else {
            return binarySearch(s, startPos, mid - 1);
        }
    }
};

void sort_test() {
    typedef RichVector<node<Attr> > NodeVector;

    NodeVector list;
    for (int i = 0; i < 10; i++) {
        node<Attr> a(rand() % 1000);
        list.push_back(a);
        cout << list[i].z() << " ";
    }

    list.printAll();

    /// bubble sort
    /// find the largest unordered element and propagate to the end
    NodeVector bubbleSort = list;
    for (int i = 0; i < bubbleSort.size(); i++) {
        for (int j = 1; j < bubbleSort.size() - i; j++) {
            if (bubbleSort[j].x < bubbleSort[j-1].x) {
                node<Attr> tmp = bubbleSort[j];
                bubbleSort[j] = bubbleSort[j-1];
                bubbleSort[j-1] = tmp;
            }
        }
    }
    bubbleSort.printAll("Bubble Sort");


    /// Perform the insertion sort
    NodeVector insertionSort = list;
    /// Assume that the first k-elements are sorted
    for (int i = 1; i < insertionSort.size(); i++) {
        for (int j = i; j > 0; j--) {
            if (insertionSort[j].x < insertionSort[j - 1].x) {
                node<Attr> tmp = insertionSort[j - 1];
                insertionSort[j - 1] = insertionSort[j];
                insertionSort[j] = tmp;
            }
        }
    }
    insertionSort.printAll("Insertion Sort");


    /// Perform quick sort?
    NodeVector quickSort = list;
    quickSort.quickSort(0, quickSort.size() - 1, 5);
    quickSort.printAll("Quick Sort");

    cout << quickSort.binarySearch(node<Attr>(878), 0, 9).z() << endl;
}



void print_matrix(int mat[][5], int m, int n) {
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < 5; i++) {
            cout << mat[j][i] << " ";
        }
        cout << endl;
    }
}

void fill(int mat[][5], int x, int y, int o, int v) {
    cout << x << "," << y << endl;
    if (x < 0 || x >= 5 || y < 0 || y >= 5) {
        return;
    }
    mat[y][x] = v;
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            if (i == 0 && j == 0) {
                continue;
            } else if (mat[y+i][x+j] == o) {
                fill(mat, x+j, y+i, o, v);
            }
        }
    }
}

void fill_test() {
    int matrix[5][5] = { { 1, 2, 3, 4, 5 }, {4, 2, 5, 9, 1 }, {4, 2, 3, 5, 1}, {4, 2, 2, 2, 5}, {2, 2, 2, 2, 2} };
    print_matrix(matrix, 5, 5);
    fill(matrix, 1, 1, matrix[1][1], 10);
    print_matrix(matrix, 5, 5);
}

int main(int argc, char* argv[]) {
//    bst_test();
//    string_permutation_test();
//    const char* filename = "/Users/joohwi/ares.sh";
//    last_klines(filename);
//    run_length_encoding();
//    remove_dups();
//    big_sum();
//    find_loop();
//    sort_test();
    fill_test();
}