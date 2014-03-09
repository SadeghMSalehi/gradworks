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



template <class T>
struct listnode {
    T data;
    listnode* next;

    listnode(T d): data(d), next(NULL) {}

    listnode* remove_next() {
        if (next == NULL) {
            return;
        }
        listnode* oldNext = next;
        next = this->next->next;
        return oldNext;
    }

    void append(listnode* l) {
        if (l == NULL) {
            return;
        }
        if (next == NULL) {
            next = l;
        } else {
            listnode *p = l;
            while (p->next != NULL) {
                p = p->next;
            }
            p->next = this->next;
            this->next = l;
        }
    }
};


template <class T>
struct qlist {
    typedef listnode<T> nodetype;

    qlist(): _head(NULL), _tail(NULL) {}

    virtual ~qlist() {
        // delete every listnode
    }

    nodetype* head() {
        return _head;
    }

    nodetype* tail() {
        return _tail;
    }

    template <class U>
    void find_all(U u, qlist<T>* l) {
        nodetype* iter = _head;

        while (iter != NULL) {
            if (u(iter->data)) {
                l->push_back(iter->data);
            }
            iter = iter->next;
        }
    }

    void push_front(T t) {
        if (_head == NULL && _tail == NULL) {
            _head = _tail = new nodetype(t);
        } else {
            nodetype* n = new nodetype(t);
            n->next = _head;
            _head = n;
        }
    }

    void push_back(T t) {
        if (_tail == NULL) {
            _head = _tail = new nodetype(t);
        } else {
            _tail->next = new nodetype(t);
            _tail = _tail->next;
        }
    }

    void push_front(qlist<T>* l) {
        if (l->_tail == NULL) {
            return;
        }

        // keep the old head to connect later
        nodetype* oldhead = _head;

        // prepare to copy the given list
        nodetype* srcIter = l->_head;

        // create a new head to begin
        _head = new nodetype(srcIter->data);

        // set up a temporary tail pointer
        nodetype* tmptail = _head;

        // proceed the src iterator
        srcIter = srcIter->next;

        // prepare a new iterator
        nodetype* iter = _head;

        // while the srcIter doesn't reach at the end

        while (srcIter != NULL) {
            iter->next = new nodetype(srcIter->data);
            srcIter = srcIter->next;
            tmptail = iter;
            iter = iter->next;
        }

        // connect newly created list and the previous one;
        tmptail->next = oldhead;
    }

    void push_back(qlist<T>* l) {
        nodetype* srcIter = l->_head;
        if (srcIter == NULL) {
            // if a given list is empty, do nothing
            return;
        }

        nodetype* thisIter = this->_tail;
        if (thisIter == NULL) {
            // if this list is empty, create a head node and copy the head from anohter list.
            this->_head = this->_tail = new nodetype(srcIter->data);
        }

        srcIter = srcIter->next;
        while (srcIter != NULL) {
            // now both srcIter and thisIter is not null
            // copy the current srcIter to this iter's next
            thisIter->next = new nodetype(srcIter->data);

            // proceed this->tail
            this->_tail = thisIter->next;

            // process srcIter
            srcIter = srcIter->next;
        }
    }
    

    T pop_front() {
        if (_head == NULL) {
            return NULL;
        }

        nodetype* oldhead = _head;
        T v = _head->data;
        _head = _head->next;
        if (_head == NULL) {
            _tail = NULL;
        }

        delete oldhead;
        return v;
    }


    void remove(nodetype* node) {
        nodetype* iter = _head;
        if (iter == NULL) {
            return;
        }
        if (node == iter) {
            _head = iter->next;
            return;
        }
        while (iter != NULL && node != iter->next) {
            iter = iter->next;
        }
        if (iter == NULL) {
            // no node found
            return;
        }
        iter->next = iter->next->next;
        return;
    }

    bool is_empty() {
        return _head == NULL && _tail == NULL;
    }
    

private:
    nodetype* _head;
    nodetype* _tail;
};


template <class T>
struct treenode {
    treenode* parent;
    qlist<treenode*> children;
    T data;

    treenode<T>(T t): parent(NULL), data(t) {};
};


struct info {
    string name;
    int id;
    int pid;

    info(string n, int i, int p): name(n), id(i), pid(p) {}
};


struct pidcomp {
    int pid ;

    pidcomp(int p): pid(p) {
    }

    bool operator()(treenode<info>* node) {
        return node->data.pid == pid;
    }
};

struct infotree {
public:
    typedef treenode<info> infonode;

    infotree(): _root(NULL) {}

    infonode* find_by_id(int id) {
        qlist<treenode<info>*> queue;

        // always check the initial condition
        if (_root == NULL) {
            return NULL;
        }
        queue.push_back(_root);

        // pop the first
        while (!queue.is_empty()) {
            treenode<info>* n = queue.pop_front();

            // process n
            if (n->data.id == id) {
                return n;
            }

            // add n's children to queue
            queue.push_front(&(n->children));
        }

        return NULL;
    }

    void add_info(string name, int id, int pid) {
        infonode* newinfo = new infonode(info(name, id, pid));

        infonode* parent = find_by_id(pid);
        if (parent == NULL) {
            // add to reserve node
            tmpnodes.push_back(newinfo);
        } else {
            parent->children.push_back(newinfo);
        }

        // find available children from tmpnodes
        qlist<infonode*> children;
        tmpnodes.find_all(pidcomp(id), &children);

        // iterate over qlist and make its parent to the newinfo
        listnode<infonode*>* iter = children.head();

        while (iter != NULL) {
            // make these as children of the newinfo
            newinfo->children.push_back(iter->data);
            iter->data->parent = newinfo;
            tmpnodes.remove(iter);
            listnode<infonode*>* delnode = iter;
            iter = iter->next;
            cout << delnode->data->data.name << endl;
            delete delnode;
        }

        if (newinfo->data.pid == 0) {
            _root = newinfo;
        }
    }

    void print_tree() {
        
    }

private:
    infonode* _root;
    qlist<infonode*> tmpnodes;
};


void n_tree_test() {
    infotree tree;
    tree.add_info("Tom", 60, 10);
    tree.add_info("Jack", 65, 10);
    tree.add_info("Boss", 10, 0);
    tree.add_info("John", 100, 60);
}


void exception_test() {
    try {
        cout << "Hello!";
        return;
    } catch (...) {
        cout << " World!" << endl;
    }
    cout << "Done.." << endl;
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
//    fill_test();
    n_tree_test();
//    n_tree_test2();
//    exception_test();
}