#include "stl.h"
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iterator>


using namespace std;

ostream& operator<<(ostream& out, const list<int>& input) {
    list<int>::const_iterator iter;
    out << "{ ";
    for (iter = input.begin(); iter != input.end(); ++iter) {
        if (iter != input.begin()) {
            out << ",";
        }
        out << *iter;
    }
    out << " } " << endl;
    return out;
}

void test_list() {
    list<int> a;
    a.push_back(9);
    a.push_back(7);
    a.push_back(5);
    a.push_back(3);
    a.push_back(1);

    a.sort();
    cout << a << endl;

    list<int> b;
    b.push_back(2);
    b.push_back(4);
    b.push_back(6);

    a.merge(b);
    cout << a << endl;
    cout << b << endl;
}

class LessThan {
public:
    int n;

    LessThan(int k): n(k) {

    }

    bool operator()(int x) {
        return x < n;
    }
};

void test_find_if() {
    list<int> a;
    for (int i = 0; i < 20; i++) {
        a.push_back(i);
    }

    // find the first element that satisfies the given predicate
    list<int>::iterator iter;
    for(iter = a.begin(); iter != a.end(); ++iter) {
        iter = find_if(iter, a.end(), LessThan(10));
        if (iter == a.end()) {
            break;
        }
        cout << *iter << endl;
    }
}

template <class T>
class MeanStd {
public:
    int n;
    double sum, sum2;
    MeanStd(): n(0), sum(0), sum2(0) {}
    void operator()(T x) {
        sum += x;
        sum2 += x*x;
        n++;
    }
    double mean() {
        return sum / n;
    }
    double std() {
        double m = mean();
        return (sum2) / n - m*m;
    }
};

void test_foreach() {
    list<int> a;
    for (int i = 0; i < 20; i++) {
        a.push_back(i);
    }

    MeanStd<int> mstd;
    mstd = for_each(a.begin(), a.end(), mstd);
    cout << "Mean: " << mstd.mean() << endl;
    cout << "Std: " << mstd.std() << endl;
}


class MemberFun {
public:
    int n;
    MemberFun(int x) : n(x) {}
    int value() {
        cout << n << ",";
        return n;
    }
};

void test_mem_fun() {
    list<MemberFun> l;
    for (int i = 0; i < 10; i++) {
        l.push_back(MemberFun(i));
    }
    for_each(l.begin(), l.end(), mem_fun_ref(&MemberFun::value));
}


void test_stream_iterator() {
    ostream_iterator<float> os(cout);
    *os = 10.1;
    ++os;
    *os = 20.5;
}

void test_accumulation() {
    vector<string> x;
    x.push_back("hello");
    x.push_back("world");
    x.push_back("~~ New World");
    cout << accumulate(x.begin(), x.end(), string("")) << endl;
    partial_sum(x.begin(), x.end() - 1, back_inserter(x));
    cout << x.size() << endl;
}

void test_insertion_sort() {
    list<int> l;
    
    int data[] = { 80, 33, 49, 101, 72, 55 };
    
    list<int> o;
    
    for (int i = 0; i < sizeof(data) / sizeof(int); i++) {
        list<int>::iterator iter = o.begin();
        while (*iter < data[i] && iter != o.end()) {
            ++iter;
        }
        
        o.insert(iter, data[i]);
    }
    
    list<int>::iterator iter = o.begin();
    for (; iter != o.end(); ++iter) {
        cout << *iter << endl;
    }

}

void find_magic_index() {
    
    int data[] = { -40, -20, -1, 1, 2, 3, 5, 7, 9, 12, 13 };
    int ndata = sizeof(data) / sizeof(int);
    
    bool stop = false;
    int idx = ndata / 2;
    
    while (!stop) {
        if (data[idx] == idx) {
            cout << idx << endl;
            break;
        }
        int pidx = idx;
        if (data[idx] > idx) {
            idx = (idx / 2);
        } else if (data[idx] < idx) {
            idx = idx + (ndata - idx) / 2;
        }
        if (idx == pidx) {
            cout << "can't find" << endl;
            break;
        }
    }
}

typedef pair<int,int> person;
bool person_comp(const person& p1, const person& p2) {
    return p1.first < p2.first;
}

void find_longest_subsequences() {

    vector<person> list;
    
    list.push_back(make_pair(2, 160));
    list.push_back(make_pair(3, 180));
    for (int i = 441; i < 485; i++) {
        list.push_back(make_pair(i, i));
    }
    for (int i = 0; i < 30; i++) {
        int x = rand() % 1000;
        int y = rand() % 1000;
        list.push_back(make_pair(x, y));
    }
    
    // sort for x
    std::sort(list.begin(), list.end(), person_comp);
    
    vector<person>::iterator iter = list.begin();
    for (int i = 0; iter != list.end(); ++iter, i++) {
        cout << iter->first << "," << iter->second << endl;
    }
    
    
    // find longest increasing interval in y
    int lastConvexIdx = 0;
    int lastConcaveIdx = 0;
    int sequenceLength = 0, maxSequenceLength = -1;
    
    for (int i = 2; i < list.size(); i++) {
        if (list[i].second > list[i-1].second && list[i-1].second <= list[i-2].second) {
            // starting point
            lastConvexIdx = i - 1;
            sequenceLength = 0;
        } else if (list[i].second <= list[i-1].second && list[i-1].second > list[i-2].second) {
            lastConcaveIdx = i - 1;
            if (sequenceLength > maxSequenceLength) {
                cout << "update max!!: " << sequenceLength << "<=" << maxSequenceLength << endl;
                maxSequenceLength = sequenceLength;
            }
        } else if (list[i].second > list[i-1].second && list[i-1].second > list[i-2].second) {
            if (sequenceLength == 0) {
                sequenceLength = 3;
            } else {
                // strictly increasing order
                sequenceLength++;
            }
        } else {
            sequenceLength = 0;
        }
    }
    
    if (sequenceLength > maxSequenceLength) {
        maxSequenceLength = sequenceLength;
    }

    cout << "max: " << maxSequenceLength << endl;
}

template <typename T>
struct ListNode {
    typedef shared_ptr<ListNode> Pointer;
    
    T value;
    Pointer next;
    
    ListNode(T t): value(t) {}
    virtual ~ListNode() {
        cout << "delete " << value << endl;
    }
    
    static Pointer New(T t) {
        return Pointer(new ListNode<T>(t));
    }
};

void test_shared_ptr() {
    typedef ListNode<int> IntNode;
    
    IntNode::Pointer head = IntNode::New(0);
    IntNode::Pointer iter = head;
    for (int i = 1; i < 10; i++) {
        IntNode::Pointer next = IntNode::New(i);
        iter->next = next;
        iter = iter->next;
    }
    
    iter = head;
    while (iter.get() != NULL) {
        cout << iter->value << endl;
        iter = iter->next;
    }
}

void test_bst() {
    
}

template <typename K, typename V>
struct BinaryNode {
    typedef shared_ptr<BinaryNode> Pointer;
    
    K key;
    V value;
    
    Pointer smaller;
    Pointer larger;
    
    BinaryNode(K k, V v): key(k), value(v) {
    }
    
    static Pointer New(K k, V v) {
        return Pointer(new BinaryNode<K,V>(k,v));
    }
};

template <typename K, typename V>
struct BinarySearchTree {
    typedef BinaryNode<K,V> NodeType;
    typedef typename NodeType::Pointer NodePointer;
    
    NodePointer root;
    
    void Insert(K k, V v) {
        
    }
    
    void Delete(K k) {
        
        
    }

    NodePointer Search(K k) {
        return Search(k, root);
    }
    
    NodePointer Search(K k, NodePointer r) {
        if (r->key == k) {
            return r;
        } else if (r->key > k) {
            return Search(k, r->larger);
        } else {
            return Search(k, r->smaller);
        }
    }
};


struct member {
    int pid;
    int id;
    string name;
    
    member(string s, int i, int p): pid(p), id(i), name(s) {};
};

bool member_sort_by_id(const member& a, const member& b) {
    return a.id < b.id;
}

void test_org_tree() {
    vector<member> l;
    l.push_back(member("Tom", 60, 10));
    l.push_back(member("Jack", 65, 10));
    l.push_back(member("Boss", 10, 0));
    l.push_back(member("John", 100, 60));
    
    sort(l.begin(), l.end(), member_sort_by_id);

    vector<member>::iterator iter = l.begin();
    for (; iter != l.end(); ++iter) {
        cout << iter->id << ", " << iter->name << endl;
    }
}

int main(int argc, char* argv[]) {
//    test_list();
//    test_find_if();
//    test_foreach();
//    test_mem_fun();
//    test_stream_iterator();
//    test_accumulation();
//    test_insertion_sort();
    
//    find_magic_index();
//    find_longest_subsequences();
//    test_shared_ptr();
//    test_bst();
    test_org_tree();
}
