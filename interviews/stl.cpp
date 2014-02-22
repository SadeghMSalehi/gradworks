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

int main(int argc, char* argv[]) {
//    test_list();
//    test_find_if();
//    test_foreach();
//    test_mem_fun();
//    test_stream_iterator();
    test_accumulation();
}
