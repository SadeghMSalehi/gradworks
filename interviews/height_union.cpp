//
//  gl.cpp
//  interviews
//
//  Created by Joohwi Lee on 2/22/14.
//
//

#include "height_union.h"
#include <algorithm>
#include <queue>
#include <vector>
#include <list>
#include <stack>

using namespace std;

struct Apt {
    double s;
    double e;
    double h;

    Apt(double a, double b, double c): s(a), e(b), h(c) {}
};


struct AptCompS {
    bool operator()(const Apt& a, const Apt& b) {
        return a.s < b.s;
    }
};

struct AptCompE {
    bool operator()(const Apt& a, const Apt& b) {
        return a.e < b.e;
    }
};


struct AptCompH {
    bool operator()(const Apt& a, const Apt& b) {
        return a.h < b.h;
    }
};

typedef priority_queue<Apt, vector<Apt>, AptCompH> AptHeap;

void processOpening(AptHeap& currents, Apt& a) {
    cout << "Open " << a.s << ": " << a.h << endl;
    currents.push(a);
}

void processClosing(AptHeap& currents, Apt& a) {
    cout << "Closing " << a.e << ": " << a.h << endl;
    while (currents.top().h == a.h) {
        currents.pop();
    }
}


void height_union2() {
    double data[3*7] = { 0, 10, 3, 9, 11, 0.5, 2, 7, 8, 3, 4, 1, 3.5, 4, 18, 3.7, 8, 12, 5, 10, 7 };

    deque<Apt> apts1, apts2;

    AptHeap currentIntervals;

    /// preparation
    for (int i = 0; i < 7; i++) {
        apts1.push_back(Apt(data[i*3], data[i*3+1], data[i*3+2]));
        apts2.push_back(Apt(data[i*3], data[i*3+1], data[i*3+2]));
    }

    /// sort Apt starting position list in the order of start position
    sort(apts1.begin(), apts1.end(), AptCompS());

    /// sort Apt ending position list in the order of ending position
    sort(apts2.begin(), apts2.end(), AptCompE());


    while (apts2.size() > 0){
        /// if an interval opens, store it into a list for later check
        if (apts1.front().s < apts2.front().e) {
            /// process the opening
            Apt opening = apts1.front();
            apts1.pop_front();

            /// check if the new opening is the highest among currently opened intervals
            bool highest = currentIntervals.top().h <= opening.h;
            if (highest) {
                currentIntervals.push(opening);
                cout << "Open " << opening.s << ": " << opening.h << endl;
            }
        } else if (apts1.front().s > apts2.front().e) {
            /// process the closing
            Apt closing = apts2.front();
            apts2.pop_front();

            /// check if this closing point is the highest among currently opened intervals
            bool highest = currentIntervals.top().h == closing.h && currentIntervals.top().e == closing.e;
            /// pop the closing from currentIntervals
            if (highest) {
                /// Assume the order is preserved
                currentIntervals.pop();
                cout << "Closing " << closing.e << ": " << closing.h << endl;
            }
        } else {
            /// process both (first closing, then opening)
//            handleOpening(apts1, apts2);
//            handleClosing(apts1, apts2);
        }
    }
}


struct Rect {
    double s;
    double e;
    double h;
};

struct RectCompS {
    bool operator()(const Rect* a, const Rect* b) {
        return a->s <= b->s;
    }
};

struct RectCompE {
    bool operator()(const Rect* a, const Rect* b) {
        return a->e <= b->e;
    }
};

struct RectCompH {
    bool operator()(const Rect* a, const Rect* b) {
        return a->h > b->h;
    }
};

struct RectPrinter {
    void operator()(const Rect* a) {
        cout << a->s << ", " << a->e << ": " << a->h << endl;
    }
};

void triggerStartingEvent(const Rect* e) {
    cout << e->s << ": " << e->h << endl;
}

void triggerEndingEvent(const Rect* e) {
    cout << e->e << ": " << e->h << endl;
}

void height_union() {
    using namespace std;

    Rect data[7] = {
        { 0, 10, 3 },
        { 9, 11, 0.5 },
        { 2, 7, 8 },
        { 3, 4, 1 },
        { 3.5, 4, 18 },
        { 3.7, 8, 12 },
        { 5, 10, 7 }
    };

    std::vector<Rect*> startingPos, endingPos;
    for (int i = 0; i < 7; i++) {
        startingPos.push_back(&data[i]);
        endingPos.push_back(&data[i]);
    }

    std::sort(startingPos.begin(), startingPos.end(), RectCompS());
    std::sort(endingPos.begin(), endingPos.end(), RectCompE());

    cout << "Input1: " << endl;
    std::for_each(startingPos.begin(), startingPos.end(), RectPrinter());
    cout << "Input2: " << endl;
    std::for_each(endingPos.begin(), endingPos.end(), RectPrinter());

    stack<Rect*, vector<Rect*> > startingStack(startingPos);
    stack<Rect*, vector<Rect*> > endingStack(endingPos);
    priority_queue<Rect*, vector<Rect*>, RectCompH> openIntervalStack;

    while (!endingStack.empty()) {
        if (startingStack.top()->s < endingStack.top()->e) {
            Rect* rect = startingStack.top();
            if (rect->h > openIntervalStack.top()->h) {
                triggerStartingEvent(rect);
            }
            openIntervalStack.push(rect);
            startingStack.pop();
        } else {
            Rect* rect = endingStack.top();
            if (rect->h == openIntervalStack.top()->h) {
                if (rect == openIntervalStack.top()) {
                    triggerEndingEvent(rect);
                    openIntervalStack.pop();
                }
            } else {
//                openIntervalStack.remove(rect);
            }
            endingStack.pop();
        }
    }
}