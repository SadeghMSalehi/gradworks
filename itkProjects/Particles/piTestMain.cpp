#include "piTestMain.h"

#ifndef Q_MOC_RUN
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#endif

namespace pi {
    TestMain::TestMain() {
    }

    bool TestMain::main(Options& opts, StringVector& args) {
        if (opts.GetBool("--boost")) {
            testBoost(opts, args);
            return true;
        }
        return false;
    }

    bool TestMain::testBoost(Options &opts, StringVector &args) {
        using namespace std;
        typedef boost::numeric::ublas::vector<int> IntArray;


        // resize test
        IntArray testarr;

        testarr.resize(10);
        for (int i = 0; i < testarr.size(); i++) {
            testarr[i] = i;
        }

        for (int i = 0; i < testarr.size(); i++) {
            cout << testarr[i] << endl;
        }

        testarr.resize(20);
        for (int j = 10; j < testarr.size(); j++) {
            testarr[j] = j;
        }

        for (int i = 0; i < testarr.size(); i++) {
            cout << testarr[i] << endl;
        }

        return true;
    }
}
