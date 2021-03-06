#ifndef PITESTMAIN_H
#define PITESTMAIN_H

#include "piOptions.h"

namespace pi {
    class TestMain
    {
    public:
        TestMain();

        bool main(Options& opts, StringVector& args);
        bool testBoost(Options& opts, StringVector& args);
        bool testNEWUOA(Options& opts, StringVector& args);
        bool testJSON(Options& opts, StringVector& args);
        bool testConfig(Options& opts, StringVector& args);
        bool testDisplacement(Options& opts, StringVector& args);
    };
}
#endif // PITESTMAIN_H
