#include "piTestMain.h"
#include "piPowellOpti.h"
#include "piImageDef.h"
#include "piImageIO.h"

#ifndef Q_MOC_RUN
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#endif

#include "piOptionParser.h"
#include "piConfigFile.h"
#include <libconfig.h++>

namespace pi {
    TestMain::TestMain() {
    }

    bool TestMain::main(Options& opts, StringVector& args) {
        if (opts.GetBool("--boost")) {
            testBoost(opts, args);
            return true;
        } else if (opts.GetBool("--newuoa")) {
            testNEWUOA(opts, args);
            return true;
        } else if (opts.GetBool("--testjson")) {
            testJSON(opts, args);
            return true;
        } else if (opts.GetBool("--testconfig")) {
            testConfig(opts, args);
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


    class X2 {
    public:
        double operator()(int n, double* x) {
            return (x[0]-3)*(x[0]-3) + 1;
        }
    };


    class Median {
    public:
        double operator()(int n, double* x) {
            static const int d[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double f = 0;
            for (int i = 0; i < 10; i++) {
                f += abs(x[0] - d[i]);
            }
            return f;
        }
    };

    bool TestMain::testNEWUOA(pi::Options &opts, StringVector &args) {
        using namespace std;

        typedef Median F;
        PowellOpti<F> fMin;
        PowellParams fInitial;
        fInitial.resize(1);
        fInitial[0] = (args.size() > 0 ? atof(args[0].c_str()) : 0);

        F f;
        double minValue = fMin.minimizeNEWUOA(f, fInitial);
        cout << "f(x_min) = " << minValue << endl;
        for (int i = 0; i < fInitial.size(); i++) {
            cout << fInitial[i] << endl;
        }

        return true;
    }



    bool TestMain::testJSON(pi::Options &opts, StringVector &args) {
        using namespace std;
        cout << "Parsing: " << opts.GetString("--json") << endl;
        OptionParser jsonParser;
        jsonParser.read(opts.GetString("--json"));

        cout << jsonParser.stringValue() << endl;
        return true;
    }

    bool TestMain::testConfig(pi::Options &opts, StringVector &args) {
        using namespace libconfig;
        using namespace std;

        Config config;
        config.readFile(opts.GetString("--config").c_str());

        Setting& images = config.lookup("images");
        for (int i = 0; i < images.getLength(); i++) {
            cout << images[i].c_str() << endl;
        }


        return true;
    }

    bool TestMain::testDisplacement(pi::Options &opts, StringVector &args) {
        return false;
    }
}