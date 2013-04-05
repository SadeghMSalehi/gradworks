//
//  piQTestModule.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/4/13.
//
//

#include "piqTestModule.h"
#include "piParticleTrace.h"
#include "piFitCurve.h"
#include "piBSplineBasis.h"

using namespace std;
using namespace pi;

namespace piq {
    TestModule::TestModule(QObject* parent): QObject(parent) {

    }

    TestModule::~TestModule() {
        
    }


    void TestModule::BSplineTest() {
        BSplineBasis::TRealVector knots;
        int k = 4;
        for (int i = 0; i < k + 1; i++) {
            knots.push_back(i);
        }
        pi::BSplineBasis basis(&knots[0], knots.size(), 0, k);
        for (float t = 0; t < k; t += 0.1) {
            cout << basis.eval(0, 3, t) << "\t" << basis.eval(1, 3, t) << endl;
        }
    }
    
    void TestModule::FitTest(pi::StringVector &args) {
        ifstream is(args[0].c_str());
        ParticleVector particles;
        ParticleTrace::Read(is, particles);
        cout << "Read " << particles.size() << " particles ..." << endl;

        CurveFitting fit;
        fit.FitCurve(particles);
        ParticleVector& result = fit.GetResult();

        ofstream os("/tmpfs/curvefit.txt");
        ParticleTrace::Write(os, result);
        os.close();
    }

    void TestModule::Run(Options *opts, StringVector &args) {
        if (opts->GetBool("--fitTest")) {
            FitTest(args);
        } else if (opts->GetBool("--bsplineBasisTest")) {
            BSplineTest();
        }
    }
}