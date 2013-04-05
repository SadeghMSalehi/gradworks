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
        int k = 3;
        int np = 3;
        BSplineBasis::TReal knots[7] = { 0, 0, 0, 1, 2, 2, 2 };
        pi::BSplineBasis b1(knots, 7, 0, 3);
        pi::BSplineBasis b2(knots, 7, 1, k);
        pi::BSplineBasis b3(knots, 7, 2, k);
        pi::BSplineBasis b4(knots, 7, 3, k);
        for (float t = 0; t < 2.1; t += 0.1) {
            cout << t << "\t" <<  b1(t) << "\t" << b2(t) << "\t" << b3(t) << "\t" << b4(t) << endl;
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