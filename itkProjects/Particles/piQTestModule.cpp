//
//  piQTestModule.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/4/13.
//
//

#include "piQTestModule.h"
#include "piParticleTrace.h"
#include "piFitCurve.h"

using namespace std;

namespace pi {
    QTestModule::QTestModule(QObject* parent): QObject(parent) {

    }

    QTestModule::~QTestModule() {
        
    }


    void QTestModule::FitTest(pi::StringVector &args) {
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

    void QTestModule::Run(Options *opts, StringVector &args) {
        if (opts->GetBool("--fitTest")) {
            FitTest(args);
        }
    }
}