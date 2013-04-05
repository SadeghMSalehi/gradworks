//
//  piQTestModule.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/4/13.
//
//

#ifndef __ParticleGuidedRegistration__piQTestModule__
#define __ParticleGuidedRegistration__piQTestModule__

#include <iostream>
#include <QObject>
#include "piOptions.h"

namespace piq {
    class TestModule: public QObject {
        Q_OBJECT
    public:
        TestModule(QObject* parent = NULL);
        virtual ~TestModule();

        void BSplineTest();
        void FitTest(pi::StringVector &args);
        void Run(pi::Options* opts, pi::StringVector& args);
    };
}
#endif /* defined(__ParticleGuidedRegistration__piQTestModule__) */
