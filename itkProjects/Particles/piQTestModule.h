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

namespace pi {
    class QTestModule: public QObject {
        Q_OBJECT
    public:
        QTestModule(QObject* parent = NULL);
        virtual ~QTestModule();

        void FitTest(pi::StringVector &args);
        void Run(Options* opts, StringVector& args);
    };
}
#endif /* defined(__ParticleGuidedRegistration__piQTestModule__) */
