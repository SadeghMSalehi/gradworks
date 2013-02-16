//
//  myParticleForces.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 1/26/13.
//
//

#ifndef __ParticlesGUI__myParticleForces__
#define __ParticlesGUI__myParticleForces__

#include <iostream>
#include "piParticleCore.h"

namespace pi {
    class InternalForce {
    public:
        InternalForce() {}
        ~InternalForce() {}
        void ComputeForce(ParticleSubject& subj);
        void ComputeForce(ParticleSubjectArray& shapes);
        void ComputeForce(Particle& a, Particle& b);
    };

    class EntropyInternalForce {
    public:
        EntropyInternalForce() {}
        ~EntropyInternalForce() {}

        void ComputeForce(ParticleSubject& subj);
        void ComputeForce(ParticleSubjectArray& subjs);
        void ComputeForce(Particle& a, Particle& b);
    };

    class EnsembleForce {
    public:
        EnsembleForce(double coeff);
        ~EnsembleForce();
        void SetImageContext(ImageContext* context);
        void ComputeEnsembleForce(ParticleSubjectArray& shapes);
        void ComputeImageForce(ParticleSubjectArray& shapes);
    private:
        ImageContext* m_ImageContext;
        ParticleSubject m_MeanShape;
        void ComputeMeanShape(ParticleSubjectArray& shapes);
        double m_Coeff;
    };

    class IntensityForce {
    public:
        IntensityForce(double coeff);
        ~IntensityForce();
        void SetImageContext(ImageContext* context);
        void ComputeIntensityForce(ParticleSystem* system);
    private:
        ImageContext* m_ImageContext;
        double m_Coeff;
    };
}

#endif /* defined(__ParticlesGUI__myParticleForces__) */
