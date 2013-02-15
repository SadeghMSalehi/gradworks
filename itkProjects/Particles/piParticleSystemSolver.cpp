//
//  piParticleSystemSolver.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/15/13.
//
//

#include "piParticleSystemSolver.h"
#include "piParticleCollision.h"
#include "piParticleForces.h"
#include "piOptions.h"

namespace pi {
    void ParticleSystemSolver::Preprocessing(ParticleSystem& system) {
        const int nSubz = system.GetNumberOfSubjects();
        const int nPoints = system.GetNumberOfParticles();
        Options& systemOptions = system.GetSystemOptions();
        ParticleSubjectArray& subs = system.GetSubjects();

        EntropyInternalForce internalForce;
        std::vector<ParticleCollision> collisionHandlers;
        collisionHandlers.resize(nSubz);

        for (int n = 0; n < nSubz; n++) {
            collisionHandlers[n].LoadBinaryMask(systemOptions.GetStringVectorValue("BinaryMask", n));
            collisionHandlers[n].UseBinaryMaskSmoothing();
            collisionHandlers[n].UseBinaryMaskSmoothingCache(systemOptions.GetStringVectorValue("BinaryMaskSmoothingCache", n).c_str());
            collisionHandlers[n].UseBinaryMaskSmoothingCache(systemOptions.GetStringVectorValue("BinaryMaskDistanceMapCache", n).c_str());
            collisionHandlers[n].UpdateImages();
        }

        for (double t = m_t0; t < m_t1; t += m_dt) {
            cout << "t: " << t << endl;
            for (int n = 0; n < nSubz; n++) {
                ParticleSubject& sub = subs[n];
                for (int i = 0; i < nPoints; i++) {
                    Particle& pi = sub[i];
                    forset(pi.x, pi.w);
                    forfill(pi.f, 0);
                }

                internalForce.ComputeForce(subs);
                collisionHandlers[n].HandleCollision(subs);

                for (int n = 0; n < nSubz; n++) {
                    ParticleSubject& sub = subs[n];
                    for (int i = 0; i < nPoints; i++) {
                        Particle& p = sub[i];
                        LabelImage::IndexType pIdx;
                        fordim (k) {
                            p.f[k] -= p.v[k];
                            p.v[k] += m_dt*p.f[k];
                            p.x[k] += m_dt*p.v[k];
                            pIdx[k] = p.x[k] + 0.5;
                        }

                        if (collisionHandlers[n].IsBufferInside(pIdx)) {
                            cout << "Stop system: out of region" << endl;
                            goto quit;
                        }
                    }
                }
            }
        }

    quit:
        return;
    }
}