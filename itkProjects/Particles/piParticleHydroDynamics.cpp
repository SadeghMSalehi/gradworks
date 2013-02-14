//
//  piParticleHydroDynamics.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/14/13.
//
//

#include "piParticleHydroDynamics.h"


namespace pi {
    const static double h = 1;
    
    double ComputeDensityKernel(Particle& pi, Particle& pj, double rij) {
        return 0;
    }

    double ComputeDensityKernelDerivative(Particle& pi, Particle& pj, double rij) {
        return 0;
    }

    double ComputePressureKernel(Particle& pi, Particle& pj, double rij) {
        return 0;
    }

    double ComputePressureKernelDerivative(Particle& pi, Particle& pj, double rij) {
        return 0;
    }

    void HydroDynamics::ComputeValues(ParticleSubject& subj) {
        const int nPoints = subj.GetNumberOfPoints();
        for (int i = 0; i < nPoints; i++) {
            subj.m_Particles[i].density = 0;
        }
        for (int i = 0; i < nPoints; i++) {
            Particle& pi = subj.m_Particles[i];
            for (int j = i + 1; j < nPoints; j++) {
                Particle& pj = subj.m_Particles[j];
                double rij = sqrt(pi.Dist2(pj));
                pi.density += ComputeDensityKernel(pi, pj, rij);
            }
        }
    }
}