//
//  myParticleForces.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 1/26/13.
//
//

#include "myParticleCore.h"
#include "myParticleBSpline.h"

namespace my {
    void InternalForce::ComputeForce(ParticleShapeArray& shapes) {
        int nPoints = shapes[0].m_nPoints;
        for (int k = 0; k < shapes.size(); k++) {
            ParticleArray& particles = shapes[k].m_Particles;
            for (int i = 0; i < nPoints; i++) {
                for (int j = i+1; j < nPoints; j++) {
                    ComputeForce(particles[i], particles[j]);
                }
            }
        }
    }
    
    void InternalForce::ComputeForce(my::Particle &a, my::Particle &b) {
        double dist = std::sqrt(a.Dist2(b));
        const double sigma = 7 * 5;
        const double coeff = M_PI_2 / sigma;
        
        double f[4] = { 0, };
        double dx[4];
        a.Sub(b, dx);
        
        if (dist > sigma) {
            return;
        }
        
        double rij = dist * coeff;
        if (rij > 0) {
            double sin2rij = std::sin(rij);
            sin2rij *= sin2rij;
            for3(i) {
                f[i] = dx[i] * (coeff * (1 - 1 / sin2rij));
            }
        }
        b.AddForce(f);
        a.SubForce(f);
    }
    

    EnsembleForce::EnsembleForce() {
        
    }
    
    EnsembleForce::~EnsembleForce() {
        
    }
    
    void EnsembleForce::ComputeForce(ParticleShapeArray& shapes) {
        ComputeMeanShape();
        for (int i = 0; i < shapes.size(); i++) {
            if (i > 0) {
                ParticleBSpline transform;
                transform.EstimateTransform(m_MeanShape, shapes[i]);
                FieldTransformType::Pointer fieldTransform = transform.GetTransform();
                shapes[i].TransformX2Y(fieldTransform.GetPointer());
                
            } else {
                shapes[i].TransformX2Y(NULL);
            }
        }
    }
}
