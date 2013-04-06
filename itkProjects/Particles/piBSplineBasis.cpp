//
//  piBSplineBasis.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/4/13.
//
//

#include "piBSplineBasis.h"


namespace pi {
    static void ComputeChordLength(pi::ParticleVector &particles) {
        const int nParticles = particles.size();
        ParticleVector::iterator iter = particles.begin() + 1;
        for (int i = 1; iter != particles.end(); i++, iter++) {
            Particle& prev = *(iter-1);
            Particle& curr = *(iter);
            double dist = std::sqrt(curr.Dist2(prev));
            particles[i].t = particles[i-1].t + dist;
        }

        double sum = particles.back().t + std::sqrt(particles.back().Dist2(particles.front()));
        for (int i = 0; i < nParticles; i++) {
            particles[i].t /= sum;
        }
    }

    void BSplineBasis::CubicContour(ParticleVector& controls, std::vector<TReal>& params, ParticleVector& contourOut) {
        const int nControls = controls.size();
        vnl_matrix<TReal> cubicBasis(1,4);
        vnl_matrix<TReal> controlPos(4,2);

        // reparameterize with number of controls
        vnl_vector<TReal> t(params.size());
        vnl_vector<int> s(params.size());

        contourOut.resize(t.size());
        for (int i = 0; i < t.size(); i++) {
            double tt = controls.size() * params[i];
            s[i] = int(tt);
            t[i] = tt - s[i];

            for (int j = 0; j < 4; j++) {
                controlPos[j][0] = controls[(s[i]+j)%nControls].x[0];
                controlPos[j][1] = controls[(s[i]+j)%nControls].x[1];
            }
            BSplineBasis::CubicBasis(t[i], cubicBasis.data_block());
            vnl_matrix<TReal> contourPos = cubicBasis * controlPos;

            contourOut[i].x[0] = contourPos[0][0];
            contourOut[i].x[1] = contourPos[0][1];
        }
    }

    void BSplineBasis::CubicContourFitting(ParticleVector &particles, int nControls, ParticleVector& controlsOut) {
        assert(particles.size() >= nControls);
        ComputeChordLength(particles);

        const int nPoints = particles.size();
        vnl_vector<TReal> cubicBasis(4);
        vnl_matrix<TReal> basisMatrix(nPoints, nControls);
        vnl_matrix<TReal> pos(nPoints, 2);

        basisMatrix.fill(0);
        for (int i = 0; i < particles.size(); i++) {
            Particle& p = particles[i];
            pos[i][0] = p.x[0];
            pos[i][1] = p.x[1];

            double s = p.t * nControls;
            p.label = int(s);
            p.t = s - p.label;
            BSplineBasis::CubicBasis(p.t, cubicBasis.data_block());

            for (int j = 0; j < 4; j++) {
                int offset = (p.label + j) % nControls;
                basisMatrix[i][offset] = cubicBasis[j];
            }
        }
        vnl_matrix<TReal> invBasis = vnl_matrix_inverse<TReal>(basisMatrix);
        vnl_matrix<TReal> controlsSolution = invBasis * pos;

        
        controlsOut.resize(nControls);
        for (int i = 0; i < nControls; i++) {
            controlsOut[i].x[0] = controlsSolution[i][0];
            controlsOut[i].x[1] = controlsSolution[i][1];
        }
    }
}