//
//  gemFitCurves.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/3/13.
//
//

#ifndef ParticleGuidedRegistration_gemFitCurves_h
#define ParticleGuidedRegistration_gemFitCurves_h

namespace gem {
    int FitCurve(double* d, int nPts, double error);
    void CopyResults(double* vector);
}

#endif
