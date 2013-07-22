//
//  piParticleTrainer.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 7/21/13.
//
//

#ifndef __ParticleGuidedRegistration__piParticleTrainer__
#define __ParticleGuidedRegistration__piParticleTrainer__

#include <iostream>
#include <vector>
#include "piParticleCore.h"

namespace pi {

class ParticleTrainer {
public:
    ParticleTrainer(ParticleSystem* system);
    ~ParticleTrainer();

    void setNumberOfParticles(IntVector& counts);
    void trainParticles();

private:
    void sampleParticles();
    void spreadParticles(int labelId);
    LabelImage::Pointer extractLabel(LabelImage::Pointer labelImage, int labelId);

private:
    ParticleSystem* _system;
    int _numLabels;
    int _numSubjs;
    IntVector _counts;
    
    typedef std::vector<LabelImage::PointType> PointVector;
    typedef std::vector<PointVector> LabelPoints;

};

}
#endif /* defined(__ParticleGuidedRegistration__piParticleTrainer__) */
