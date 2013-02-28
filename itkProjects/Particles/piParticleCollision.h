//
//  piParticleCollision.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#ifndef __ParticlesGUI__piParticleCollision__
#define __ParticlesGUI__piParticleCollision__

#include <iostream>

#include "piImageDef.h"
#include "piParticleCore.h"

namespace pi {

    static const int NO_CROSSING = 4;
    static const int STARTING_CONTACT = 1;
    static const int CROSSING_CONTACT = 2;
    static const int ENDING_CONTACT = 3;

    struct ContactPoint {
        DataReal cp[__Dim];
        int status;
    };

    class ParticleCollision {
    public:
        std::string binaryMaskCache;
        std::string distanceMapCache;
        bool applyMaskSmoothing;
        ParticleSubject* subject;

    public:
        ParticleCollision(): applyMaskSmoothing(false), subject(NULL) {}
        ~ParticleCollision() {}

        inline bool IsCrossing(IntIndex& xm) {
            return m_CrossingPicker->EvaluateAtIndex(xm) > 0;
        }

        inline bool IsRegionInside(IntIndex& xp) {
            return m_RegionPicker->EvaluateAtIndex(xp) > 0;
        }

        inline bool IsBufferInside(IntIndex& xp) {
            return m_RegionPicker->IsInsideBuffer(xp);
        }

        inline LabelImage::Pointer GetBinaryMask() {
            return m_BinaryMask;
        }

        bool ComputeContactPoint(DataReal* x0, DataReal* x1, ContactPoint& cp);
        bool ComputeNormal(DataReal* cp, DataReal* normal);
        DataReal ComputeDistance(DataReal* x1, DataReal* cp);
        void ComputeClosestBoundary(Particle& p, DataReal* x1, DataReal* cp);

        void SetLabelImage(LabelImage::Pointer labelImage);
        void SetBinaryMask(LabelImage::Pointer labelImage);

        void UpdateImages();

        bool LoadBinaryMask(std::string smoothingCache);
        bool LoadDistanceMap(std::string cache);

        void HandleCollision(ParticleSubject& subj);
        void HandleCollision(ParticleSubjectArray& subj);

        void ConstrainPoint(ParticleSubject& subj);
        void ProjectForceAndVelocity(ParticleSubject& subj);

    private:
        LabelImage::Pointer m_LabelImage;
        LabelImage::Pointer m_BinaryMask;
        LabelImage::Pointer m_ZeroCrossing;
        VectorImage::Pointer m_DistanceMap;
        GradientImage::Pointer m_Gradient;
        NNLabelInterpolatorType::Pointer m_CrossingPicker;
        NNLabelInterpolatorType::Pointer m_RegionPicker;
        GradientInterpolatorType::Pointer m_NormalPicker;
        NNVectorImageInterpolatorType::Pointer m_DistOffsetPicker;

    };
}
#endif /* defined(__ParticlesGUI__piParticleCollision__) */
