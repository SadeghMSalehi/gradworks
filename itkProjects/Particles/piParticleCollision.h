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

namespace pi {

    static const int NO_CROSSING = 4;
    static const int STARTING_CONTACT = 1;
    static const int CROSSING_CONTACT = 2;
    static const int ENDING_CONTACT = 3;

    struct ContactPoint {
        double cp[__Dim];
        int status;
    };

    class ParticleCollision {
    private:
        typedef NNLabelInterpolatorType::IndexType IntIndex;
        typedef LinearVectorImageInterpolatorType::ContinuousIndexType RealIndex;
        
    public:
        ParticleCollision(): m_ApplySmoothing(false) {}
        ~ParticleCollision() {}

        inline bool IsCrossing(IntIndex& xm) {
            return m_CrossingPicker->EvaluateAtIndex(xm) > 0;
        }

        inline bool IsRegionInside(IntIndex& xp) {
            return m_RegionPicker->EvaluateAtIndex(xp) > 0;
        }

        
        bool ComputeContactPoint(double* x0, double* x1, ContactPoint& cp);
        bool ComputeNormal(double* cp, double* normal);
        double ComputeDistance(double* x1, double* cp);
        void ComputeClosestBoundary(double* x1, double* cp);
        void SetBinaryMask(LabelImage::Pointer binary);
        void UpdateImages();
        
        bool LoadDistanceMap(const char* filename);
        void SaveDistanceMap(const char* filename);
        
        void Write(std::string b, std::string c);
        
    private:
        bool m_ApplySmoothing;
        LabelImage::Pointer m_BinaryMask;
        LabelImage::Pointer m_AntiAliasedMask;
        LabelImage::Pointer m_InvertedBinaryMask;
        LabelImage::Pointer m_ZeroCrossing;
        VectorImage::Pointer m_DistanceMap;
        VectorImage::Pointer m_Gradient;
        NNLabelInterpolatorType::Pointer m_CrossingPicker;
        NNLabelInterpolatorType::Pointer m_RegionPicker;
        LinearVectorImageInterpolatorType::Pointer m_NormalPicker;
        NNVectorImageInterpolatorType::Pointer m_DistOffsetPicker;

    };
}
#endif /* defined(__ParticlesGUI__piParticleCollision__) */
