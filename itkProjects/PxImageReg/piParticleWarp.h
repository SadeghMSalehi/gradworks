//
//  ParticleWarp.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/31/13.
//
//

#ifndef __ParticleGuidedRegistration__ParticleWarp__
#define __ParticleGuidedRegistration__ParticleWarp__

#include <iostream>
#include "piParticleRunner.h"

namespace pi {
    class ParticleMesh {
    public:
        void constructNeighbors(int regionId, int nPx, PxSubj& subj, double cutoff, PxGlobal::Neighbors& neighbors);
    };


    /// @brief A class estimates the warp between two set of corresponding points. The set of displacement vectors is constructed from sets of points. The reference label image and the control point spacing must be set before calling ParticleWarp::estimateBsplineWarp.
    class ParticleWarp {
    public:
        /// The control point spacing of B-spline grid in pixels
        int controlSpacing;

        /// The reference image to define the displacement field as well as the B-spline grid. __This should be set before calling ParticleWarp::estimateBsplineWarp__.
        LabelImage::Pointer reference;

        /// The resulting displacement fields
        DisplacementFieldType::Pointer displacementField;

        ParticleWarp(): controlSpacing(4) {}

        /// @brief Estimate the B-spline deformation fields. The number of points should be the same and corresponding between the source and the destination point set.
        /// @param src The source point set
        /// @param dst The destination point set
        void estimateBsplineWarp(Px::Vector& src, Px::Vector& dst);


        /// @brief Warp the given label image with the deformation estimated
        /// @param input LabelImage type
        LabelImage::Pointer warpLabel(LabelImage::Pointer input);

        /// @brief Warp the given intensity image with the deformation estimated
        /// @param input RealImage type
        RealImage::Pointer warpImage(RealImage::Pointer input);

    private:
        DisplacementFieldPointSetType::Pointer m_FieldPoints;
        FieldTransformType::Pointer transform;
    };
}
#endif /* defined(__ParticleGuidedRegistration__ParticleWarp__) */
