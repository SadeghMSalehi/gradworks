//
//  myParticleBSpline.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 1/25/13.
//
//

#ifndef __ParticlesGUI__myParticleBSpline__
#define __ParticlesGUI__myParticleBSpline__

#include <iostream>
#include "myParticleCore.h"
#include "itkDisplacementFieldTransform.h"

namespace my {


    class ParticleBSpline {
    public:
        ParticleBSpline() : m_SplineOrder(3), m_SplineLevel(1), m_ControlPoints(8) {
        }

        void EstimateTransform(const ParticleShape a, const ParticleShape b);
        void ApplyTransform(ParticleShape a);
        DoubleImage::Pointer WarpImage(DoubleImage::Pointer image);
        TransformType::Pointer GetTransform();

    private:
        DoubleImage::Pointer m_RefImage;
        DisplacementFieldPointSetType::Pointer m_FieldPoints;
        DisplacementFieldType::Pointer m_DisplacementField;

        int m_SplineOrder;
        int m_SplineLevel;
        int m_ControlPoints;
    };
}
#endif /* defined(__ParticlesGUI__myParticleBSpline__) */
