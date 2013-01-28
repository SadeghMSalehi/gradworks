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

namespace pi {
    class ParticleBSpline {
    public:
        ParticleBSpline() : m_SplineOrder(3), m_SplineLevel(1), m_ControlPoints(8) {
        }

        LabelImage::Pointer GetReferenceImage();
        void SetReferenceImage(LabelImage::Pointer img);
        void EstimateTransform(const ParticleSubject a, const ParticleSubject b);
        void ApplyTransform(ParticleSubject a);
        DoubleImage::Pointer WarpImage(DoubleImage::Pointer image);
        LabelImage::Pointer WarpLabel(LabelImage::Pointer srcImage);
        FieldTransformType::Pointer GetTransform();

    private:
        LabelImage::Pointer m_RefImage;
        DisplacementFieldPointSetType::Pointer m_FieldPoints;
        DisplacementFieldType::Pointer m_DisplacementField;

        int m_SplineOrder;
        int m_SplineLevel;
        int m_ControlPoints;
    };
}
#endif /* defined(__ParticlesGUI__myParticleBSpline__) */
