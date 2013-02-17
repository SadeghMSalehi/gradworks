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
#include "piParticleCore.h"
#include "itkDisplacementFieldTransform.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkBSplineTransform.h"
#include "itkWarpImageFilter.h"

namespace pi {
    typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DisplacementFieldPointSetType, DisplacementFieldType> BSplineFilterType;
    typedef BSplineFilterType::WeightsContainerType WeightsContainerType;
    typedef itk::BSplineTransform<double,__Dim,3> BSplineTransform;
    typedef itk::WarpImageFilter<DoubleImage, DoubleImage, DisplacementFieldType> WarpImageFilterType;


    class ParticleBSpline {
    public:
        ParticleBSpline() : m_SplineOrder(3), m_SplineLevel(1), m_ControlPoints(8) {
        }

        LabelImage::Pointer GetReferenceImage();
        void SetReferenceImage(LabelImage::Pointer img);
        void EstimateTransform(const ParticleSubject& a, const ParticleSubject& b);
        template <class T> void EstimateTransform(const T& src, const T& dst, const int nPoints);
        void ApplyTransform(ParticleSubject& a);
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

    template <class T>
    void ParticleBSpline::EstimateTransform(const T& src, const T& dst, const int nPoints) {
        if (m_FieldPoints.IsNull()) {
            m_FieldPoints = DisplacementFieldPointSetType::New();
        }
        m_FieldPoints->Initialize();

        // create point structures
        IntPointSetType::Pointer srcPoints = IntPointSetType::New();
        IntPointSetType::Pointer dstPoints = IntPointSetType::New();

        srcPoints->Initialize();
        dstPoints->Initialize();

        int n = nPoints;
        for (int i = 0; i < n; i++) {
            IntPointSetType::PointType iPoint;
            fordim(j) {
                iPoint[j] = src[i].x[j];
            }
            VectorType vector;
            fordim(j) {
                vector[j] = dst[i].x[j] - src[i].x[j];
            }
            m_FieldPoints->SetPoint(i, iPoint);
            m_FieldPoints->SetPointData(i, vector);
        }

        int splineOrder = m_SplineOrder;
        int numOfLevels = m_SplineLevel;
        int nSize = m_ControlPoints;

        BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
        BSplineFilterType::ArrayType numControlPoints;
        numControlPoints.Fill(nSize + splineOrder);

        DoubleImage::SizeType imageSize = m_RefImage->GetBufferedRegion().GetSize();
        DoubleImage::SpacingType imageSpacing = m_RefImage->GetSpacing();
        DoubleImage::PointType imageOrigin = m_RefImage->GetOrigin();

        try {
            // debug: reparameterized point component is outside
            bspliner->SetOrigin(imageOrigin);
            bspliner->SetSpacing(imageSpacing);
            bspliner->SetSize(imageSize);
            bspliner->SetGenerateOutputImage(true);
            bspliner->SetNumberOfLevels(numOfLevels);
            bspliner->SetSplineOrder(splineOrder);
            bspliner->SetNumberOfControlPoints(numControlPoints);
            bspliner->SetInput(m_FieldPoints);
            bspliner->Update();
            m_DisplacementField = bspliner->GetOutput();
        } catch (itk::ExceptionObject& e) {
            e.Print(std::cout);
        }
    }
    

}
#endif /* defined(__ParticlesGUI__myParticleBSpline__) */
