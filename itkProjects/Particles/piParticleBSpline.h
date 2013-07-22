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

#ifndef NCONTROL
#define NCONTROL 16
#endif

#ifndef NORDER
#define NORDER 3
#endif

namespace pi {
    typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DisplacementFieldPointSetType, DisplacementFieldType> BSplineFilterType;
    typedef BSplineFilterType::WeightsContainerType WeightsContainerType;
    typedef itk::BSplineTransform<PointReal,__Dim,3> BSplineTransform;
    typedef itk::WarpImageFilter<RealImage, RealImage, DisplacementFieldType> WarpImageFilterType;

    static int __SPLINE_ORDER__ = NORDER;
    static int __SPLINE_CONTROL_POINTS__ = NCONTROL;

    class ParticleBSpline {
    public:
        ParticleBSpline() : m_SplineOrder(__SPLINE_ORDER__), m_SplineLevel(1), m_ControlPoints(__SPLINE_CONTROL_POINTS__) {
            m_ControlPointSpacing = 4;
        }

        int m_ControlPointSpacing;
        int m_SplineOrder;
        int m_SplineLevel;
        int m_ControlPoints;

        void SetControlPointSpacing(int spacing) {
            m_ControlPointSpacing = spacing;
        }

        void SetWeights(VNLVector weights) {
            m_UseWeights = true;
            m_Weights = weights;
        }

        LabelImage::Pointer GetReferenceImage();
        void SetReferenceImage(LabelImage::Pointer img);
        void EstimateTransform(const ParticleSubject& a, const ParticleSubject& b);
        void EstimateTransformY(const ParticleSubject& a, const ParticleSubject& b);
        void EstimateTransformZ(const ParticleSubject& a, const ParticleSubject& b);
        void EstimateTransformYZ(const ParticleSubject& a, const ParticleSubject& b);



        template <class C, class R, class T> void EstimateTransform(const T& src, const T& dst, const int nPoints, typename R::Pointer refImage);
        void ApplyTransform(ParticleSubject& a);
        RealImage::Pointer WarpImage(RealImage::Pointer image);
        LabelImage::Pointer WarpLabel(LabelImage::Pointer srcImage);
        FieldTransformType::Pointer GetTransform();

    private:
        LabelImage::Pointer m_RefImage;
        DisplacementFieldPointSetType::Pointer m_FieldPoints;
        DisplacementFieldType::Pointer m_DisplacementField;
        bool m_UseWeights;
        VNLVector m_Weights;

    };

    template <class C,class R,class T>
    void ParticleBSpline::EstimateTransform(const T& src, const T& dst, const int nPoints, typename R::Pointer refImage) {
        if (refImage.IsNull()) {
            cout << "Cannot estimate grid size without reference image!" << endl;
            return;
        }
        if (m_FieldPoints.IsNull()) {
            m_FieldPoints = DisplacementFieldPointSetType::New();
        }
        m_FieldPoints->Initialize();


        // casting object C;
        C caster;
        int n = nPoints;
        for (int i = 0; i < n; i++) {
            IntPointSetType::PointType srcPoint;
            IntPointSetType::PointType dstPoint;
            fordim(j) {
                srcPoint[j] = caster.castSource(src[i],j);
                dstPoint[j] = caster.castTarget(dst[i],j);
            }
            VectorType vector;
            fordim(j) {
                vector[j] = caster.castTarget(dst[i],j) - caster.castSource(src[i],j);
            }
            m_FieldPoints->SetPoint(i, srcPoint);
            m_FieldPoints->SetPointData(i, vector);

//            cout << "src[" << i << "] = " << srcPoint[0] << "," << srcPoint[1] << "," << srcPoint[2] << endl;
//            cout << "dst[" << i << "] = " << dstPoint[0] << "," << dstPoint[1] << "," << dstPoint[2] << endl;
        }

        int splineOrder = m_SplineOrder;
        int numOfLevels = m_SplineLevel;
        int nSize = m_ControlPoints;


        BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
        BSplineFilterType::ArrayType numControlPoints;
        numControlPoints.Fill(nSize + splineOrder);

        typename R::SizeType imageSize = refImage->GetBufferedRegion().GetSize();
        typename R::SpacingType imageSpacing = refImage->GetSpacing();
        typename R::PointType imageOrigin = refImage->GetOrigin();
        typename R::DirectionType imageDirection = refImage->GetDirection();



        for (int i = 0; i < imageSize.GetSizeDimension(); i++) {
            numControlPoints[i] = imageSize[i] / m_ControlPointSpacing + splineOrder;
        }

//        cout << "# control points: " << numControlPoints << endl;

        try {
            // debug: reparameterized point component is outside
            bspliner->SetOrigin(imageOrigin);
            bspliner->SetSpacing(imageSpacing);
            bspliner->SetSize(imageSize);
            bspliner->SetDirection(imageDirection);
            bspliner->SetGenerateOutputImage(true);
            bspliner->SetNumberOfLevels(numOfLevels);
            bspliner->SetSplineOrder(splineOrder);
            bspliner->SetNumberOfControlPoints(numControlPoints);
            bspliner->SetInput(m_FieldPoints);

            if (m_UseWeights && n == m_Weights.size()) {
                BSplineFilterType::WeightsContainerType::Pointer weights = BSplineFilterType::WeightsContainerType::New();

                weights->Reserve(n);

                for (int i = 0; i < n; i++) {
                    // how to choose src or dst weight?
                    weights->SetElement(i, m_Weights[i]);
                }
                bspliner->SetPointWeights(weights.GetPointer());
            }
            
            bspliner->Update();
            m_DisplacementField = bspliner->GetOutput();
        } catch (itk::ExceptionObject& e) {
            e.Print(std::cout);
        }
    }
    

}
#endif /* defined(__ParticlesGUI__myParticleBSpline__) */
