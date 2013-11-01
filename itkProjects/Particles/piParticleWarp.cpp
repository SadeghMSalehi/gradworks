//
//  ParticleWarp.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/31/13.
//
//

#include "piParticleWarp.h"
#include <itkDisplacementFieldTransform.h>
#include <itkBSplineScatteredDataPointSetToImageFilter.h>
#include <itkBSplineTransform.h>
#include <itkWarpImageFilter.h>

namespace pi {

    const int __SplineOrder = 3;

    typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DisplacementFieldPointSetType, DisplacementFieldType> BSplineFilterType;
    typedef BSplineFilterType::WeightsContainerType WeightsContainerType;
    typedef itk::BSplineTransform<PointReal,__Dim, __SplineOrder> BSplineTransform;
    typedef itk::WarpImageFilter<RealImage, RealImage, DisplacementFieldType> WarpImageFilterType;
    typedef itk::WarpImageFilter<LabelImage, LabelImage, DisplacementFieldType> WarpLabelFilterType;

    void ParticleWarp::setParameters(pi::ConfigFile &config) {
        controlSpacing = config["particles.bspline-transform.control-point-spacing"];
    }

    void ParticleWarp::estimateBsplineWarp(Px::Vector &src, Px::Vector &dst) {
        if (reference.IsNull()) {
            cout << "Cannot estimate grid size without reference image!" << endl;
            return;
        }

        // how itk::PointSet manage the number of points?
        if (m_FieldPoints.IsNull()) {
            m_FieldPoints = DisplacementFieldPointSetType::New();
            m_FieldPoints->Initialize();
        }

        int n = src.size();
        for (int i = 0; i < n; i++) {
            IntPointSetType::PointType srcPoint;
            IntPointSetType::PointType dstPoint;
            fordim(j) {
                srcPoint[j] = src[i][j];
                dstPoint[j] = dst[i][j];
            }
            VectorType vector;
            fordim(j) {
                vector[j] = dst[i][j] - src[i][j];
            }

            // if there is no warp, check below
//            cout << srcPoint << " => " << dstPoint << " : " << vector << endl;


            m_FieldPoints->SetPoint(i, srcPoint);
            m_FieldPoints->SetPointData(i, vector);

        }

        LabelImage::SizeType imageSize = reference->GetBufferedRegion().GetSize();
        LabelImage::SpacingType imageSpacing = reference->GetSpacing();
        LabelImage::PointType imageOrigin = reference->GetOrigin();
        LabelImage::DirectionType imageDirection = reference->GetDirection();


        BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
        BSplineFilterType::ArrayType numControlPoints;

        // control point is given by its number, not in physical unit!
        int numOfLevels = 3;
        for (int i = 0; i < imageSize.GetSizeDimension(); i++) {
            numControlPoints[i] = imageSize[i] / controlSpacing + __SplineOrder;
        }

        try {
            // debug: reparameterized point component is outside
            bspliner->SetOrigin(imageOrigin);
            bspliner->SetSpacing(imageSpacing);
            bspliner->SetSize(imageSize);
            bspliner->SetDirection(imageDirection);
            bspliner->SetGenerateOutputImage(true);
            bspliner->SetNumberOfLevels(numOfLevels);
            bspliner->SetSplineOrder(__SplineOrder);
            bspliner->SetNumberOfControlPoints(numControlPoints);
            bspliner->SetInput(m_FieldPoints);

            /*
            if (m_UseWeights && n == m_Weights.size()) {
                BSplineFilterType::WeightsContainerType::Pointer weights = BSplineFilterType::WeightsContainerType::New();

                weights->Reserve(n);

                for (int i = 0; i < n; i++) {
                    // how to choose src or dst weight?
                    weights->SetElement(i, m_Weights[i]);
                }
                bspliner->SetPointWeights(weights.GetPointer());
            }
             */

            bspliner->Update();
            displacementField = bspliner->GetOutput();
        } catch (itk::ExceptionObject& e) {
            e.Print(std::cout);
        }
    }


    LabelImage::Pointer ParticleWarp::warpLabel(LabelImage::Pointer input) {
        if (displacementField.IsNull()) {
            return input;
        }
        WarpLabelFilterType::Pointer warpFilter = WarpLabelFilterType::New();
        warpFilter->SetInput(input);
        warpFilter->SetDisplacementField(displacementField);
        warpFilter->SetOutputParametersFromImage(input);
        warpFilter->Update();
        return warpFilter->GetOutput();
    }

    RealImage::Pointer ParticleWarp::warpImage(RealImage::Pointer input) {
        if (displacementField.IsNull()) {
            return input;
        }
        WarpImageFilterType::Pointer warpFilter = WarpImageFilterType::New();
        warpFilter->SetInput(input);
        warpFilter->SetDisplacementField(displacementField);
        warpFilter->SetOutputParametersFromImage(input);
        warpFilter->Update();
        return warpFilter->GetOutput();
    }

}