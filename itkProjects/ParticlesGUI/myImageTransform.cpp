//
//  myImageTransform.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 12/3/12.
//
//

#include "myImageTransform.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "itkElasticBodySplineKernelTransform.h"
#include "itkResampleImageFilter.h"

const static int nDim = 2;

typedef itk::ThinPlateR2LogRSplineKernelTransform<double, nDim> R2LogRTPSTransformType;
typedef itk::ThinPlateSplineKernelTransform<double, nDim> TPSTransformType;
typedef itk::ElasticBodySplineKernelTransform<double, nDim> EBSTransformType;
typedef itk::ResampleImageFilter<SliceType, SliceType> ResampleImageFilterType;

myImageTransform::myImageTransform() {

}

myImageTransform::~myImageTransform() {

}


myImageTransform::KernelTransformPointer myImageTransform::CreateKernelTransform(int type, int n, double* src, double* dst) {
    int nPoints = n;
    TPSTransformType::PointsContainer::Pointer targetPoints = TPSTransformType::PointsContainer::New();

    double* pSrc = src;
    for (int i = 0; i < nPoints; i++) {
        TPSTransformType::PointsContainer::Element iPoint;
        iPoint[0] = pSrc[0];
        iPoint[1] = pSrc[1];
        targetPoints->InsertElement(i, iPoint);
        pSrc += nDim;
    }
    TPSTransformType::PointSetPointer targetPointSet = TPSTransformType::PointSetType::New();
    targetPointSet->SetPoints(targetPoints.GetPointer());

    double *pDst = dst;
    KernelTransformType::PointsContainer::Pointer sourcePoints = KernelTransformType::PointsContainer::New();
    for (int i = 0; i < nPoints; i++) {
        KernelTransformType::PointsContainer::Element iPoint;
        iPoint[0] = pDst[0];
        iPoint[1] = pDst[1];
        sourcePoints->InsertElement(i, iPoint);
        pDst += nDim;
    }
    KernelTransformType::PointSetPointer srcPointSet = TPSTransformType::PointSetType::New();
    srcPointSet->SetPoints(sourcePoints.GetPointer());

    // set up TPS transform
    KernelTransformType::Pointer kernelTransform(NULL);
    if (type == 0) {
        TPSTransformType::Pointer tps = TPSTransformType::New();
        tps->SetSourceLandmarks(srcPointSet.GetPointer());
        tps->SetTargetLandmarks(targetPointSet.GetPointer());
        tps->ComputeWMatrix();
        kernelTransform = tps;
    } else if (type == 1) {
        EBSTransformType::Pointer ebs = EBSTransformType::New();
        ebs->SetSourceLandmarks(srcPointSet.GetPointer());
        ebs->SetTargetLandmarks(targetPointSet.GetPointer());
        ebs->ComputeWMatrix();
        kernelTransform = ebs;
    } else if (type == 2) {
        R2LogRTPSTransformType::Pointer tps = R2LogRTPSTransformType::New();
        tps->SetSourceLandmarks(srcPointSet.GetPointer());
        tps->SetTargetLandmarks(targetPointSet.GetPointer());
        tps->ComputeWMatrix();
    }
    return kernelTransform;
}

myImageTransform::BSplineTransform::Pointer myImageTransform::CreateBSplineTransform(int type, int n, double *src, double* dst, SliceType::Pointer spaceImage) {
    BSplineTransform::Pointer transform = BSplineTransform::New();
    transform->SetGridRegion(spaceImage->GetBufferedRegion());
    transform->SetGridSpacing(spaceImage->GetSpacing());
    transform->SetGridOrigin(spaceImage->GetOrigin());
    
    BSplineTransform::ParametersType params(transform->GetNumberOfParameters());
    transform->SetParameters(params);

    return transform;
}

SliceType::Pointer myImageTransform::ResampleSlice(SliceType::Pointer sourceImage, myImageTransform::KernelTransformPointer  transform) {
    ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
    resampler->SetInput(sourceImage);
    resampler->SetTransform(transform.GetPointer());
    resampler->SetReferenceImage(sourceImage);
    resampler->UseReferenceImageOn();
    resampler->Update();
    return resampler->GetOutput();
}
