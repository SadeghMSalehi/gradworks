//
//  myImplicitSurfaceConstraint.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/23/12.
//
//

#include "myImplicitSurfaceConstraint.h"
#include "itkUnaryFunctorImageFilter.h"
#include "myImageContainer.h"
#include "itkVectorMagnitudeImageFilter.h"


typedef itk::VectorMagnitudeImageFilter<GradientImageType, SliceType> GradientMagnitudeFilterType;

template <class TIn, class TOut>
class InvertLabel {
public:
    InvertLabel() {}
    ~InvertLabel() {}
    bool operator!=(const InvertLabel&) const { return false; }
    bool operator==(const InvertLabel& other) const { return !(*this != other); }
    inline TOut operator()(const TIn& A) const { return !(A>0); }
};

template <class T1, class T2>
class BinaryThreshold {
    public :
    BinaryThreshold() {}
    ~BinaryThreshold() {}
    bool operator!=(const BinaryThreshold&) const { return false; }
    bool operator==(const BinaryThreshold& other) const { return !(*this != other); }
    inline T2 operator()(const T1& A) const { return (A>0)?1:0; }
};

void myImplicitSurfaceConstraint::Clear() {
    m_DistanceMaps.clear();
    m_DistanceMapInterpolators.clear();
    m_InsideDistanceVectorMaps.clear();
    m_OutsideDistanceVectorMaps.clear();
}

void myImplicitSurfaceConstraint::SetImageList(ImageContainer::List* imageList) {
    int nSubj = imageList->size();
    for (int i = 0; i < nSubj; i++) {
        LabelSliceType::Pointer labelMap = imageList->at(i)->GetLabelSlice();
        SignedDistanceMapFilterType::Pointer distmapFilter = SignedDistanceMapFilterType::New();
        distmapFilter->SetInput(labelMap);
        distmapFilter->Update();
        m_DistanceMaps.push_back(distmapFilter->GetOutput());
        m_OutsideDistanceVectorMaps.push_back(distmapFilter->GetVectorDistanceMap());

        typedef itk::UnaryFunctorImageFilter<LabelSliceType, LabelSliceType, InvertLabel<LabelSliceType::PixelType, LabelSliceType::PixelType> > InvertImageFilterType;
        InvertImageFilterType::Pointer invertFilter = InvertImageFilterType::New();
        invertFilter->SetInput(labelMap);
        invertFilter->Update();

        DistanceMapFilterType::Pointer insideDistmapFilter = DistanceMapFilterType::New();
        insideDistmapFilter->SetInput(invertFilter->GetOutput());
        insideDistmapFilter->Update();
        m_InsideDistanceVectorMaps.push_back(insideDistmapFilter->GetVectorDistanceMap());

        InterpolatorType::Pointer interpol = InterpolatorType::New();
        interpol->SetInputImage(m_DistanceMaps.back());
        m_DistanceMapInterpolators.push_back(interpol);

        typedef itk::UnaryFunctorImageFilter<LabelSliceType, LabelSliceType, BinaryThreshold<LabelSliceType::PixelType, LabelSliceType::PixelType> > BinaryLabelFilterType;

        BinaryLabelFilterType::Pointer binaryFilter = BinaryLabelFilterType::New();
        binaryFilter->SetInput(labelMap);
        binaryFilter->Update();

        GradientFilterType::Pointer gradient = GradientFilterType::New();
        gradient->SetInput(binaryFilter->GetOutput());
        gradient->SetSigma(.5);
        gradient->Update();
        m_GradientMaps.push_back(gradient->GetOutput());

        GradientInterpolatorType::Pointer gradientInterpolator = GradientInterpolatorType::New();
        gradientInterpolator->SetInputImage(gradient->GetOutput());
        m_GradientInterpolators.push_back(gradientInterpolator);

        GradientMagnitudeFilterType::Pointer magnitudeFilter = GradientMagnitudeFilterType::New();
        magnitudeFilter->SetInput(gradient->GetOutput());
        magnitudeFilter->Update();

        SliceToRGBAImageFilterType::Pointer rgbFilter = SliceToRGBAImageFilterType::New();
        rgbFilter->SetInput(magnitudeFilter->GetOutput());
        rgbFilter->SetUseInputImageExtremaForScaling(true);
        rgbFilter->Update();
        imageList->at(i)->AddDerivedView(imageList->at(i)->GetName() + "/gradientMagnitude", rgbFilter->GetOutput());

        SliceToRGBAImageFilterType::Pointer rgbFilter2 = SliceToRGBAImageFilterType::New();
        rgbFilter2->SetInput(m_DistanceMaps.at(i));
        rgbFilter2->SetUseInputImageExtremaForScaling(true);
        rgbFilter2->Update();
        imageList->at(i)->AddDerivedView(imageList->at(i)->GetName() + "/distanceMap", rgbFilter2->GetOutput());

    }
}

double myImplicitSurfaceConstraint::GetDistance(int subjId, SliceInterpolatorType::ContinuousIndexType &idx) const {
    return m_DistanceMapInterpolators[subjId]->EvaluateAtContinuousIndex(idx);
}

bool myImplicitSurfaceConstraint::IsInsideRegion(int subjId, SliceInterpolatorType::ContinuousIndexType& idx) const {
    if (subjId < m_DistanceMaps.size()) {
        return m_DistanceMapInterpolators[subjId]->IsInsideBuffer(subjId);
    }
    return false;
}

bool myImplicitSurfaceConstraint::IsInsideRegion(int subjId, SliceInterpolatorType::IndexType& idx) const {
    if (subjId < m_DistanceMaps.size()) {
        return m_DistanceMapInterpolators[subjId]->IsInsideBuffer(subjId);
    }
    return false;
}

myImplicitSurfaceConstraint::DistanceVectorImageType::PixelType myImplicitSurfaceConstraint::GetInsideOffset(int subjId, SliceType::IndexType& idx) const {
    return m_InsideDistanceVectorMaps[subjId]->GetPixel(idx);
}

myImplicitSurfaceConstraint::DistanceVectorImageType::PixelType myImplicitSurfaceConstraint::GetOutsideOffset(int subjId, SliceType::IndexType& idx) const {
    return m_OutsideDistanceVectorMaps[subjId]->GetPixel(idx);
}

myImplicitSurfaceConstraint::GradientPixelType myImplicitSurfaceConstraint::GetGradient(int subjId, GradientInterpolatorType::ContinuousIndexType& idx) const {
//    return m_GradientMaps[subjId]->GetPixel(idx);
    return m_GradientInterpolators[subjId]->EvaluateAtContinuousIndex(idx);
}

void myImplicitSurfaceConstraint::ApplyConstraint(OptimizerParametersType& params) const {
    const int nSubj = GetNumberOfSubjects();
    const int nDims = SliceType::GetImageDimension();
    const int nVars = params.GetSize() / nSubj;

    for (int n = 0; n < nSubj; n++) {
        for (int j = 0; j < nVars; j += nDims) {
            SliceType::IndexType idx;
            idx[0] = params[j];
            idx[1] = params[j+1];
            DistanceVectorImageType::PixelType offset = m_OutsideDistanceVectorMaps[n]->GetPixel(idx);
            params[j] += offset[j];
            params[j+1] += offset[j+1];
        }
    }
}