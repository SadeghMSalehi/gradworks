//
//  ImplicitSurfaceConstraint.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/23/12.
//
//

#include "myImageDef.h"
#include "myImplicitSurfaceConstraint3D.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"


namespace my {
typedef itk::BinaryThresholdImageFilter<LabelImage, LabelImage> BinaryThresholdFilterType;
typedef itk::VectorMagnitudeImageFilter<VectorImage, LabelImage> GradientMagnitudeFilterType;

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

void ImplicitSurfaceConstraint::Clear() {
    m_DistanceMaps.clear();
    m_DistanceMapInterpolators.clear();
    m_InsideDistanceVectorMaps.clear();
    m_OutsideDistanceVectorMaps.clear();
    m_GradientInterpolators.clear();
    m_GradientMaps.clear();
}

void ImplicitSurfaceConstraint::SetImageList(LabelVectors& imageList) {
    int nSubj = imageList.size();
    Clear();
    for (int i = 0; i < nSubj; i++) {
        LabelImage::Pointer labelMap = imageList[i];

        // create binary image for a mask for a correct distance map
        BinaryThresholdFilterType::Pointer binThreshFilter = BinaryThresholdFilterType::New();
        binThreshFilter->SetInput(labelMap);
        binThreshFilter->SetInsideValue(1);
        binThreshFilter->SetOutsideValue(0);
        binThreshFilter->SetLowerThreshold(1);
        binThreshFilter->SetUpperThreshold(255);
        LabelImage::Pointer binaryMap = binThreshFilter->GetOutput();

        // construct signed distance filter
        SignedDistanceMapFilterType::Pointer distmapFilter = SignedDistanceMapFilterType::New();
        distmapFilter->SetInput(binaryMap);
        distmapFilter->Update();
        m_DistanceMaps.push_back(distmapFilter->GetOutput());
        m_OutsideDistanceVectorMaps.push_back(distmapFilter->GetVectorDistanceMap());

        // to compute inside offset correctly, invert the label map
        typedef itk::UnaryFunctorImageFilter<LabelImage, LabelImage, InvertLabel<LabelImage::PixelType, LabelImage::PixelType> > InvertImageFilterType;
        InvertImageFilterType::Pointer invertFilter = InvertImageFilterType::New();
        invertFilter->SetInput(binaryMap);
        invertFilter->Update();

        // compute inside offset using inverted image
        DistanceMapFilterType::Pointer insideDistmapFilter = DistanceMapFilterType::New();
        insideDistmapFilter->SetInput(invertFilter->GetOutput());
        insideDistmapFilter->Update();
        m_InsideDistanceVectorMaps.push_back(insideDistmapFilter->GetVectorDistanceMap());

        LinearImageInterpolatorType::Pointer interpol = LinearImageInterpolatorType::New();
        interpol->SetInputImage(m_DistanceMaps.back());
        m_DistanceMapInterpolators.push_back(interpol);

        // compute gradient for outside boundary
        GradientFilterType::Pointer gradient = GradientFilterType::New();
        gradient->SetInput(binaryMap);
        gradient->SetSigma(.5);
        gradient->Update();
        m_GradientMaps.push_back(gradient->GetOutput());

        GradientInterpolatorType::Pointer gradientInterpolator = GradientInterpolatorType::New();
        gradientInterpolator->SetInputImage(gradient->GetOutput());
        m_GradientInterpolators.push_back(gradientInterpolator);

        GradientMagnitudeFilterType::Pointer magnitudeFilter = GradientMagnitudeFilterType::New();
        magnitudeFilter->SetInput(gradient->GetOutput());
        magnitudeFilter->Update();
    }
}

double ImplicitSurfaceConstraint::GetDistance(int subjId, LinearImageInterpolatorType::ContinuousIndexType &idx) const {
    return m_DistanceMapInterpolators[subjId]->EvaluateAtContinuousIndex(idx);
}

bool ImplicitSurfaceConstraint::IsInsideRegion(int subjId, LinearImageInterpolatorType::ContinuousIndexType& idx) const {
    if (subjId < m_DistanceMaps.size()) {
        return m_DistanceMapInterpolators[subjId]->IsInsideBuffer(subjId);
    }
    return false;
}

bool ImplicitSurfaceConstraint::IsInsideRegion(int subjId, LinearImageInterpolatorType::IndexType& idx) const {
    if (subjId < m_DistanceMaps.size()) {
        return m_DistanceMapInterpolators[subjId]->IsInsideBuffer(subjId);
    }
    return false;
}

ImplicitSurfaceConstraint::DistanceVectorImageType::PixelType ImplicitSurfaceConstraint::GetInsideOffset(int subjId, DoubleImage::IndexType& idx) const {
    return m_InsideDistanceVectorMaps[subjId]->GetPixel(idx);
}

ImplicitSurfaceConstraint::DistanceVectorImageType::PixelType ImplicitSurfaceConstraint::GetOutsideOffset(int subjId, DoubleImage::IndexType& idx) const {
    return m_OutsideDistanceVectorMaps[subjId]->GetPixel(idx);
}

ImplicitSurfaceConstraint::GradientPixelType ImplicitSurfaceConstraint::GetGradient(int subjId, GradientInterpolatorType::ContinuousIndexType& idx) const {
//    return m_GradientMaps[subjId]->GetPixel(idx);
    return m_GradientInterpolators[subjId]->EvaluateAtContinuousIndex(idx);
}

    /*
void ImplicitSurfaceConstraint::ApplyConstraint(OptimizerParametersType& params) const {
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
     */

}