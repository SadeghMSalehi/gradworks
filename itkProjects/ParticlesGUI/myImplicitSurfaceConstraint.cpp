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
#include "itkBinaryThresholdImageFilter.h"


typedef itk::BinaryThresholdImageFilter<LabelSliceType, LabelSliceType> BinaryThresholdFilterType;
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
    m_GradientInterpolators.clear();
    m_GradientMaps.clear();
}

void myImplicitSurfaceConstraint::SetImageList(ImageContainer::List* imageList) {
    int nSubj = imageList->size();
    Clear();
    for (int i = 0; i < nSubj; i++) {
        // every image must have corresponding label images
        if (!imageList->at(i)->HasLabel()) {
            return;
        }

        LabelSliceType::Pointer labelMap = imageList->at(i)->GetLabelSlice();

        // create binary image for a mask for a correct distance map
        BinaryThresholdFilterType::Pointer binThreshFilter = BinaryThresholdFilterType::New();
        binThreshFilter->SetInput(labelMap);
        binThreshFilter->SetInsideValue(1);
        binThreshFilter->SetOutsideValue(0);
        binThreshFilter->SetLowerThreshold(1);
        binThreshFilter->SetUpperThreshold(255);
        LabelSliceType::Pointer binaryMap = binThreshFilter->GetOutput();

        // construct signed distance filter
        SignedDistanceMapFilterType::Pointer distmapFilter = SignedDistanceMapFilterType::New();
        distmapFilter->SetInput(binaryMap);
        distmapFilter->Update();
        m_DistanceMaps.push_back(distmapFilter->GetOutput());
        m_OutsideDistanceVectorMaps.push_back(distmapFilter->GetVectorDistanceMap());

        // to compute inside offset correctly, invert the label map
        typedef itk::UnaryFunctorImageFilter<LabelSliceType, LabelSliceType, InvertLabel<LabelSliceType::PixelType, LabelSliceType::PixelType> > InvertImageFilterType;
        InvertImageFilterType::Pointer invertFilter = InvertImageFilterType::New();
        invertFilter->SetInput(binaryMap);
        invertFilter->Update();

        // compute inside offset using inverted image
        DistanceMapFilterType::Pointer insideDistmapFilter = DistanceMapFilterType::New();
        insideDistmapFilter->SetInput(invertFilter->GetOutput());
        insideDistmapFilter->Update();
        m_InsideDistanceVectorMaps.push_back(insideDistmapFilter->GetVectorDistanceMap());

        InterpolatorType::Pointer interpol = InterpolatorType::New();
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

        // prepare gradient magnitude attribute image
        SliceToRGBAImageFilterType::Pointer rgbFilter = SliceToRGBAImageFilterType::New();
        rgbFilter->SetInput(magnitudeFilter->GetOutput());
        rgbFilter->SetUseInputImageExtremaForScaling(true);
        rgbFilter->Update();
        imageList->at(i)->AddDerivedView(imageList->at(i)->GetName() + "/gradientMagnitude", rgbFilter->GetOutput());

        // prepare distance map as an attribute image
        SliceToRGBAImageFilterType::Pointer rgbFilter2 = SliceToRGBAImageFilterType::New();
        rgbFilter2->SetInput(m_DistanceMaps.at(i));
        rgbFilter2->SetUseInputImageExtremaForScaling(true);
        rgbFilter2->Update();
        imageList->at(i)->AddDerivedView(imageList->at(i)->GetName() + "/distanceMap", rgbFilter2->GetOutput());
    }

    std::cout << "surface constraint updated: " << imageList->at(0)->GetSliceIndex() << std::endl;
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