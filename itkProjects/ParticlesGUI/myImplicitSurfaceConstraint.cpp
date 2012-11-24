//
//  myImplicitSurfaceConstraint.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/23/12.
//
//

#include "myImplicitSurfaceConstraint.h"

void myImplicitSurfaceConstraint::AddSurface(SliceType::Pointer labelMap) {
    DistanceMapFilterType::Pointer distmapFilter = DistanceMapFilterType::New();
    distmapFilter->SetInput(labelMap);
    distmapFilter->Update();
    m_DistanceMaps.push_back(distmapFilter->GetOutput());
    m_DistanceVectorMaps.push_back(distmapFilter->GetVectorDistanceMap());

    SliceInterpolatorType::Pointer interpol = SliceInterpolatorType::New();
    interpol->SetInputImage(m_DistanceMaps.back());
    m_DistanceMapInterpolators.push_back(interpol);
}

double myImplicitSurfaceConstraint::GetDistance(int subjId, SliceInterpolatorType::ContinuousIndexType &idx) {
    return m_DistanceMapInterpolators[subjId]->EvaluateAtContinuousIndex(idx);
}

myImplicitSurfaceConstraint::DistanceVectorImageType::PixelType myImplicitSurfaceConstraint::GetDistanceVector(int subjId, SliceType::IndexType& idx) {
    return m_DistanceVectorMaps[subjId]->GetPixel(idx);
}

void myImplicitSurfaceConstraint::ApplyConstraint(OptimizerParametersType& params) {
    const int nSubj = GetNumberOfSubjects();
    const int nDims = SliceType::GetImageDimension();
    const int nVars = params.GetSize() / nSubj;

    for (int n = 0; n < nSubj; n++) {
        for (int j = 0; j < nVars; j += nDims) {
            SliceType::IndexType idx;
            idx[0] = params[j];
            idx[1] = params[j+1];
            DistanceVectorImageType::PixelType offset = m_DistanceVectorMaps[n]->GetPixel(idx);
            params[j] += offset[j];
            params[j+1] += offset[j+1];
        }
    }
}