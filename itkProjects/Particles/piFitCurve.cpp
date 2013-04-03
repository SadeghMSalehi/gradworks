//
//  piFitCurve.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/3/13.
//
//

#include "piFitCurve.h"

#include <cmath>

namespace pi {
    CurveFitting::CurveFitting() {

    }

    CurveFitting::~CurveFitting() {

    }

    ParticleVector& CurveFitting::GetResult() {
        return _result;
    }
    
    void CurveFitting::FitCurve(pi::ParticleVector &particles, bool closed) {
        ComputeChordLengthParameterization(particles);
        ApproximateCurve(particles, closed);
    }

    CurveFitting::DenseScalarImage::Pointer CurveFitting::OneDimensionalBsplineFitting(pi::ParticleVector &particles, int dim, bool closed) {
        SparseScalarSet::Pointer scalarSet = SparseScalarSet::New();
        const int nParticles = particles.size();
        for (int i = 0; i < nParticles; i++) {
            scalarSet->SetPoint(i, _chordLength[i]);
            ScalarPoint point;
            scalarSet->SetPointData(i, particles[i].x[dim]);
        }

        BSplineFitter::ArrayType closedDimensions;
        closedDimensions[0] = 1;

        BSplineFitter::PointType origin;
        origin[0] = 0;

        BSplineFitter::SpacingType spacing;
        spacing[0] = 0.01;

        BSplineFitter::SizeType size;
        size[0] = 1;

        BSplineFitter::ArrayType controlPoints;
        controlPoints[0] = 25;

        BSplineFitter::Pointer fitter = BSplineFitter::New();
        fitter->SetCloseDimension(closedDimensions);
        fitter->SetOrigin(origin);
        fitter->SetSpacing(spacing);
        fitter->SetSize(size);
        fitter->SetNumberOfControlPoints(controlPoints);
        fitter->SetNumberOfLevels(1);
        fitter->SetSplineOrder(3);
        fitter->SetGenerateOutputImage(true);
        fitter->SetInput(scalarSet);
        fitter->Update();

        DenseScalarImage::Pointer result = fitter->GetOutput();
        return result;
    }

    void CurveFitting::ApproximateCurve(ParticleVector &particles, bool closed) {
        assert(_chordLength.size() == particles.size());

        DenseScalarImage::Pointer xFit = OneDimensionalBsplineFitting(particles, 0, closed);
        DenseScalarImage::Pointer yFit = OneDimensionalBsplineFitting(particles, 1, closed);

        const int nElems = xFit->GetBufferedRegion().GetSize(0);
        _result.resize(nElems);

        ScalarPoint* xSrc = xFit->GetBufferPointer();
        ScalarPoint* ySrc = yFit->GetBufferPointer();
        for (int i = 0; i < nElems; i++) {
            _result[i].x[0] = xSrc[i][0];
            _result[i].x[1] = ySrc[i][0];
        }
    }

    void CurveFitting::ComputeChordLengthParameterization(pi::ParticleVector &particles) {
        const int nParticles = particles.size();
        _chordLength.set_size(nParticles);
        ParticleVector::iterator iter = particles.begin() + 1;
        for (int i = 0; iter != particles.end(); i++, iter++) {
            Particle& prev = *(iter-1);
            Particle& curr = *(iter);
            double dist = std::sqrt(curr.Dist2(prev));
            _chordLength[i] = dist;
        }
        _chordLength /= _chordLength.sum();
    }
}