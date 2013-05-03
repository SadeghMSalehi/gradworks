//
//  piFitCurve.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/3/13.
//
//

#include "piFitCurve.h"
#include "piBSplineBasis.h"
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
        ParticleVector controls;
        BSplineBasis::CubicContourFitting(particles, 25, controls);

        std::vector<float> params(50);
        for (int i = 0; i < params.size(); i++) {
            params[i] = i / float(params.size());
        }
        BSplineBasis::CubicContour(controls, params, _result);
    }

    CurveFitting::DenseScalarImage::Pointer CurveFitting::OneDimensionalBsplineFitting(pi::ParticleVector &particles, int dim, bool closed) {
        SparseScalarSet::Pointer scalarSet = SparseScalarSet::New();
        const int nParticles = particles.size();
        for (int i = 0; i < nParticles; i++) {
            SparseScalarSet::PointType point;
            point[0] = particles[i].t;
            scalarSet->SetPoint(i, point);
            ScalarPoint data;
            data[0] = particles[i].x[dim];
            scalarSet->SetPointData(i, data);
            cout << particles[i].t << "," << particles[i].x[dim] << endl;
        }

        BSplineFitter::ArrayType closedDimensions;
        closedDimensions[0] = 1;

        BSplineFitter::PointType origin;
        origin[0] = -1;

        BSplineFitter::SpacingType spacing;
        spacing[0] = 0.1;

        BSplineFitter::SizeType size;
        size[0] = 30;

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
        DenseScalarImage::Pointer xFit = OneDimensionalBsplineFitting(particles, 0, closed);
        DenseScalarImage::Pointer yFit = OneDimensionalBsplineFitting(particles, 1, closed);


        BSplineInterpolator::Pointer xIntp = BSplineInterpolator::New();
        xIntp->SetInputImage(xFit);

        BSplineInterpolator::Pointer yIntp = BSplineInterpolator::New();
        yIntp->SetInputImage(yFit);

        BSplineInterpolator::PointType point;
        pi::createParticles(_result, 0, particles.size());
        for (int i = 0; i < particles.size(); i++) {
            point[0] = particles[i].t;
            _result[i].x[0] = xIntp->Evaluate(point)[0];
            _result[i].x[1] = yIntp->Evaluate(point)[0];
        }
    }

    void CurveFitting::ComputeChordLengthParameterization(pi::ParticleVector &particles) {
        const int nParticles = particles.size();
        ParticleVector::iterator iter = particles.begin() + 1;
        for (int i = 1; iter != particles.end(); i++, iter++) {
            Particle& prev = *(iter-1);
            Particle& curr = *(iter);
            double dist = std::sqrt(curr.Dist2(prev));
            particles[i].t = particles[i-1].t + dist;
        }

        for (int i = 0; i < nParticles; i++) {
            particles[i].t /= particles.back().t;
        }
    }
}