//
//  myEnsembleEntropy.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/24/12.
//
//

#include "myEnsembleEntropy.h"


myEnsembleEntropy::myEnsembleEntropy() {
    
}

myEnsembleEntropy::~myEnsembleEntropy() {
    
}

void myEnsembleEntropy::SetImageList(ImageContainer::List *imageList) {
    m_ImageList = imageList;
}

/**
 * must provide exactly matching number of subjects and points with the already provided image list
 */
void myEnsembleEntropy::SetInitialPositions(const OptimizerParametersType &params, int nSubject, int nPoints, int nParams) {
    m_nSubjects = nSubject;
    m_nPoints = nPoints;
    m_nParams = nParams;
}

void myEnsembleEntropy::EstimateRigidParameters(MatrixType &transformParams, const OptimizerParametersType &params, int target, int source) const {
    const int Dim = 2;
    
    MatrixType targetPoints(&params[target * m_nParams], m_nPoints, Dim);
    MatrixType sourcePoints(&params[source * m_nParams], m_nPoints, Dim);
    
    // make those centroid match
    VectorType targetCentroid(2), sourceCentroid(2);
    targetCentroid.fill(0);
    sourceCentroid.fill(0);
    
    for (int i = 0; i < m_nPoints; i++) {
        for (int k = 0; k < Dim; k++) {
            targetCentroid[k] += targetPoints[i][k];
            sourceCentroid[k] += sourcePoints[i][k];
        }
    }
    targetCentroid /= m_nPoints;
    sourceCentroid /= m_nPoints;

    for (int i = 0; i < m_nPoints; i++) {
        for (int k = 0; k < Dim; k++) {
            targetPoints[i][k] = targetPoints[i][k] - targetCentroid[k];
            sourcePoints[i][k] = sourcePoints[i][k] - sourceCentroid[k];
        }
    }
    
    double rotationParam;
    VectorType translationParam = targetCentroid - sourceCentroid;
    
    MatrixType cov = targetPoints.transpose() * sourcePoints;
    vnl_svd<double> svd(cov);
    MatrixType rotationMatrix = svd.V() * svd.U().transpose();

    rotationParam = acos(rotationMatrix[0][0]);
    
    transformParams[source][0] = translationParam[0];
    transformParams[source][1] = translationParam[1];
    transformParams[source][2] = rotationParam;
}

void myEnsembleEntropy::ComputePositionalEnsemble(const OptimizerParametersType &params, double &cost, MatrixType &deriv) const {
    
    MatrixType transformParams(m_nSubjects, m_nTransformParams);
    transformParams.fill(0);
    
    for (int i = 1; i < m_nSubjects; i++) {
        EstimateRigidParameters(transformParams, params, 0, i);
    }

    const int Dim = 2;
    MatrixType targetPoints(&params[0], m_nPoints, Dim);
    for (int n = 1; n < m_nSubjects; n++) {
        MatrixType rotationMatrix(2,2);
        rotationMatrix[0][0] = ::cos(transformParams[n][2]);
        rotationMatrix[0][1] = ::sin(transformParams[n][2]);
        rotationMatrix[1][0] = -::sin(transformParams[n][2]);
        rotationMatrix[1][1] = ::cos(transformParams[n][2]);

        MatrixType sourcePoints(&params[n * m_nParams], m_nPoints, Dim);
        MatrixType transformedPoints = sourcePoints * rotationMatrix;
        for (int i = 0; i < m_nPoints; i++) {
            transformedPoints[i][0] -= transformParams[n][0];
            transformedPoints[i][1] -= transformParams[n][1];
        }
        
    }
}