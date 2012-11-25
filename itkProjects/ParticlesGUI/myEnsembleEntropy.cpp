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
    
    const MatrixType targetPoints(&params[target * m_nParams], m_nPoints, Dim);
    const MatrixType sourcePoints(&params[source * m_nParams], m_nPoints, Dim);
    
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
    
    MatrixType targetCentered(m_nPoints, Dim), sourceCentered(m_nPoints, Dim);
    for (int i = 0; i < m_nPoints; i++) {
        for (int k = 0; k < Dim; k++) {
            targetCentered[i][k] = targetPoints[i][k] - targetCentroid[k];
            sourceCentered[i][k] = sourcePoints[i][k] - sourceCentroid[k];
        }
    }
    
    double rotationParam;
    VectorType translationParam = targetCentroid - sourceCentroid;
    
    MatrixType cov = targetCentered.transpose() * sourceCentered;
    vnl_svd<double> svd(cov);
    MatrixType rotationMatrix = svd.V() * svd.U().transpose();

    rotationParam = acos(rotationMatrix[0][0]);
    
    transformParams[source][0] = translationParam[0];
    transformParams[source][1] = translationParam[1];
    transformParams[source][2] = rotationParam * 180 / 3.141592;
}

void myEnsembleEntropy::ComputePositionalEnsemble(const OptimizerParametersType &params, double &cost, MatrixType &deriv) const {
    
    MatrixType transformParams(m_nSubjects, m_nTransformParams);
    transformParams.fill(0);
    
    for (int i = 1; i < m_nSubjects; i++) {
        EstimateRigidParameters(transformParams, params, 0, i);
        std::cout << "Estimated Transform: " << transformParams << std::endl;
    }
}