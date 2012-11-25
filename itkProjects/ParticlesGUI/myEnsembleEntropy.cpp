//
//  myEnsembleEntropy.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/24/12.
//
//

#include "myEnsembleEntropy.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

const int Dim = 2;

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

void myEnsembleEntropy::SetVariableCounts(int s, int p, int np) {
    m_nSubjects = s;
    m_nPoints = p;
    m_nParams = np;
}

void myEnsembleEntropy::EstimateRigidParameters(VNLMatrix &transformParams, const OptimizerParametersType &params, int target, int source) const {
    const int Dim = 2;
    
    VNLMatrix targetPoints(&params[target * m_nParams], m_nPoints, Dim);
    VNLMatrix sourcePoints(&params[source * m_nParams], m_nPoints, Dim);
    
    // make those centroid match
    VNLVector targetCentroid(2), sourceCentroid(2);
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
    VNLVector translationParam = targetCentroid - sourceCentroid;
    
    VNLMatrix cov = targetPoints.transpose() * sourcePoints;
    vnl_svd<double> svd(cov);
    VNLMatrix rotationMatrix = svd.V() * svd.U().transpose();

    rotationParam = acos(rotationMatrix[0][0]);
    
    transformParams[source][0] = translationParam[0];
    transformParams[source][1] = translationParam[1];
    transformParams[source][2] = rotationParam;
}

static void ComputeEnsembleEntropy(VNLMatrix& data, const VNLMatrixArray& jacobianList, double &cost, VNLMatrix& deriv) {
    VNLMatrix cov = data * data.transpose();
    
    // cost function is the sum of log of eigenvalues
    vnl_symmetric_eigensystem<double> eigen(cov);
    cost = 0;
    for (int i = 0; i < eigen.D.size(); i++) {
        if (eigen.D[i] > 0) {
            cost += log(eigen.D[i]);
        }
    }

    // relaxation parameter
    double alpha = 1;
    // gradient computation
    for (int i = 0; i < cov.rows(); i++) {
        cov[i][i] += alpha;
    }
    VNLMatrix grad = cov * data;
    
    // multiply jacobian of the function
    const int nDim = jacobianList[1].cols();
    const int nPoints = data.cols() / nDim;
    for (int n = 1; n < data.rows(); n++) {
        VNLVector v(nDim);
        for (int i = 0; i < nPoints; i++) {
            const int iOffset = i * nDim;
            for (int k = 0; k < nDim; k++) {
                v[k] = data[n][iOffset + k];
            }
            VNLVector w = jacobianList[n] * v;
            for (int k = 0; k < nDim; k++) {
                deriv[n][iOffset+k] = w[k];
            }
        }
    }
}

void myEnsembleEntropy::ComputePositionalEnsemble(const OptimizerParametersType &params, double &cost, VNLMatrix &deriv) const {
    
    VNLMatrix transformParams(m_nSubjects, m_nTransformParams);
    transformParams.fill(0);

    const int nDim = 2;
    VNLMatrixArray jacobianList;
    jacobianList.push_back(VNLMatrix());
    for (int i = 1; i < m_nSubjects; i++) {
        EstimateRigidParameters(transformParams, params, 0, i);
        VNLMatrix jacobian(nDim, nDim);
        jacobian[0][0] = cos(transformParams[i][2]);
        jacobian[0][1] = -sin(transformParams[i][2]);
        jacobian[1][0] = sin(transformParams[i][2]);
        jacobian[1][1] = cos(transformParams[i][2]);
        jacobianList.push_back(jacobian);
    }

    VNLMatrix targetPoints(&params[0], m_nSubjects, m_nParams);
    vnl_center(targetPoints);
    
    ComputeEnsembleEntropy(targetPoints, jacobianList, cost, deriv);
}