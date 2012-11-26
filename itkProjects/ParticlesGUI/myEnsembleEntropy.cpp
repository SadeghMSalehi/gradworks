//
//  myEnsembleEntropy.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/24/12.
//
//

#include "myEnsembleEntropy.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_matrix_inverse.h"


using namespace std;

const int Dim = 2;

myEnsembleEntropy::myEnsembleEntropy() {
    
}

myEnsembleEntropy::~myEnsembleEntropy() {
    
}

unsigned int myEnsembleEntropy::GetNumberOfParameters() const {
    return m_nSubjects * m_nParams;
}

/** This method returns the value of the cost function corresponding
 * to the specified parameters.    */
double myEnsembleEntropy::GetValue(const ParametersType & parameters) const {
    MeasureType value;
    DerivativeType derivative;
    derivative.SetSize(GetNumberOfParameters());
    GetValueAndDerivative(parameters, value, derivative);
    return value;
}

/** This method returns the derivative of the cost function corresponding
 * to the specified parameters.   */
void myEnsembleEntropy::GetDerivative(const ParametersType & parameters,
                                            DerivativeType & derivative) const {
    MeasureType value;
    GetValueAndDerivative(parameters, value, derivative);
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
    
    VNLMatrix targetPoints(params.data_block() + target * m_nParams, m_nPoints, Dim);
    VNLMatrix sourcePoints(params.data_block() + source * m_nParams, m_nPoints, Dim);
    
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
    
    cout << sourcePoints << endl;
    cout << targetPoints << endl;
    
    double rotationParam;
    VNLVector translationParam = targetCentroid - sourceCentroid;
    VNLMatrix cov = targetPoints.transpose() * sourcePoints;
    
    // debug: nan appears with unknown causes
    if (cov.has_nans()) {
        std::cout << "Params: " << params << endl;
        std::cout << "Target: " << targetPoints << std::endl;
        std::cout << "Source: " << sourcePoints << std::endl;
    }

    std::cout << "COV: " << cov << std::endl;
    vnl_svd<double> svd(cov);
    VNLMatrix rotationMatrix = svd.V() * svd.U().transpose();
//    std::cout << rotationMatrix << std::endl;

    // debug: atan2 is more stable than acos
    rotationParam = atan2(rotationMatrix[1][0], rotationMatrix[0][0]);
//    std::cout << rotationParam << std::endl;
//    
//    std::cout << "acos(1) = " << acos(1) << std::endl;
    
    transformParams[source][0] = translationParam[0];
    transformParams[source][1] = translationParam[1];
    transformParams[source][2] = rotationParam;
}

static void ComputeEnsembleEntropy(VNLMatrix& data, const VNLMatrixArray& jacobianList, double &cost, VNLVector& deriv, int nSubj, int nParams) {
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
    
    VNLMatrix covInverse = vnl_matrix_inverse<double>(cov);
    VNLMatrix grad = covInverse * data;
    
    if (grad.has_nans()) {
        cout << "Gradient has Nans: " << grad << endl;
    }
    // multiply jacobian of the function
    const int nDim = jacobianList[1].cols();
    const int nPoints = data.cols() / nDim;
    for (int n = 1; n < data.rows(); n++) {
        int nOffset = n * nParams;
        VNLVector v(nDim);
        for (int i = 0; i < nPoints; i++) {
            const int iOffset = i * nDim;
            for (int k = 0; k < nDim; k++) {
                v[k] = data[n][iOffset + k];
            }
            VNLVector w = jacobianList[n] * v;
            for (int k = 0; k < nDim; k++) {
                deriv[nOffset+iOffset+k] = w[k];
            }
        }
    }
}

void myEnsembleEntropy::GetValueAndDerivative(const OptimizerParametersType &params, double &cost, DerivativeType &deriv) const {
    deriv.set_size(m_nParams * m_nSubjects);
    deriv.fill(0);
    
    if (vnl_has_nan(params)) {
        cout << "NaN detected: " << params << endl;
        return;
    }
    
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
    
//    std::cout << "Transformation: " << transformParams << std::endl;
//    std::cout << "Params: " << params << std::endl;

    VNLMatrix targetPoints(params.data_block(), m_nSubjects, m_nParams);
    vnl_center(targetPoints);
 
    cout << "Points: " << targetPoints << endl;
    ComputeEnsembleEntropy(targetPoints, jacobianList, cost, deriv, m_nSubjects, m_nParams);
}