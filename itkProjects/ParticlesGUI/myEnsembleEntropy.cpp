//
//  myEnsembleEntropy.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/24/12.
//
//

#include "myEnsembleEntropy.h"
#include "itkGradientImageFilter.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_matrix_inverse.h"


using namespace std;

const int Dim = 2;

myEnsembleEntropy::myEnsembleEntropy() {
    m_nImageParamDims = 0;
    m_gradientScale = 0;
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

    typedef itk::GradientImageFilter<SliceType,double,double> GradientFilterType;
    for (int n = 0; n < m_ImageList->size(); n++) {
        SliceType::Pointer image = m_ImageList->at(n)->GetSlice();

        // compute gradient image
        GradientFilterType::Pointer filter = GradientFilterType::New();
        filter->SetInput(image);
        filter->Update();
        GradientImageType::Pointer gradImg = filter->GetOutput();
        m_GradientImageList.push_back(gradImg);

        // add slice interpolator
        SliceInterpolatorType::Pointer interp = SliceInterpolatorType::New();
        interp->SetInputImage(image);
        m_ImageInterpolatorList.push_back(interp);
    }
}

/**
 * must provide exactly matching number of subjects and points with the already provided image list
 */
void myEnsembleEntropy::SetInitialPositions(const OptimizerParametersType &params, int nSubject, int nPoints, int nParams) {
    m_nSubjects = nSubject;
    m_nPoints = nPoints;
    m_nParams = nParams;
    m_nImageParams = m_nPoints * m_nImageParamDims;
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
//  debug: sourcePoints are transformed onto corresponding targetPoints
//    cout << sourcePoints << endl;
//    cout << targetPoints << endl;

    double rotationParam;
    VNLVector translationParam = targetCentroid - sourceCentroid;
    VNLMatrix cov = sourcePoints.transpose() * targetPoints;
    
    // debug: nan appears with unknown causes
    if (cov.has_nans()) {
        std::cout << "Params: " << params << endl;
        std::cout << "Target: " << targetPoints << std::endl;
        std::cout << "Source: " << sourcePoints << std::endl;
    }

    // debug: if the covariance is computed correctly
//    std::cout << "COV: " << cov << std::endl;
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


void myEnsembleEntropy::ComputeEnsembleEntropy(VNLMatrix& data, const VNLMatrixArray& jacobianList, const OptimizerParametersType& params, const std::vector<GradientImageType::Pointer>& gradientList, double &cost, VNLVector& deriv, const int nSubj, const int nPoints, const int nParams, const int nImageParams) const {
    // debug: positional entropy should work without image params
//    cout << "nPaarms: " << nParams << endl;
//    cout << "Data: " << data << endl;
    VNLMatrix cov = data * data.transpose();

    // relaxation parameter for singular matrix
    // this produce pseudo-inverse matrix,
    // if alpha is zero, the inverse matrix will have really high numbers
    double alpha = 1;
    for (int i = 0; i < cov.rows(); i++) {
        cov[i][i] += alpha;
    }

    // cost function is the sum of log of eigenvalues
    vnl_symmetric_eigensystem<double> eigen(cov);

    // debug: eigenvalues for singular matrix; still producing eigenvalue with zero
//    cout << "Eigenvalues: " << eigen.D << endl;
//    cout << "Eigenvectors: " << eigen.V << endl;

    cost = 0;
    for (int i = 0; i < eigen.D.size(); i++) {
        if (eigen.D[i] > 1) {
            cost += log(eigen.D[i]);
        }
    }

    // debug: let's see difference between image gradient and position gradient
//    cout << "Covariance:" << cov << endl;

    VNLMatrix covInverse = vnl_matrix_inverse<double>(cov);
    VNLMatrix grad = covInverse * data;
    // debug: gradient is the direction minimizing covariance matrix
//    cout << "Ensemble COV:" << cov << endl;
//    cout << "Ensemble InverseCOV: " << covInverse << endl;
//    cout << "gradient: " << grad << endl;
    if (grad.has_nans()) {
        cout << "Gradient has Nans: " << grad << endl;
    }
    // multiply jacobian of the function
    const int nDim = jacobianList[1].cols();
    for (int n = 1; n < data.rows(); n++) {
        VNLVector v(nDim);
        for (int i = 0; i < nPoints; i++) {
            // compute positional jacobian gradient
            const int iOffset = i * nDim;
            for (int k = 0; k < nDim; k++) {
                v[k] = grad[n][iOffset + k];
            }
            VNLVector w = jacobianList[n] * v;
            if (m_nImageParamDims > 0) {
                GradientImageType::IndexType idx;
                idx[0] = params[n*nParams + i*nDim];
                idx[1] = params[n*nParams + i*nDim + 1];
                if (gradientList[n]->GetBufferedRegion().IsInside(idx)) {
                    for (int k = 0; k < nDim; k++) {
                        w[k] = (1 - m_gradientRatio) * w[k] + (m_gradientRatio * m_gradientScale) * gradientList[n]->GetPixel(idx)[k] * grad[n][m_nParams + i * m_nImageParamDims];
                    }
                }
            }
            for (int k = 0; k < nDim; k++) {
                deriv[n*nParams+iOffset+k] = -w[k];
            }
        }
    }

    // debug: if derivative is correct
//    cout << "Gradient1: " << grad << endl;
//    cout << "Gradient2: " << deriv << endl;
}

static inline void rotate(double x0, double y0, double theta, double &x1, double &y1) {
    x1 = cos(theta)*x0 - sin(theta)*y0;
    y1 = sin(theta)*x0 + cos(theta)*y0;
}

static void TransformData(VNLMatrix& transformParams, VNLMatrix& data, int nPoints) {
    // debug: if transformation is correct
//    cout << transformParams << endl;

    for (int i = 1; i < data.rows(); i++) {
//        cout << "Before Transform: " << data.get_row(i) << endl;
        for (int j = 0; j < nPoints; j++) {
            double x, y;
            rotate(data[i][2*j] - transformParams[i][0], data[i][2*j+1] - transformParams[i][1], transformParams[i][2], x, y);
            data[i][2*j] = x;
            data[i][2*j+1] = y;
        }
//        cout << "After Transform: " << data.get_row(i) << endl;
    }
}


// what about orientation & scale for patch?
// just sample a pixel;
void myEnsembleEntropy::SampleImage(int subj, ContinuousIndexType& pos, double* out) const {
    out[0] = m_ImageInterpolatorList[subj]->EvaluateAtContinuousIndex(pos);
}

// sample intensity data for every point (original point, not transformed one)
void myEnsembleEntropy::SampleImageData(const OptimizerParametersType& params, VNLMatrix& data) const {
    for (int i = 0; i < m_nSubjects; i++) {
        for (int j = 0; j < m_nPoints * m_nImageParamDims; j++) {
            ContinuousIndexType pos;
            pos[0] = params[m_nSubjects*i + j*2];
            pos[1] = params[m_nSubjects*i + j*2 + 1];
            SampleImage(i, pos, &data[i][j*m_nImageParamDims+m_nParams]);
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
    VNLMatrix pointData(m_nSubjects, m_nParams + m_nImageParams);
    for (int i = 0; i < m_nSubjects; i++) {
        for (int j = 0; j < m_nParams; j++) {
            pointData[i][j] = params[i * m_nParams + j];
        }
    }

    // debug: data must be centered in the shape space
//    cout << "Before Transform:" << pointData << endl;
    TransformData(transformParams, pointData, m_nPoints);
//    cout << "After Transform: " << pointData << endl;

    if (m_nImageParamDims > 0) {
        SampleImageData(params, pointData);
    }


    // debug: data must be centered in the shape space
//    cout << "Before Center:" << pointData << endl;
    vnl_center(pointData);
//    cout << "After Center: " << pointData << endl;
    // debug: avoid NaN targetPoints;
    // targetPoints is transformed with estimation
//    cout << "Points: " << targetPoints << endl;
    ComputeEnsembleEntropy(pointData, jacobianList, params, m_GradientImageList, cost, deriv, m_nSubjects, m_nPoints, m_nParams, m_nImageParams);
    // debug: to see if numerical computation can be stabilized if minor errors are ignored => numerical stability is more related to the algorithm itself
//    for (int i = 0; i < deriv.GetSize(); i++) {
//        if (deriv[i] < 1e-16 && deriv[i] > -1e-16) {
//            deriv[i] = 0;
//        }
//    }
}