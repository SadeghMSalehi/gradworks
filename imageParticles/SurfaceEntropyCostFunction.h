//
//  SurfaceEntropyCostFunction.h
//  imageParticles
//
//  Created by Joohwi Lee on 10/22/12.
//
//

#ifndef __imageParticles__SurfaceEntropyCostFunction__
#define __imageParticles__SurfaceEntropyCostFunction__

#include <iostream>
#include <itkSingleValuedCostFunction.h>
#include <armadillo>

template <unsigned int VDim>
class SurfaceEntropyCostFunction: public itk::SingleValuedCostFunction {
public:
    typedef SurfaceEntropyCostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkTypeMacro(Self, Superclass);
    itkNewMacro(Self);

    typedef double MeasureType;
    typedef Superclass::ParametersType ParametersType;
    typedef Superclass::ParametersValueType ParametersValueType;
    typedef itk::Array<ParametersValueType> DerivativeType;
    typedef itk::Point<float,VDim> PointType;
    typedef std::vector<PointType> PointContainerType;

    void Initialize(PointContainerType points) {
        m_NumberOfPoints = points.size();
        m_Points = points;
        m_SampleSigmas.SetSize(points.size());
        for (int i = 0; i < m_NumberOfPoints; i++) {
            arma::running_stat<double> stats;
            PointType a = m_Points[i];
            for (int j = 0; j < m_NumberOfPoints; j++) {
                PointType b = m_Points[j];
                PointType ba = b - a;
                stats(ba.GetNorm());
            }
            double sampleStdev = ::sqrt(stats.var());
            double sigmaEstimate = 1.06 * sampleStdev * ::pow(m_NumberOfPoints, -5);
            m_SampleSigmas[i] = sigmaEstimate;
        }
    }


    ParametersType GetInitialParameters() {
        ParametersType initial;
        initial.SetSize(m_Points.size() * 2);
        for (int i = 0; i < m_Points.size(); i++) {
            for (int j = 0; j < VDim; j++) {
                initial[VDim*i + j] = m_Points[i][j];
            }
        }
        return initial;
    }

    double G(double dxi, double si) const {
        if (dxi > 3*si || dxi < -3*si) {
            return 0;
        }
        return exp(-(dxi*dxi) / (2*si*si)) / (sqrt(2*M_PI)*si);
    }

    virtual unsigned int GetNumberOfParameters() const {
        return m_Points.size() * VDim;
    }

    /** This method returns the value of the cost function corresponding
     * to the specified parameters.    */
    virtual MeasureType GetValue(const ParametersType & parameters) const {
        MeasureType value;
        DerivativeType derivative;
        GetValueAndDerivative(parameters, value, derivative);
        return value;
    }

    /** This method returns the derivative of the cost function corresponding
     * to the specified parameters.   */
    virtual void GetDerivative(const ParametersType & parameters,
                               DerivativeType & derivative) const {
        MeasureType value;
        GetValueAndDerivative(parameters, value, derivative);
    }

    /** This method returns the value and derivative of the cost function corresponding
     * to the specified parameters    */
    virtual void GetValueAndDerivative(const ParametersType & p,
                                       MeasureType & value,
                                       DerivativeType & derivative) const {
        MeasureType cost = 0;
        int nParams = p.GetSize();
        for (int i = 0, ptId = 0; i < nParams; i += VDim, ptId++) {
            double si = m_SampleSigmas[ptId];
            double gSum = 0;
            std::vector<double> distVector;
            std::vector<double> gVector;
            distVector.resize(m_NumberOfPoints);
            for (int j = 0; j < nParams; j += VDim) {
                double dist = 0;
                for (int k = 0; k < VDim; k++) {
                    dist += (p[i+k] - p[j+k]) * (p[i+k] - p[j+k]);
                }
                dist = ::sqrt(dist);
                double g = G(dist, si);
                distVector[j] = dist;
                gVector[j] = g;
                gSum += g;
            }
            cost += ::log(gSum);
            for (int k = 0; k < m_NumberOfPoints; k++) {
                distVector[k] * gVector[j] / gSum;
            }
        }
    }

protected:
    SurfaceEntropyCostFunction() {
        m_NumberOfPoints = 0;
    }
    virtual ~SurfaceEntropyCostFunction() {
    }

private:
    PointContainerType m_Points;
    ParametersType m_SampleSigmas;
    int m_NumberOfPoints;
};

#endif /* defined(__imageParticles__SurfaceEntropyCostFunction__) */