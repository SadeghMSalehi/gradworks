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
#include <itkVector.h>
#include <armadillo>
#ifdef __APPLE__
#include "/usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/include/omp.h"
#endif
#include "imageParticleTypes.h"

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
    typedef itk::Vector<float,VDim> PointType;
    typedef std::vector<PointType> PointContainerType;

    void SetImage(ImageType::Pointer image) {
        m_Image = image;
    }
    
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
            double sigmaEstimate = 1.06 * sampleStdev * ::pow(m_NumberOfPoints, -1/5);
            m_SampleSigmas[i] = sigmaEstimate;
        }
    }

    void SetUseAdaptiveSampling(bool adaptiveSampling) {
    	m_AdaptiveSampling = adaptiveSampling;
    }

    ParametersType GetSampleSigmas() {
    	return m_SampleSigmas;
    }


    ParametersType GetInitialParameters() {
        ParametersType initial;
        initial.SetSize(m_Points.size() * 2);
        for (unsigned int i = 0; i < m_Points.size(); i++) {
            for (unsigned int j = 0; j < VDim; j++) {
                initial[VDim*i + j] = m_Points[i][j];
            }
        }
        return initial;
    }

    inline double G(double dxi, double si) const {
        if (dxi > 3*si || dxi < -3*si) {
            return 0;
        }
        return exp(-(dxi*dxi) / (2*si*si)) / (sqrt(2*M_PI)*si);
    }

    inline double weightedG(double dxi, double si, double x, double y, double z) const {
        if (m_Image.IsNull()) {
            cout << "Image is null" << endl;
            exit(0);
        }
        ImageType::IndexType i;
        i[0] = ::round(x);
        i[1] = ::round(y);
        ImageType::PixelType v = 1;
        if (m_Image->GetBufferedRegion().IsInside(i)) {
            v = m_Image->GetPixel(i);
        } else {
            cout << "(" << x << "," << y << ") => " << i << endl;
            cout << m_Image->GetBufferedRegion().GetSize() << endl;
        }

        if (v != v) {
            v = 1;
        }
        v = (v < 1) ? 1 : v;

        if (dxi > 3*si || dxi < -3*si) {
            return 0;
        }

        if (v > 1) {
           // cout << v << endl;
        }

        dxi = dxi * (v / 150);


        double g = exp(-(dxi*dxi) / (2*si*si)) / (sqrt(2*M_PI)*si);
        return g;
    }

    virtual unsigned int GetNumberOfParameters() const {
        return m_Points.size() * VDim;
    }

    /** This method returns the value of the cost function corresponding
     * to the specified parameters.    */
    virtual MeasureType GetValue(const ParametersType & parameters) const {
        MeasureType value;
        DerivativeType derivative;
        derivative.SetSize(GetNumberOfParameters());
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
        int iPtId = -1;

        for (int i = 0; i < nParams; i += VDim) {
            iPtId ++;
            double si = 3; //m_SampleSigmas[iPtId];
            double gSum = 0;
            arma::vec weights;
            weights.zeros(m_NumberOfPoints);
            for (int j = 0, jPtId = 0; j < nParams; j += VDim, jPtId++) {
            	if (j == i) {
            		continue;
            	}
                double dist_ij = 0;
                for (unsigned int k = 0; k < VDim; k++) {
                	dist_ij += (p[i+k] - p[j+k]) * (p[i+k] - p[j+k]);
                }
                dist_ij = ::sqrt(dist_ij);
                double g = 0;
                if (m_AdaptiveSampling) {
                	g = weightedG(dist_ij, si, p[j], p[j+1], 0);
                } else {
                	g = G(dist_ij, si);
                }
                weights[jPtId] = g;
            }
            gSum = arma::sum(weights);
            cost += gSum;

            if (gSum == 0) {
                weights.fill(0);
            } else {
                weights /= gSum;
            }
            

            arma::vec deriv_i;
            deriv_i.zeros(VDim);
            int jPtId = -1;
            for (int j = 0; j < nParams; j += VDim) {
                jPtId ++;
            	if (j == i) {
            		continue;
            	}
            	for (unsigned int k = 0; k < VDim; k++) {
            		deriv_i[k] += (p[i+k] - p[j+k]) * weights[jPtId];
            	}
            }
            deriv_i = - deriv_i / (si*si);
            for (unsigned int k = 0; k < VDim; k++) {
            	derivative[i + k] = deriv_i[k];
            }
        }
        value = cost;
        //cout << "Cost: " << value << "; #: " << p << endl;
    }

protected:
    SurfaceEntropyCostFunction() {
        m_NumberOfPoints = 0;
        m_AdaptiveSampling = true;
    }

    virtual ~SurfaceEntropyCostFunction() {
    }

private:
    PointContainerType m_Points;
    ParametersType m_SampleSigmas;
    int m_NumberOfPoints;
    ImageType::Pointer m_Image;
    bool m_AdaptiveSampling;
};

#endif /* defined(__imageParticles__SurfaceEntropyCostFunction__) */
