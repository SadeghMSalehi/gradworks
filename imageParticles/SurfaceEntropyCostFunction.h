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
#include "itkRescaleIntensityImageFilter.h"

#define __CLAMP(x,s)((x<-3*si?0:(x>3*si?0:x)))

template<unsigned int VDim>
class SurfaceEntropyCostFunction: public itk::SingleValuedCostFunction {
public:
	typedef SurfaceEntropyCostFunction Self;
	typedef itk::SingleValuedCostFunction Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkTypeMacro(Self, Superclass)
	;itkNewMacro(Self)
	;

	typedef double MeasureType;
	typedef Superclass::ParametersType ParametersType;
	typedef Superclass::ParametersValueType ParametersValueType;
	typedef itk::Array<ParametersValueType> DerivativeType;
	typedef itk::Vector<float, VDim> PointType;
	typedef std::vector<PointType> PointContainerType;

    itkSetMacro(CutoffDistance, double);
    itkSetMacro(Sigma, double);
    itkSetMacro(PhantomCutoffDistance, double);
    itkSetMacro(MaxKappa, double);

	void SetDistanceVectorImage(DistanceVectorImageType::Pointer distVector) {
		m_distVectorMap = distVector;
	}

	void SetDistanceVectorMagnitudeImage(ImageType::Pointer distVectorMag) {
		m_distVectorMagnitudeMap = distVectorMag;
		m_Interpolator = InterpolatorType::New();
		m_Interpolator->SetInputImage(m_distVectorMagnitudeMap);
	}

	void SetUsePhantomParticles(bool usePhantoms) {
		m_UsePhantomParticles = usePhantoms;
	}

	void SetPhantomParticles(ListOfPointVectorType* phantomParticles) {
		m_PhantomParticles = phantomParticles;
	}

	void SetImage(ImageType::Pointer image) {
		typedef itk::RescaleIntensityImageFilter<ImageType> RescaleFilter;
		RescaleFilter::Pointer filter = RescaleFilter::New();
		filter->SetInput(image);
		filter->SetOutputMinimum(1);
		filter->SetOutputMaximum(m_MaxKappa);
		filter->Update();
		m_Image = filter->GetOutput();
	}

	ImageType::Pointer GetImage() {
		return m_Image;
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
			double sigmaEstimate = 1.06 * sampleStdev
					* ::pow(m_NumberOfPoints, -1 / 5);
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
				initial[VDim * i + j] = m_Points[i][j];
			}
		}
		cout << "Initial Parameters: " << initial << endl;
		return initial;
	}

	inline double G(double dxi, double si) const {
		if (dxi > 3 * si || dxi < -3 * si) {
			return 0;
		}
		return exp(-(dxi * dxi) / (2 * si * si)) / (sqrt(2 * M_PI) * si);
	}

	inline double weightedG(double dxi, double si, double x, double y,
			double z) const {
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
//            cout << "(" << x << "," << y << ") => " << i << endl;
//            cout << m_Image->GetBufferedRegion().GetSize() << endl;
		}

		if (dxi > 3 * si || dxi < -3 * si) {
			return 0;
		}

		dxi = dxi * (2 * (v - 1) + 1);

		double g = exp(-(dxi * dxi) / (2 * si * si)) / (sqrt(2 * M_PI) * si);
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

	void GetUnguardedValueAndDerivative(const ParametersType & p,
			MeasureType & value, DerivativeType & derivative) const {
		MeasureType cost = 0, rCost = 0;
		int nParams = p.GetSize();
		int iPtId = -1;

		derivative.SetSize(GetNumberOfParameters());
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		interpolator->SetInputImage(m_distVectorMagnitudeMap);

		//        cout << "Parameters: " << p << endl;
		//        cout << "Derivatives: " << derivative << "[" << derivative.GetSize() << "]" << endl;
		//        cout << "Number of Points: " << m_NumberOfPoints << endl;

		for (int i = 0; i < nParams; i += VDim) {
			iPtId++;
			double si = 15; //m_SampleSigmas[iPtId];
			double gSum = 0;
			arma::vec weights;

			weights.zeros(m_NumberOfPoints);

			ImageType::IndexType ii;
			ii[0] = ::round(p[i]);
			ii[1] = ::round(p[i + 1]);
			ImageType::PixelType v = 1;
			if (m_Image->GetBufferedRegion().IsInside(ii)) {
				v = m_Image->GetPixel(ii);
			}

			for (int j = 0, jPtId = 0; j < nParams; j += VDim, jPtId++) {
				if (j == i) {
					continue;
				}
				ImageType::IndexType jj;
				jj[0] = ::round(p[j]);
				jj[1] = ::round(p[j + 1]);
				ImageType::PixelType w = 1;
				if (m_Image->GetBufferedRegion().IsInside(ii)) {
					w = m_Image->GetPixel(ii);
				}

				double dist_ij = 0;
				for (unsigned int k = 0; k < VDim; k++) {
					dist_ij += (p[i + k] - p[j + k]) * (p[i + k] - p[j + k]);
				}
				dist_ij = ::sqrt(dist_ij);
				double g = 0;
				if (m_AdaptiveSampling) {
					//g = weightedG(dist_ij, si, p[i], p[i+1], 0);
					g = G(w * dist_ij, si);
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

			//            cout << "Computing derivatives ..." << endl;

			if (derivative.GetSize() > 0) {
				arma::vec deriv_i;
				deriv_i.zeros(VDim);
				int jPtId = -1;
				for (int j = 0; j < nParams; j += VDim) {
					jPtId++;
					if (j == i) {
						continue;
					}
					for (unsigned int k = 0; k < VDim; k++) {
						deriv_i[k] += (p[i + k] - p[j + k]) * weights[jPtId];
					}
				}
				deriv_i = -deriv_i / (si * si);

				ImageType::IndexType ii;
				ii[0] = p[i];
				ii[1] = p[i + 1];
				ImageType::OffsetType distOffset = m_distVectorMap->GetPixel(
						ii);
				for (unsigned int k = 0; k < VDim; k++) {
					derivative[i + k] = deriv_i[k]; // + abs(ri) / 100 * distOffset[k];
				}
			}
		}
		value = cost; // + rCost;
		//        cout << "Cost: " << value << "; #: " << p << endl;

	}

	inline double computeEntropy(const ParametersType& p, int i, int j,
			double kappa = 1) const {
		double dist2 = 0;
		double si = m_Sigma / kappa;
		for (int k = 0; k < 2; k++) {
			dist2 += (p[2 * i + k] - p[2 * j + k])
					* (p[2 * i + k] - p[2 * j + k]);
		}
		if (dist2 > (m_CutoffDistance * m_CutoffDistance)) {
			return 0;
		} else {
			double force = exp(-dist2 / (si * si));
			if (force != force) {
				cout << "Force is NaN; " << i << ", " << j << "; " << dist2
						<< ";" << p << endl;
				exit(0);
			}
			return force;
		}
	}

	inline double computePhantomParticleEntropy(const ParametersType& p, int i,
			int j, double kappa = 1) const {
		double dist2 = 0;
		double si = m_Sigma / kappa;

		arma::vec phantom;
		phantom << m_PhantomParticles->at(1).at(2 * j)
				<< m_PhantomParticles->at(1).at(2 * j + 1);

		for (int k = 0; k < 2; k++) {
			dist2 += (p[2 * i + k] - phantom[k]) * (p[2 * i + k] - phantom[k]);
		}

		if (dist2 > (m_PhantomCutoffDistance * m_PhantomCutoffDistance)) {
			return 0;
		} else {
			double force = exp(-dist2 / (si * si));
			if (force != force) {
				cout << "Force is NaN; " << i << ", " << j << "; " << dist2
						<< ";" << p << endl;
				exit(0);
			}
			return force;
		}
	}

	inline double computePhantomEntropy(const ParametersType& p, int i) const {
		ContinuousIndexType idx;
		idx[0] = p[2 * i];
		idx[1] = p[2 * i + 1];
		if (!m_Interpolator->IsInsideBuffer(idx)) {
			return exp(5 * 5 / 9);
		}
		ImageType::PixelType dist = m_Interpolator->EvaluateAtContinuousIndex(
				idx);
		if (dist != dist) {
			cout << "dist is NaN" << endl;
			return exp(5 * 5);
		}
		if (dist > 0) {
			if (dist > 5) {
				return exp(5 * 5);
			} else {
				return exp(dist * dist / (m_Sigma * m_Sigma));
			}
		} else {
			return exp(-dist * dist / (m_Sigma * m_Sigma));
		}
	}

	arma::vec normalizeGradient(ImageType::OffsetType gradient) const {
		arma::vec norm;
		norm << gradient[0] << gradient[1];
		if (arma::norm(norm, 2) == 0) {
			norm.zeros(2);
			return norm;
		}
		return norm / arma::norm(norm, 2);
	}

	void GetPhantomValueAndDerivative(const ParametersType& p,
			MeasureType & value, DerivativeType & derivative) const {
		int nParams = p.GetSize();
		int nPoints = nParams / 2;
        if (m_PhantomParticles->size() == 0) {
            cout << "Phantom Particles is not initialized!!" << endl;
            exit(0);
        }
		int nPhantoms = m_PhantomParticles->at(0).size() / 2;

		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		interpolator->SetInputImage(m_Image);

		arma::mat entropies;
		entropies.zeros(nPoints, nPoints + nPhantoms);

		double cost = 0;
		for (int i = 0; i < nPoints; i++) {
			double sum = 0;
			ImageType::ValueType kappa = 1;
			if (m_AdaptiveSampling) {
				ContinuousIndexType idx;
				idx[0] = p[2 * i];
				idx[1] = p[2 * i + 1];
				kappa = interpolator->EvaluateAtContinuousIndex(idx);
			}
			for (int j = 0; j < nPoints; j++) {
				if (i == j) {
					continue;
				}
				entropies.at(i, j) = computeEntropy(p, i, j, kappa);
				sum += entropies.at(i, j);
			}
			for (int j = 0; j < nPhantoms; j++) {
				entropies.at(i, j + nPoints) = computePhantomParticleEntropy(p,
						i, j, kappa);
				sum += entropies.at(i, j + nPoints);
			}
			if (sum > 0) {
				entropies.row(i) = entropies.row(i) / sum;
			}
			cost += sum;
		}

		for (int i = 0; i < nPoints; i++) {
			for (int j = 0; j < nPoints; j++) {
				if (i == j) {
					continue;
				}
				for (int k = 0; k < 2; k++) {
					derivative[i * 2 + k] -= entropies.at(i, j)
							* (p[2 * i + k] - p[2 * j + k]);
				}
			}
			for (int j = 0; j < nPhantoms; j++) {
				for (int k = 0; k < 2; k++) {
					derivative[i * 2 + k] -= entropies.at(i, j + nPoints)
							* (p[2 * i + k] - m_PhantomParticles->at(0).at(2 * j + k));
				}
			}
		}
	}

	void GetGuardedValueAndDerivative(const ParametersType & p,
			MeasureType & value, DerivativeType & derivative) const {
		int nParams = p.GetSize();
		int nPoints = nParams / 2;

		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		interpolator->SetInputImage(m_Image);

		arma::mat entropies;
		arma::vec phantomEntropies, entropiesSum;
		entropies.zeros(nPoints, nPoints + 1);
		phantomEntropies.zeros(nPoints);
		entropiesSum.zeros(nPoints);

		double cost = 0;
		for (int i = 0; i < nPoints; i++) {
			double sum = 0;

			for (int j = 0; j < nPoints; j++) {
				if (i == j) {
					continue;
				}
				ImageType::ValueType kappa = 1;
				if (m_AdaptiveSampling) {
					ContinuousIndexType idx;
					idx[0] = p[2 * j];
					idx[1] = p[2 * j + 1];
					kappa = interpolator->EvaluateAtContinuousIndex(idx);
				}
				entropies.at(i, j) = computeEntropy(p, i, j, kappa);
				sum += entropies.at(i, j);
			}
			entropies.at(i, nPoints) = computePhantomEntropy(p, i);
			sum += entropies.at(i, nPoints);
			if (sum > 0) {
				entropies.row(i) = entropies.row(i) / sum;
			}
			if (sum != sum) {

				cout << "sum is NaN; " << "phantom = "
						<< entropies.at(i, nPoints) << "; " << p[2 * i] << ","
						<< p[2 * i + 1] << endl;
			}
			cost += sum;
		}
		for (int i = 0; i < nPoints; i++) {
			for (int j = 0; j < nPoints; j++) {
				if (i == j) {
					continue;
				}
				for (int k = 0; k < 2; k++) {
					derivative[i * 2 + k] -= entropies.at(i, j)
							* (p[2 * i + k] - p[2 * j + k]);
					if (derivative[i * 2 + k] != derivative[i * 2 + k]) {
						cout << "Derivative is NaN; " << i << "; "
								<< entropies.at(i, j) << endl;
					}
				}
			}
			ContinuousIndexType idx;
			idx[0] = p[2 * i];
			idx[1] = p[2 * i + 1];
			if (m_Interpolator->IsInsideBuffer(idx)) {
				double dist = m_Interpolator->EvaluateAtContinuousIndex(idx);
				if (dist > 0) {
					ImageType::IndexType gIdx;
					gIdx[0] = ::round(p[2 * i]);
					gIdx[1] = ::round(p[2 * i + 1]);
					ImageType::OffsetType gradientPixel =
							m_distVectorMap->GetPixel(gIdx);
					arma::vec normalizedGradient = normalizeGradient(
							gradientPixel);
					for (int k = 0; k < 2; k++) {
						derivative[i * 2 + k] -= (entropies.at(i, nPoints)
								* normalizedGradient[k]);
						if (derivative[i * 2 + k] != derivative[i * 2 + k]) {
							cout << "Derivative is NaN; " << i
									<< "; phantomForce = "
									<< entropies.at(i, nPoints)
									<< "; Gradient = " << normalizedGradient
									<< "; " << endl;
						}
					}
				} else if (dist < 0) {
					ImageType::IndexType gIdx;
					gIdx[0] = ::round(p[2 * i]);
					gIdx[1] = ::round(p[2 * i + 1]);
					ImageType::OffsetType gradientPixel =
							m_distVectorMap->GetPixel(gIdx);
					arma::vec normalizedGradient = normalizeGradient(
							gradientPixel);
					for (int k = 0; k < 2; k++) {
						derivative[i * 2 + k] -= (entropies.at(i, nPoints)
								* normalizedGradient[k]);
						if (derivative[i * 2 + k] != derivative[i * 2 + k]) {
							cout << "Derivative is NaN; " << i
									<< "; phantomForce = "
									<< entropies.at(i, nPoints)
									<< "; Gradient = " << normalizedGradient
									<< "; " << endl;
						}
					}
				}
			}
		}
		value = cost;
		if (value != value) {
			cout << "Value is NAN; " << p << endl;
			exit(0);
		}
		for (int i = 0; i < nParams; i++) {
			if (derivative[i] != derivative[i]) {
				cout << i << " th derivative is nan; " << p << ";" << derivative
						<< endl;
			}
		}
	}

	/** This method returns the value and derivative of the cost function corresponding
	 * to the specified parameters    */
	virtual void GetValueAndDerivative(const ParametersType & p,
			MeasureType & value, DerivativeType & derivative) const {
		derivative.SetSize(GetNumberOfParameters());
		derivative.Fill(0);
		if (m_PhantomParticles != NULL) {
            cout << "Using phantom particles ..." << endl;
            GetPhantomValueAndDerivative(p, value, derivative);
        } else {
            cout << "Using boundary guarded particles ..." << endl;            
            GetGuardedValueAndDerivative(p, value, derivative);
        }
	}

    void Print() {
        cout << "Sigma: " << m_Sigma << endl;
        cout << "Max Kappa: " << m_MaxKappa << endl;
        cout << "Cutoff Distance: " << m_CutoffDistance << endl;
        cout << "Cutoff Distance for boundary: " << m_PhantomCutoffDistance << endl;
        cout << "Adaptive Sampling: " << m_AdaptiveSampling << endl;
    }

protected:
	SurfaceEntropyCostFunction() {
		m_NumberOfPoints = 0;
		m_AdaptiveSampling = true;
        m_PhantomParticles = NULL;

        m_Sigma = 7.0;
        m_MaxKappa = sqrt(2);
        m_CutoffDistance = 15;
        m_PhantomCutoffDistance = 3;
	}

	virtual ~SurfaceEntropyCostFunction() {
	}

private:
	PointContainerType m_Points;
	ParametersType m_SampleSigmas;
	int m_NumberOfPoints;
	ImageType::Pointer m_Image;
	ImageType::Pointer m_distVectorMagnitudeMap;
	InterpolatorType::Pointer m_Interpolator;
	DistanceVectorImageType::Pointer m_distVectorMap;
	ListOfPointVectorType* m_PhantomParticles;
	bool m_AdaptiveSampling;
	bool m_UsePhantomParticles;

    double m_Sigma;
    double m_MaxKappa;
    double m_CutoffDistance;
    double m_PhantomCutoffDistance;
};

#endif /* defined(__imageParticles__SurfaceEntropyCostFunction__) */
