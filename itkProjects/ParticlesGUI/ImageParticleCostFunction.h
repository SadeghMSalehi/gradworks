//
//  ImageParticleCostFunction.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/24/12.
//
//

#ifndef ParticlesGUI_ImageParticleCostFunction_h
#define ParticlesGUI_ImageParticleCostFunction_h

class ImageEntropyCostFunction: public itk::SingleValuedCostFunction {
public:
	typedef ImageEntropyCostFunction Self;
	typedef itk::SingleValuedCostFunction Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
    typedef vnl_vector<itk::OffsetValueType> OffsetVectorType;
    typedef vnl_vector<double> VectorType;
	typedef vnl_matrix<double> MatrixType;

    typedef itk::CovariantVector<double,2> GradientType;
    typedef itk::Image<GradientType,2> GradientImageType;
    typedef itk::GradientRecursiveGaussianImageFilter<SliceType,GradientImageType> GradientImageFilter;
    typedef itk::VectorMagnitudeImageFilter<GradientImageType,SliceType> VectorMagnitudeImageFilter;

	typedef itk::SignedDanielssonDistanceMapImageFilter<LabelSliceType, SliceType> DistanceMapFilter;
	typedef DistanceMapFilter::VectorImageType DistanceVectorImageType;
	typedef itk::LinearInterpolateImageFunction<SliceType,float> InterpolatorType;
	typedef InterpolatorType::ContinuousIndexType ContinuousIndexType;
	typedef std::vector<SliceType::Pointer> ImageList;
	typedef std::vector<DistanceVectorImageType::Pointer> DistanceVectorList;
	typedef std::vector<InterpolatorType::Pointer> InterpolatorList;

	itkTypeMacro(Self, Superclass);
	itkNewMacro(Self);

	typedef double MeasureType;
	typedef Superclass::ParametersType ParametersType;
	typedef Superclass::ParametersValueType ParametersValueType;
	typedef itk::Array<ParametersValueType> DerivativeType;

    itkSetMacro(CutoffDistance, double);
    itkSetMacro(Sigma, double);
    itkSetMacro(PhantomCutoffDistance, double);
    itkSetMacro(MaxKappa, double);
    itkSetMacro(EnsembleFactor, double);
    itkSetMacro(GradientSigma, double);


    void Clear() {
        m_KappaMaps.clear();
        m_KappaMapInterpolators.clear();
        m_DistanceMaps.clear();
        m_DistanceMapInterpolators.clear();
        m_DistanceVectorMaps.clear();
        m_Points.clear();
    }

    void ConstructDistanceMap(LabelSliceType::Pointer shapeMask, SliceType::Pointer& shapeDistanceMap, DistanceVectorImageType::Pointer& shapeDistanceVectorMap) {
        DistanceMapFilter::Pointer distmapFilter = DistanceMapFilter::New();
        distmapFilter->SetInput(shapeMask);
        distmapFilter->Update();
        shapeDistanceMap = distmapFilter->GetOutput();
        shapeDistanceVectorMap = distmapFilter->GetVectorDistanceMap();
    }

    bool AddSubjects(ImageContainer::Pointer image) {
        cout << "adding subject begin!" << endl;

        // generate kappa map from gradient magnitude image
        GradientImageFilter::Pointer gradFilter = GradientImageFilter::New();
        gradFilter->SetInput(image->GetSlice());
        gradFilter->SetSigma(m_GradientSigma);
        gradFilter->Update();

        VectorMagnitudeImageFilter::Pointer magFilter = VectorMagnitudeImageFilter::New();
        magFilter->SetInput(gradFilter->GetOutput());
        magFilter->Update();

        // rescale kappa map
		typedef itk::RescaleIntensityImageFilter<SliceType> RescaleFilter;
		RescaleFilter::Pointer filter = RescaleFilter::New();
		filter->SetInput(magFilter->GetOutput());
		filter->SetOutputMinimum(1);
		filter->SetOutputMaximum(m_MaxKappa);
		filter->Update();
        SliceType::Pointer kappaMap = filter->GetOutput();
        image->AddDerivedView(image->GetName() + "/kappaMap", ImageContainer::CreateBitmap(kappaMap));

        // adding attribute map
        InterpolatorType::Pointer kappaInterpolator = InterpolatorType::New();
        kappaInterpolator->SetInputImage(kappaMap);
        m_KappaMapInterpolators.push_back(kappaInterpolator);
        m_KappaMaps.push_back(kappaMap);


        // adding shape distance map
        SliceType::Pointer shapeDistanceMap;
        DistanceVectorImageType::Pointer shapeDistanceVectorMap;
        ConstructDistanceMap(image->GetLabelSlice(), shapeDistanceMap, shapeDistanceVectorMap);
        image->AddDerivedView(image->GetName() + "/distanceMap", ImageContainer::CreateBitmap(shapeDistanceMap));


        // compute distance map bound for an object
        typedef itk::ThresholdImageFilter<SliceType> ThresholdFilterType;
        ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
        threshold->SetInput(shapeDistanceMap);
        threshold->ThresholdBelow(0);
        threshold->SetOutsideValue(0);
        threshold->Update();

        ThresholdFilterType::Pointer threshold2 = ThresholdFilterType::New();
        threshold2->SetInput(threshold->GetOutput());
        threshold2->ThresholdAbove(0);
        threshold2->SetOutsideValue(1);
        threshold2->Update();
        image->AddDerivedView(image->GetName() + "/distanceMapBound", ImageContainer::CreateBitmap(threshold2->GetOutput()));


        // create interpolator for distance map and distance vector map
        InterpolatorType::Pointer shapeDistanceInterpolator = InterpolatorType::New();
        shapeDistanceInterpolator->SetInputImage(shapeDistanceMap);
        m_DistanceMaps.push_back(shapeDistanceMap);
        m_DistanceVectorMaps.push_back(shapeDistanceVectorMap);
        m_DistanceMapInterpolators.push_back(shapeDistanceInterpolator);

        m_nSubjects = m_DistanceMaps.size();

        cout << "adding subject done!" << endl;

        return true;
    }


	void SetUseAdaptiveSampling(bool adaptiveSampling) {
		m_AdaptiveSampling = adaptiveSampling;
	}


	inline double G(double dxi, double si) const {
		if (dxi > 3 * si || dxi < -3 * si) {
			return 0;
		}
		return exp(-(dxi * dxi) / (2 * si * si)) / (sqrt(2 * M_PI) * si);
	}

    void SetNumberOfParameters(int nSubjects, int nPoints, int nTotalVars) {
        m_nVars = nPoints * Dimensions;
        m_nPoints = nPoints;
        m_nSubjects = nSubjects;
        m_nTotalVars = nTotalVars;

    }

	virtual unsigned int GetNumberOfParameters() const {
		return m_nTotalVars;
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

	inline double computeEntropy(const ParametersType& p, int nOffset, int i, int j,
                                 double kappa = 1) const {
		double dist2 = 0;
		double si = m_Sigma / kappa;
		for (int k = 0; k < 2; k++) {
			dist2 += (p[nOffset + 2 * i + k] - p[nOffset + 2 * j + k])
            * (p[nOffset + 2 * i + k] - p[nOffset + 2 * j + k]);
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

	inline double computePhantomEntropy(const ParametersType& p, int n, int i) const {
        const int nOffset = n * m_nVars;
		ContinuousIndexType idx;
		idx[0] = p[nOffset + 2 * i];
		idx[1] = p[nOffset + 2 * i + 1];
		if (!m_DistanceMapInterpolators[n]->IsInsideBuffer(idx)) {
			return exp(5 * 5 / 9);
		}
		ImageType::PixelType dist = m_DistanceMapInterpolators[n]->EvaluateAtContinuousIndex(idx);
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
		} else if (dist > -3){
			return exp(-dist * dist / (m_Sigma * m_Sigma));
		}
        return 0;
	}

	void GetGuardedValueAndDerivative(const ParametersType & p,
                                      MeasureType & value, DerivativeType & derivative) const {
        // cout << "# points: " << m_nPoints << "; # vars: " << m_nVars << "; # subjects: " << m_nSubjects << endl;
        value = 0;

        MatrixType ensembleDeriv;
        MeasureType ensembleCost = 0;

        if (m_EnsembleFactor > 0) {
        	ComputeEnsembleEntropies(p, ensembleCost, ensembleDeriv);
            int k = 0;
            for (int i = 0; i < ensembleDeriv.rows(); i++) {
                for (int j = 0; j < ensembleDeriv.cols(); j++) {
                    derivative[k] = m_EnsembleFactor * ensembleDeriv[i][j];
                    k++;
                }
            }
        }

        for (int n = 0; n < m_nSubjects; n++) {
            const unsigned nOffset = m_nVars * n;
            MatrixType entropies(m_nPoints, m_nPoints + 1);
            double cost = 0;
            //#pragma omp parallel for
            for (int i = 0; i < m_nPoints; i++) {
                double sum = 0;
                for (int j = 0; j < m_nPoints; j++) {
                    if (i == j) {
                        continue;
                    }
                    ImageType::ValueType kappa = 1;
                    if (m_AdaptiveSampling) {
                        ContinuousIndexType idx;
                        idx[0] = p[nOffset + 2 * j];
                        idx[1] = p[nOffset + 2 * j + 1];
                        kappa = m_KappaMapInterpolators[n]->EvaluateAtContinuousIndex(idx);
                    }
                    entropies[i][j] = computeEntropy(p, nOffset, i, j, kappa);
                    sum += entropies[i][j];
                }
                entropies[i][m_nPoints] = computePhantomEntropy(p, n, i);
                sum += entropies[i][m_nPoints];
                if (sum > 0) {
                    for (int j = 0; j < entropies.cols(); j++) {
                        entropies[i][j] /= sum;
                    }
                }
                //#pragma omp critical
                {
                    cost += sum;
                }
            }
            //#pragma omp parallel for
            for (int i = 0; i < m_nPoints; i++) {
                for (int j = 0; j < m_nPoints; j++) {
                    if (i == j) {
                        continue;
                    }
                    for (int k = 0; k < 2; k++) {
                        derivative[nOffset + i * 2 + k] -= entropies[i][j] * (p[nOffset + 2 * i + k] - p[nOffset + 2 * j + k]);
                    }
                }
                ContinuousIndexType idx;
                idx[0] = p[nOffset + 2 * i];
                idx[1] = p[nOffset + 2 * i + 1];
                if (m_KappaMapInterpolators[n]->IsInsideBuffer(idx)) {
                    double dist = m_DistanceMapInterpolators[n]->EvaluateAtContinuousIndex(idx);
                    SliceType::IndexType gIdx;
                    gIdx[0] = ::round(p[nOffset + 2 * i]);
                    gIdx[1] = ::round(p[nOffset + 2 * i + 1]);
                    DistanceVectorImageType::PixelType gOffset = m_DistanceVectorMaps[n]->GetPixel(gIdx);
                    OffsetVectorType gradientPixel(gOffset.GetOffset(), DistanceVectorImageType::OffsetType::GetOffsetDimension());
                    OffsetVectorType normalizedGradient = gradientPixel.normalize();
                    if (dist > 0.5) {
                        for (int k = 0; k < 2; k++) {
                            derivative[nOffset + i * 2 + k] = - (entropies[i][m_nPoints] * normalizedGradient[k]);
                        }
                    }
                }
            }
            value += cost;
        }
	}

    void ComputeEnsembleEntropies(const ParametersType& p, MeasureType& ensembleCost, MatrixType& ensembleDeriv) const {
        MatrixType P(m_nSubjects, m_nVars);
        P.fill(0);
        for (int i = 0; i < m_nVars; i++) {
            for (int j = 0; j  < m_nSubjects; j++) {
                P[j][i] = p[j * m_nVars + i];
            }
        }
        VectorType meanP(P.cols());
        for (int i = 0; i < P.cols(); i++) {
            meanP[i] = P.get_column(i).mean();
        }

        // SxV matrix
        MatrixType Y(P);
        for (int j = 0; j < m_nSubjects; j++) {
            VectorType y = P.get_row(j) - meanP;
            Y.set_row(j, y);
        }

        // SxS matrix
        MatrixType YYt = Y * Y.transpose();

        vnl_symmetric_eigensystem<double> eigen(YYt);
        double sum = 0;
        for (int i = 0; i < YYt.cols(); i++) {
            sum += eigen.get_eigenvalue(i);
        }
        ensembleCost = sum;

        // S*V matrix
        ensembleDeriv = YYt * Y;
    }

	/** This method returns the value and derivative of the cost function corresponding
	 * to the specified parameters    */
	virtual void GetValueAndDerivative(const ParametersType & p,
                                       MeasureType & value, DerivativeType & derivative) const {
		derivative.SetSize(GetNumberOfParameters());
		derivative.Fill(0);
        GetGuardedValueAndDerivative(p, value, derivative);
	}

    void Print() {
        cout << "Sigma: " << m_Sigma << endl;
        cout << "Max Kappa: " << m_MaxKappa << endl;
        cout << "Cutoff Distance: " << m_CutoffDistance << endl;
        cout << "Cutoff Distance for boundary: " << m_PhantomCutoffDistance << endl;
        cout << "Adaptive Sampling: " << m_AdaptiveSampling << endl;
    }

protected:
	ImageEntropyCostFunction() {
		m_NumberOfPoints = 0;
		m_AdaptiveSampling = true;
        m_PhantomParticles = NULL;

        m_Sigma = 7.0;
        m_MaxKappa = sqrt(2);
        m_CutoffDistance = 15;
        m_PhantomCutoffDistance = 3;
        m_EnsembleFactor = 0.5;

        m_nSubjects = 0;
        m_nVars = 0;
        m_nPoints = 0;
        m_nTotalVars = 0;

	}

	virtual ~ImageEntropyCostFunction() {
	}

private:

	ParametersType m_SampleSigmas;
	int m_NumberOfPoints;
    // int m_NumberOfSubjects;

    // NumberOfSubjects * NumberOfPoints matrix
    ImageParticlesAlgorithm::ParametersList m_Points;

    int m_nSubjects;
    int m_nVars;
    int m_nPoints;
    int m_nTotalVars;

    ImageList m_KappaMaps;
    ImageList m_DistanceMaps;
    DistanceVectorList m_DistanceVectorMaps;
    InterpolatorList m_KappaMapInterpolators;
    InterpolatorList m_DistanceMapInterpolators;
    
	ImageParticlesAlgorithm::ParametersList* m_PhantomParticles;
	bool m_AdaptiveSampling;
	// bool m_UsePhantomParticles;
    
    double m_Sigma;
    double m_MaxKappa;
    double m_CutoffDistance;
    double m_PhantomCutoffDistance;
    double m_EnsembleFactor;
    double m_GradientSigma;
};


#endif
