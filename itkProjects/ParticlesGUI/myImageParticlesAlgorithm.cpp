/*
 * myImageParticlesAlgorithm.cpp
 *
 *  Created on: 2012. 11. 15.
 *      Author: joohwi
 */

#include "myImageParticlesAlgorithm.h"
#include "iostream"
#include "vnl/vnl_matrix.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkFRPROptimizer.h"
#include "itkLBFGSOptimizer.h"
#include "itkImageIO.h"
#include "algorithm"

using namespace std;
const static int Dimensions = 2;

typedef itk::RegularStepGradientDescentOptimizer GDOptimizerType;
typedef itk::FRPROptimizer CGOptimizerType;
typedef itk::LBFGSOptimizer LBFGSOptimizerType;


class ImageOptimizerProgress: public itk::Command {
private:
	int m_Counter;

public:
	/** Standard class typedefs. */
	typedef ImageOptimizerProgress Self;
	typedef itk::Command Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
    typedef itk::SingleValuedNonLinearOptimizer OptimizerType;

	/** Run-time type information (and related methods). */
	itkTypeMacro(ImageOptimizerProgress, itk::Command)
	;

	itkNewMacro(ImageOptimizerProgress)
	;

	/** Abstract method that defines the action to be taken by the command. */
	virtual void Execute(itk::Object *caller, const itk::EventObject & event) {
		this->Execute((const itk::Object*) caller, event);
	}

	/** Abstract method that defines the action to be taken by the command.
	 * This variant is expected to be used when requests comes from a
	 * const Object */
	virtual void Execute(const itk::Object *caller, const itk::EventObject & event) {
		const OptimizerType* realCaller = dynamic_cast<const OptimizerType*>(caller);
		if (realCaller == NULL) {
			return;
		}
		if (++m_Counter % 1 == 0) {
			double value = 0;
			if (dynamic_cast<const LBFGSOptimizerType*>(caller) != NULL) {
				const LBFGSOptimizerType* opti = dynamic_cast<const LBFGSOptimizerType*>(caller);
                value = opti->GetCachedValue();
                //                g_costHistory.push_back(opti->GetCachedValue());
			} else if (dynamic_cast<const CGOptimizerType*>(caller) != NULL) {
				const CGOptimizerType* opti = dynamic_cast<const CGOptimizerType*>(caller);
                value = opti->GetValue();
                //                g_costHistory.push_back(opti->GetValue());
			} else if (dynamic_cast<const GDOptimizerType*>(caller) != NULL) {
				const GDOptimizerType* opti = dynamic_cast<const GDOptimizerType*>(caller);
                value = opti->GetValue();
                //                g_costHistory.push_back(opti->GetValue());
			}
            //            if (realCaller->GetCurrentPosition().GetSize() > 0) {
            //                m_ParametersHistory->push_back(realCaller->GetCurrentPosition());
            //            }
            cout << "Iteration: " << m_Counter << "; Cost: " << value << endl;
            if (m_Algo != NULL) {
                m_Algo->ReportParameters(realCaller->GetCurrentPosition(), m_Counter, value);
            }
		}
	}

    void SetReportCallback(ImageParticlesAlgorithm* algo) {
        m_Algo = algo;
    }

protected:
	ImageOptimizerProgress() {
        m_Counter = 0;
        m_Algo = NULL;
	}
	virtual ~ImageOptimizerProgress() {
	}
private:
	ImageOptimizerProgress(const Self &);        //purposely not implemented
	void operator=(const Self &); //purposely not implemented

    ImageParticlesAlgorithm* m_Algo;
};


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

    void SetNumberOfParameters(int m) {
        m_nVars = m;
    }

    void SetNumberOfPoints(int n) {
        m_nPoints = n;
    }

    void SetNumberOfSubjects(int n) {
        m_nSubjects = n;
    }

	virtual unsigned int GetNumberOfParameters() const {
		return m_nVars;
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
		} else {
			return exp(-dist * dist / (m_Sigma * m_Sigma));
		}
	}

	void GetGuardedValueAndDerivative(const ParametersType & p,
			MeasureType & value, DerivativeType & derivative) const {
        // cout << "# points: " << m_nPoints << "; # vars: " << m_nVars << "; # subjects: " << m_nSubjects << endl;
        value = 0;
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
                    if (dist > 0) {
                        SliceType::IndexType gIdx;
                        gIdx[0] = ::round(p[nOffset + 2 * i]);
                        gIdx[1] = ::round(p[nOffset + 2 * i + 1]);
                        OffsetVectorType gradientPixel(m_DistanceVectorMaps[n]->GetPixel(gIdx).GetOffset(), ImageType::OffsetType::GetOffsetDimension());
                        OffsetVectorType normalizedGradient = gradientPixel.normalize();
                        for (int k = 0; k < 2; k++) {
                            derivative[nOffset + i * 2 + k] = -(entropies[i][m_nPoints] * normalizedGradient[k]);
                            // cout << derivative[nOffset+i*2+k] << "[" << i << "; " << dist << "]";
                        }
                    } else if (dist < 0) {
                        SliceType::IndexType gIdx;
                        gIdx[0] = ::round(p[nOffset + 2 * i]);
                        gIdx[1] = ::round(p[nOffset + 2 * i + 1]);
                        OffsetVectorType gradientPixel(m_DistanceVectorMaps[n]->GetPixel(gIdx).GetOffset(), ImageType::OffsetType::GetOffsetDimension());
                        OffsetVectorType normalizedGradient = gradientPixel.normalize();
                        for (int k = 0; k < 2; k++) {
                            derivative[nOffset + i * 2 + k] -= (entropies[i][m_nPoints] * normalizedGradient[k]);
                        }
                    }
                }
            }
            value += cost;
        }

        if (m_EnsembleFactor > 0) {
        	MeasureType ensembleCost = 0;
        	const MatrixType ensembleDeriv = ComputeEnsembleEntropies(p, ensembleCost);
        	int k = 0;
        	for (int j = 0; j < (int) ensembleDeriv.rows(); j++) {
        		for (int i = 0; i < (int) ensembleDeriv.cols(); i++) {
        			derivative[k] = ((1 - m_EnsembleFactor) * derivative[k] + m_EnsembleFactor * ensembleDeriv[j][i]);
        			k++;
        		}
        	}

        	value += (m_EnsembleFactor * ensembleCost);
        }
	}

    MatrixType ComputeEnsembleEntropies(const ParametersType& p, MeasureType& entropy) const {
        MatrixType P;
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

        // S*V matrix
        MatrixType YYtY = YYt * Y;
        return YYtY;
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
	}

	virtual ~ImageEntropyCostFunction() {
	}

private:

	ParametersType m_SampleSigmas;
	int m_NumberOfPoints;
    int m_NumberOfSubjects;

    // NumberOfSubjects * NumberOfPoints matrix
    ImageParticlesAlgorithm::ParametersList m_Points;

    int m_nSubjects;
    int m_nVars;
    int m_nPoints;

    ImageList m_KappaMaps;
    ImageList m_DistanceMaps;
    DistanceVectorList m_DistanceVectorMaps;
    InterpolatorList m_KappaMapInterpolators;
    InterpolatorList m_DistanceMapInterpolators;

	ImageParticlesAlgorithm::ParametersList* m_PhantomParticles;
	bool m_AdaptiveSampling;
	bool m_UsePhantomParticles;

    double m_Sigma;
    double m_MaxKappa;
    double m_CutoffDistance;
    double m_PhantomCutoffDistance;
    double m_EnsembleFactor;
    double m_GradientSigma;
};

ImageParticlesAlgorithm::ImageParticlesAlgorithm() : m_ImageList(NULL), m_Running(false), m_EventCallback(NULL) {
	// TODO Auto-generated constructor stub
    m_ViewingDimension = 0;
}

ImageParticlesAlgorithm::~ImageParticlesAlgorithm() {
	// TODO Auto-generated destructor stub
}


/**
 * Assume that images are already set up by SetImageList
 * Actual subjects to be processed are determined by this initial points
 * The number of points should be maintained as same with previous initial points
 */
void ImageParticlesAlgorithm::AddInitialPoints(OptimizerParametersType& points) {
//    m_InitialPoints.push_back(points);
//    m_nSubjects = m_InitialPoints.size();
//    m_nParams = points.GetSize();
//    m_nPoints = m_nParams / Dims;
//    m_nTotalParams = m_nParams * m_nSubjects;
}

/**
 * Assuming the overlapping regions are large, randomly picks points inside the overlaps
 * First, compute the intersection across subjects
 * Second, list up every point in the domain
 * Third, randomly picks nPoints from the list
 *
 */
void ImageParticlesAlgorithm::CreateRandomInitialPoints(int nPoints) {
    if (m_ImageList == NULL || m_ImageList->size() == 0) {
        return;
    }

    m_nSubjects = m_ImageList->size();
    m_nParams = nPoints*Dims;
    m_nPoints = nPoints;
    m_nTotalParams = m_nSubjects*m_nParams;

    itkcmds::itkImageIO<LabelSliceType> io;
    LabelSliceType::Pointer intersection = io.NewImageT(m_ImageList->at(0)->GetLabelSlice(m_ViewingDimension));
    LabelSliceType::RegionType region = intersection->GetBufferedRegion();

    // compute intersection by looping over region
    std::vector<LabelSliceType::IndexType> indexes;
    LabelSliceIteratorType iter(intersection, region);
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
        LabelSliceType::IndexType idx = iter.GetIndex();
        LabelSliceType::PixelType pixel = 1;
        for (int i = 0; i < m_ImageList->size(); i++) {
            if (m_ImageList->at(i)->GetLabelSlice(m_ViewingDimension)->GetPixel(idx) == 0) {
                pixel = 0;
                break;
            }
        }
        if (pixel > 0) {
            iter.Set(pixel);
            indexes.push_back(idx);
        }
    }

    std::random_shuffle(indexes.begin(), indexes.end());

    // random pick up
    OptimizerParametersType initial;
    initial.SetSize(m_nTotalParams);
    std::random_shuffle(indexes.begin(), indexes.end());
    for (int l = 0; l < m_ImageList->size(); l++) {
        for (int i = 0; i < nPoints; i++) {
            int k = i * Dims;
            for (int j = 0; j < Dims; j++) {
                initial[k+j] = indexes[i][j];
            }
        }
    }

    m_CurrentParams = initial;
}


/**
 * Compute distance map for labels and set up initial parameters
 *
 */
void ImageParticlesAlgorithm::PrepareOptimization() {
    m_CurrentParams.SetSize(m_nTotalParams);
    for (int i = 0; i < m_nSubjects; i++) {
        int j = i * m_nParams;
        for (int k = 0; k < m_nParams; k++) {
            m_CurrentParams[j+k] = m_InitialPoints[i][k];
        }
    }
}

/**
 * Create the instance of the optimizer will be used and set up appropriate parameters from a user
 */
ImageParticlesAlgorithm::OptimizerType::Pointer ImageParticlesAlgorithm::CreateOptimizer() {
    if (m_Props.GetBool("useCGOpti", false)) {
        CGOptimizerType::Pointer opti = CGOptimizerType::New();
        opti->SetUseUnitLengthGradient(m_Props.GetBool("useUnitLengthGradient", true));
        opti->SetMaximumIteration(m_Props.GetInt("numberOfIterations", 100));
        return OptimizerType::Pointer(dynamic_cast<OptimizerType*>(opti.GetPointer()));
    } else if (m_Props.GetBool("useGDOpti", true)) {
        GDOptimizerType::Pointer opti = GDOptimizerType::New();
        opti->SetNumberOfIterations(m_Props.GetInt("numberOfIterations", 100));
        return OptimizerType::Pointer(dynamic_cast<OptimizerType*>(opti.GetPointer()));
    } else if (m_Props.GetBool("useLBFGSOpti", false)) {

    }
    return OptimizerType::Pointer(NULL);
}

/**
 * Set up cost function and execute optimizer
 */
void ImageParticlesAlgorithm::RunOptimization() {
    ImageEntropyCostFunction::Pointer costFunc = ImageEntropyCostFunction::New();
    m_CostFunc = CostFunctionType::Pointer(dynamic_cast<CostFunctionType*>(costFunc.GetPointer()));
    m_iters = 0;

    costFunc->SetCutoffDistance(m_Props.GetDouble("cutoffDistance", 15));
    costFunc->SetMaxKappa(m_Props.GetDouble("maxKappa", 2));
    costFunc->SetPhantomCutoffDistance(m_Props.GetDouble("phantomCutoff", 3));
    costFunc->SetEnsembleFactor(m_Props.GetDouble("ensembleFactor", 0));
    costFunc->SetGradientSigma(m_Props.GetDouble("gradientSigma", 0.5));


    for (int i = 0; i < m_nSubjects; i++) {
        costFunc->AddSubjects(m_ImageList->at(i));
    }
    costFunc->SetNumberOfParameters(m_CurrentParams.GetSize());
    costFunc->SetNumberOfSubjects(m_nSubjects);
    costFunc->SetNumberOfPoints(m_nPoints);

    m_Traces.clear();
    ContinueOptimization();
	m_Running = true;
}

/**
 * Continue optimization process with current parameter values
 * current parameters can be the result of previous optimization
 * , and may be different from initial parameters
 */
void ImageParticlesAlgorithm::ContinueOptimization() {
    try {
        ImageOptimizerProgress::Pointer progress = ImageOptimizerProgress::New();
        progress->SetReportCallback(this);

        cout << "Starting optimization..." << endl;
        cout << "# of params: " << m_CurrentParams.GetSize() << endl;
        OptimizerType::Pointer opti = CreateOptimizer();
        OptimizerType::ScalesType scales;
        scales.SetSize(m_nTotalParams);
        scales.Fill(1);
        opti->SetScales(scales);
        opti->SetCostFunction(m_CostFunc);
        opti->SetInitialPosition(m_CurrentParams);
        opti->AddObserver(itk::IterationEvent(), progress);
        opti->StartOptimization();
        m_CurrentParams = opti->GetCurrentPosition();
    } catch (itk::ExceptionObject& e) {
        cout << e << endl;
    }
    cout << "Optimization done..." << endl;

}

bool ImageParticlesAlgorithm::IsRunning() {
	return m_Running;
}

void ImageParticlesAlgorithm::ReportParameters(const OptimizerParametersType &params, int iterNo, double cost) {
    m_Traces.push_back(params);
    double iterCost[] = { m_iters++, cost };
    m_EventCallback->EventRaised(0xADDCEC, 0, this, iterCost);
}

const ImageParticlesAlgorithm::OptimizerParametersType* ImageParticlesAlgorithm::GetTraceParameters(int idx) {
    if (idx < m_Traces.size()) {
        return &m_Traces[idx];
    } else {
        return NULL;
    }
}