/*
 * myImageParticlesAlgorithm.cpp
 *
 *  Created on: 2012. 11. 15.
 *      Author: joohwi
 */

#include "myImageParticlesAlgorithm.h"
#include "iostream"
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
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
#include "vtkPoints.h"
#include "itkOptimizerCommon.h"
#include "myEnsembleEntropy.h"
#include "vnlCommon.h"
#include "myParticleDynamics.h"


using namespace std;
const static int Dimensions = 2;

typedef itk::GradientRecursiveGaussianImageFilter<SliceType,ImageParticlesAlgorithm::GradientImageType> GradientImageFilter;
typedef itk::VectorMagnitudeImageFilter<ImageParticlesAlgorithm::GradientImageType,SliceType> VectorMagnitudeImageFilter;


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
			} else if (dynamic_cast<const FRPROptimizerType*>(caller) != NULL) {
				const FRPROptimizerType* opti = dynamic_cast<const FRPROptimizerType*>(caller);
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

ImageParticlesAlgorithm::ImageParticlesAlgorithm() : m_ImageList(NULL), m_EventCallback(NULL) {
	// TODO Auto-generated constructor stub
    m_ViewingDimension = 0;
    m_iters = 0;
    m_nSubjects = 0;
    m_nPoints = 0;
    m_nPoints = 0;
    m_nTotalParams = 0;
    m_Dim = 2;
    m_Slice = 0;
    m_ImageId = 0;
    m_Sigma = m_MaxKappa = m_CutoffDistance = m_EnsembleFactor = m_GradientSigma = 0;
}

ImageParticlesAlgorithm::~ImageParticlesAlgorithm() {
	// TODO Auto-generated destructor stub
}


// currently ignore; might not be necessary
void ImageParticlesAlgorithm::GetIndex(SliceInterpolatorType::ContinuousIndexType& idxOut, const OptimizerParametersType* params, const int subj, const int point) const {
    for (int k = 0; k < m_Dim; k++) {
        idxOut[k] = params->GetElement(subj*m_nParams + point*m_Dim + k);
    }
}



bool ImageParticlesAlgorithm::IsInsideBoundary(const OptimizerParametersType* params, int subj, int point) const {
    if (subj < m_nSubjects) {
        return false;
    }
    SliceInterpolatorType::ContinuousIndexType idx;
    GetIndex(idx, params, subj, point);
    return m_KappaMapInterpolators[subj]->IsInsideBuffer(idx) && m_Constraint.GetDistance(subj, idx) <= 0;
}

bool ImageParticlesAlgorithm::IsOutsideBoundary(const OptimizerParametersType* params, int subj, int point) const {
    if (subj < m_nSubjects) {
        return false;
    }
    SliceInterpolatorType::ContinuousIndexType idx;
    GetIndex(idx, params, subj, point);
    return m_KappaMapInterpolators[subj]->IsInsideBuffer(idx) && m_Constraint.GetDistance(subj, idx) > 0;
}

unsigned int ImageParticlesAlgorithm::GetNumberOfParameters() const {
    return m_nTotalParams;
}

/** This method returns the value of the cost function corresponding
 * to the specified parameters.    */
ImageParticlesAlgorithm::MeasureType ImageParticlesAlgorithm::GetValue(const ParametersType & parameters) const {
    MeasureType value;
    DerivativeType derivative;
    derivative.SetSize(GetNumberOfParameters());
    GetValueAndDerivative(parameters, value, derivative);
    return value;
}

/** This method returns the derivative of the cost function corresponding
 * to the specified parameters.   */
void ImageParticlesAlgorithm::GetDerivative(const ParametersType & parameters,
                                            DerivativeType & derivative) const {
    MeasureType value;
    GetValueAndDerivative(parameters, value, derivative);
}



double ImageParticlesAlgorithm::computeEntropy(const ParametersType& p, int nOffset, int i, int j,
                                               double kappa) const {

    double dist2 = 0;
    double si = m_Sigma / kappa;

    for (int k = 0; k < 2; k++) {
        dist2 += (p[nOffset + 2 * i + k] - p[nOffset + 2 * j + k])
        * (p[nOffset + 2 * i + k] - p[nOffset + 2 * j + k]);
    }

    // debug: entropy is always zero
    // m_Sigma was not set properly and dist2/(sigma)^2 was Inf (11/25)
    //    cout << "Dist2: " << dist2 << "; Sigma: " << si << "; Dist2/Sigma: " << (-dist2/(si*si)) << "; Cutoff: " << m_CutoffDistance << endl;
    //    cout << "Kappa: " << kappa << endl;

    if (dist2 > (m_CutoffDistance * m_CutoffDistance)) {
        return 0;
    } else {
        double force = exp(-dist2 / (si * si));
        // test: what happen if there's no repulsion term
        //        cout << "Force: " << force << endl;
        return force;
    }
}

double ImageParticlesAlgorithm::computePhantomEntropy(const ParametersType& p, int n, int i) const {
    const int nOffset = n * m_nParams;
    myImplicitSurfaceConstraint::ContinuousIndexType idx;
    idx[0] = p[nOffset + 2 * i];
    idx[1] = p[nOffset + 2 * i + 1];

    // debug: be sure if the point is out of the region or inside of the region
    if (!m_Constraint.IsInsideRegion(n, idx)) {
        return exp(5 * 5 / 9);
    }
    double dist = m_Constraint.GetDistance(n, idx);
    if (dist > 0) {
        return (dist + 1)*(dist+1);
    } else if (dist <= 0 && dist > -3) {
        return exp(-dist*dist/9);
    } else {
        return 0;
    }
}

void ImageParticlesAlgorithm::GetValueAndDerivative(const ParametersType & p,
                                                    MeasureType & value, DerivativeType & derivative) const {
    const bool useDomainConstraint = false;
    const bool useKappa = useDomainConstraint && true;

    // cout << "# points: " << m_nPoints << "; # vars: " << m_nVars << "; # subjects: " << m_nSubjects << endl;
    value = 0;

    DerivativeType ensembleDeriv(m_nSubjects * m_nParams);
    double ensembleCost = 0;

    bool applyEnsembleEntropy = false;
    if (applyEnsembleEntropy) {
        /**
         * This positional ensemble computes transformation parameters
         * to match the positions into the master example (the first case)
         */
        m_EnsembleEntropy->GetValueAndDerivative(p, ensembleCost, ensembleDeriv);
        value += ensembleCost;
        cout << "Ensemble Cost: " << ensembleCost << endl;
        cout << "Ensemble Factor: " << m_EnsembleFactor << endl;
        //        cout << "ensembleDeriv: " << ensembleDeriv << endl;
    }


//    
//    for (int n = 0; n < m_nSubjects; n++) {
//        for (int i = 0; i < m_nParams; i++) {
//            derivative[n * m_nParams + i] = ensembleDeriv[n][i];
//        }
//    }
    
    for (int n = 0; n < m_nSubjects; n++) {
        const unsigned nOffset = m_nParams * n;
        VNLMatrix entropies(m_nPoints, m_nPoints + 1);
        double cost = 0;
        //#pragma omp parallel for
        for (int i = 0; i < m_nPoints; i++) {
            double sum = 0;
            for (int j = 0; j < m_nPoints; j++) {
                if (i == j) {
                    continue;
                }
                ImageType::ValueType kappa = 1;
                if (useKappa) {
                    ContinuousIndexType idx;
                    idx[0] = p[nOffset + 2 * j];
                    idx[1] = p[nOffset + 2 * j + 1];
                    kappa = m_KappaMapInterpolators[n]->EvaluateAtContinuousIndex(idx);
                }
                entropies[i][j] = computeEntropy(p, nOffset, i, j, kappa);
                sum += entropies[i][j];
            }
            if (useDomainConstraint) {
                entropies[i][m_nPoints] = computePhantomEntropy(p, n, i);
            } else {
                entropies[i][m_nPoints] = 0;
            }
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

        int nBoundaryParticles = 0;
        
        //#pragma omp parallel for
        bool applySurfaceEntropy = (n == 0);
        if (applySurfaceEntropy) {
            for (int i = 0; i < m_nPoints; i++) {
                for (int j = 0; j < m_nPoints; j++) {
                    if (i == j) {
                        continue;
                    }
                    for (int k = 0; k < 2; k++) {
                        derivative[nOffset + i * 2 + k] -= entropies[i][j] * (p[nOffset + 2 * i + k] - p[nOffset + 2 * j + k]);
                    }
                }
            }
        }

        if (applyEnsembleEntropy && n > 0) {
            VNLVector deriv(&ensembleDeriv[nOffset], m_nParams);
            VNLVector normalizedDeriv = deriv.normalize();
            cout << normalizedDeriv << endl;
            for (int i = 0; i < m_nParams; i++) {
                derivative[nOffset + i] = (1 - m_EnsembleFactor) * derivative[nOffset + i] + (m_EnsembleFactor) * 10 * normalizedDeriv[i];
            }
        }

        bool applyBoundaryConstraint = false;
        if (applyBoundaryConstraint) {
            for (int i = 0; i < m_nPoints; i++) {
                SliceType::IndexType idx1;
                SliceInterpolatorType::ContinuousIndexType idx2;
                idx2[0] = p[nOffset + i * 2];
                idx2[1] = p[nOffset + i * 2 + 1];
                idx1[0] = ::round(idx2[0]);
                idx1[1] = ::round(idx2[1]);

                double boundaryDistance = m_Constraint.GetDistance(n, idx2);
                // debug: boundaryDistance should be less than zero for correct domain
                // cout << "Boundary Distance: " << boundaryDistance << endl;


                const bool boundsCheck = useDomainConstraint && m_Constraint.IsInsideRegion(n, idx2);
                if (!m_Constraint.IsInsideRegion(n, idx2)) {
                    for (int k = 0; k < m_Dim; k++) {
                        derivative[nOffset + i*2 + k] = 0;
                    }
                } else if (boundsCheck) {
                    if (boundaryDistance > 0) {
                        myImplicitSurfaceConstraint::DistanceVectorType offset = m_Constraint.GetOutsideOffset(n, idx1);
                        // what if the repulsion is not weighted; particles will move to boundary next step
                        //                    OffsetVectorType gradientPixel;
                        //                    for (int k = 0; k < m_Dim; k++) {
                        //                        gradientPixel[k] = offset[k];
                        //                    }
                        //                    OffsetVectorType normalizedGradient = gradientPixel.normalize();
                        // cout << normalizedGradient << endl;
                        for (int k = 0; k < m_Dim; k++) {
                            derivative[nOffset + i * 2 + k] = -offset[k];
                        }
                    } else if (boundaryDistance <= 0 && boundaryDistance > -1) {
                        // test: if there is no repulsion term
                        // debug: boundary check to prevent runtime exception
                        myImplicitSurfaceConstraint::DistanceVectorType offset = m_Constraint.GetInsideOffset(n, idx1);
                        OffsetVectorType gradientPixel;
                        gradientPixel[0] = offset[0];
                        gradientPixel[1] = offset[1];
                        OffsetVectorType normalizedGradient = gradientPixel.normalize();
                        //                    cout << "offset: " << offset.GetOffset()[0] << "," << offset.GetOffset()[1] << "; grad: " << gradientPixel << "; deriv:" << normalizedGradient << endl;
                        for (int k = 0; k < m_Dim; k++) {
                            derivative[nOffset + i * 2 + k] += entropies[i][m_nPoints] * normalizedGradient[k];
                        }
                    } else if (boundaryDistance < -1) {
                        myImplicitSurfaceConstraint::GradientPixelType grad = m_Constraint.GetGradient(n, idx1);
                        VectorType gradVector(grad.GetDataPointer(), m_Dim);
                        VectorType forceVector(&derivative[nOffset + i * 2], m_Dim);
                        double gf = dot_product(gradVector, forceVector);
                        if (gf > 0) {
                            VectorType resultVector = forceVector - gf * (gradVector);
                            for (int k = 0; k < m_Dim; k++) {
                                derivative[nOffset + i * 2 + k] = resultVector[k];
                            }
                            nBoundaryParticles ++;
                        }
                        // cout << boundaryDistance << "," << resultVector << "; subject: " << n << endl;
                    }
                }
                //            cout << "# of boundary particles: " << nBoundaryParticles << endl;
                
            }
        }
        

        value += cost;
    }
//    cout << "Value: " << value << endl;
//    cout << "Subjects: " << m_nSubjects << "; Points: " << m_nPoints << "; Params: " << m_nParams << "; Total Params: " << m_nTotalParams << endl;
}


void ImageParticlesAlgorithm::computeEnsembleEntropies(const ParametersType& p, MeasureType& ensembleCost, MatrixType& ensembleDeriv) const {
    MatrixType P(m_nSubjects, m_nParams);
    P.fill(0);
    for (int i = 0; i < m_nParams; i++) {
        for (int j = 0; j  < m_nSubjects; j++) {
            P[j][i] = p[j * m_nParams + i];
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
            for (int j = 0; j < Dims; j++) {
                initial[l*m_nParams+i*Dims+j] = indexes[i][j];
            }
        }
    }

    m_CurrentParams = initial;
}


void ImageParticlesAlgorithm::CreateInitialPoints(vtkPoints* pointSet) {
    m_nSubjects = m_ImageList->size();
    m_nParams = pointSet->GetNumberOfPoints()*Dims;
    m_nPoints = pointSet->GetNumberOfPoints();
    m_nTotalParams = m_nSubjects*m_nParams;

    // random pick up
    OptimizerParametersType initial;
    initial.SetSize(m_nTotalParams);
    for (int l = 0; l < m_ImageList->size(); l++) {
        for (int i = 0; i < m_nPoints; i++) {
            double *xyz = pointSet->GetPoint(i);
            initial[l*m_nParams+i*Dims+0] = xyz[0];
            initial[l*m_nParams+i*Dims+1] = xyz[1];
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
OptimizerType::Pointer ImageParticlesAlgorithm::CreateOptimizer() {
    if (m_Props.GetBool("useCGOpti", false)) {
        FRPROptimizerType::Pointer opti = FRPROptimizerType::New();
        opti->SetUseUnitLengthGradient(m_Props.GetBool("useUnitLengthGradient", true));
        opti->SetMaximumIteration(m_Props.GetInt("numberOfIterations", 100));
        return OptimizerType::Pointer(dynamic_cast<OptimizerType*>(opti.GetPointer()));
    } else if (m_Props.GetBool("useGDOpti", true)) {
        // gradient descent is the default optimizer
        GDOptimizerType::Pointer opti = GDOptimizerType::New();
        opti->SetNumberOfIterations(m_Props.GetInt("numberOfIterations", 100));
        return OptimizerType::Pointer(dynamic_cast<OptimizerType*>(opti.GetPointer()));
    } else if (m_Props.GetBool("useLBFGSOpti", false)) {

    }
    return OptimizerType::Pointer(NULL);
}


void ImageParticlesAlgorithm::SetPropertyAccess(PropertyAccess props) {
    m_Props = props;
    SetCutoffDistance(m_Props.GetDouble("cutoffDistance", 15));
    SetMaxKappa(m_Props.GetDouble("maxKappa", 2));
    SetPhantomCutoffDistance(m_Props.GetDouble("phantomCutoff", 3));
    SetEnsembleFactor(m_Props.GetDouble("ensembleFactor", 0));
    SetGradientSigma(m_Props.GetDouble("gradientSigmaX", 0.5));
    SetSigma(m_Props.GetDouble("sigma", 3));
}

void ImageParticlesAlgorithm::SetImageList(ImageContainer::List *list) {
    m_ImageList = list;
    m_nSubjects = list->size();
    m_KappaMaps.clear();
    m_KappaMapInterpolators.clear();
    m_Constraint.Clear();

    for (int i = 0; i < m_nSubjects; i++) {
        ImageContainer::Pointer image = list->at(i);
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
    }

    m_Constraint.SetImageList(list);
}

/**
 * Set up cost function and execute optimizer
 */
void ImageParticlesAlgorithm::RunOptimization() {
    m_iters = 0;
    m_Traces.clear();

    m_EnsembleEntropy = myEnsembleEntropy::New();
    m_EnsembleEntropy->SetImageList(m_ImageList);
    m_EnsembleEntropy->SetInitialPositions(m_CurrentParams, m_nSubjects, m_nPoints, m_nParams);
    // 3x3 image patch
    m_EnsembleEntropy->SetPatchSize(m_Props.GetInt("particlePatchSize", 3));
    m_EnsembleEntropy->SetTransformTypeToRigid();

    ContinueOptimization();
}

/**
 * Continue optimization process with current parameter values
 * current parameters can be the result of previous optimization
 * , and may be different from initial parameters
 */
void ImageParticlesAlgorithm::ContinueOptimization() {
    try {
        m_EnsembleEntropy->SetGradientScale(m_Props.GetDouble("ensembleGradientRatio", 0), m_Props.GetDouble("ensembleGradientScale", 0));
        ImageOptimizerProgress::Pointer progress = ImageOptimizerProgress::New();
        progress->SetReportCallback(this);

        cout << "Starting optimization..." << endl;
        cout << "# of params: " << m_CurrentParams.GetSize() << endl;
        OptimizerType::Pointer opti = CreateOptimizer();
        if (opti.IsNull()) {
            return;
        }
        OptimizerType::ScalesType scales;
        scales.SetSize(m_nTotalParams);
        scales.Fill(1);
        opti->SetScales(scales);
        opti->SetCostFunction(this);
        opti->SetInitialPosition(m_CurrentParams);
        opti->AddObserver(itk::IterationEvent(), progress);
        opti->StartOptimization();
        m_CurrentParams = opti->GetCurrentPosition();
    } catch (itk::ExceptionObject& e) {
        cout << e << endl;
    }
    cout << "Optimization done..." << endl;

}

void ImageParticlesAlgorithm::RunODE() {
    ParticleSystem system(m_nSubjects, m_nPoints);
    system.SetHistoryVector(&m_Traces);
    system.SetConstraint(&m_Constraint);
    system.SetPositions(&m_CurrentParams);
    system.Integrate();
    system.GetPositions(&m_CurrentParams);
}

void ImageParticlesAlgorithm::ContinueODE() {

}

void ImageParticlesAlgorithm::ReportParameters(const OptimizerParametersType &params, int iterNo, double cost) {
    m_Traces.push_back(params);
    double iterCost[] = { m_iters++, cost };
    m_EventCallback->EventRaised(0xADDCEC, 0, this, iterCost);
}

const VNLVector* ImageParticlesAlgorithm::GetTraceParameters(int idx) {
    if (idx < m_Traces.size()) {
        return &m_Traces[idx];
    } else {
        return NULL;
    }
}

