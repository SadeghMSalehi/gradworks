/*
 * myImageParticlesAlgorithm.h
 *
 *  Created on: 2012. 11. 15.
 *      Author: joohwi
 */

#ifndef MYIMAGEPARTICLESALGORITHM_H_
#define MYIMAGEPARTICLESALGORITHM_H_

#include <itkLightObject.h>
#include <itkSmartPointer.h>
#include <itkObjectFactory.h>
#include "myImageContainer.h"
#include "myImplicitSurfaceConstraint.h"
#include "myEventCallback.h"
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "PropertyAccess.h"
#include "itkOptimizerCommon.h"
#include "itkCommand.h"
#include "itkArray.h"
#include "myEnsembleEntropy.h"
#include "myBSplineRegistration.h"

class vtkPoints;

class ImageParticlesAlgorithm: public itk::SingleValuedCostFunction {
public:
	typedef ImageParticlesAlgorithm Self;
	typedef itk::SingleValuedCostFunction Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

    typedef double MeasureType;
    typedef Superclass::ParametersValueType ParametersValueType;
	typedef itk::Array<ParametersValueType> DerivativeType;
    // debug: type casting eleminate 0-1 double values to zero!!
    typedef vnl_vector_fixed<double, 2> OffsetVectorType;
    typedef vnl_vector<double> VectorType;
	typedef vnl_matrix<double> MatrixType;

	typedef std::vector<SliceType::Pointer> ImageList;
	typedef std::vector<SliceInterpolatorType::Pointer> InterpolatorList;
    typedef std::vector<VectorInterpolatorType::Pointer> GradientInterpolatorList;
    

    
    const static int Dims = 2;

	itkTypeMacro(ImageParticlesAlgorithm, itk::SingleValuedCostFunction);
	itkNewMacro(ImageParticlesAlgorithm);
    
    /**
     * Setters and Getters
     *
     */
    itkSetMacro(CutoffDistance, double);
    itkSetMacro(Sigma, double);
    itkSetMacro(PhantomCutoffDistance, double);
    itkSetMacro(MaxKappa, double);
    itkSetMacro(EnsembleFactor, double);
    itkSetMacro(GradientSigma, double);
    
    // Methods for Cost Function
    // total parameters for optimizer
    virtual unsigned int GetNumberOfParameters() const;
    virtual MeasureType GetValue(const ParametersType & parameters) const;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const;
	virtual void GetValueAndDerivative(const ParametersType & p,
                                       MeasureType & value, DerivativeType & derivative) const;

    double computeEntropy(const ParametersType& p, int nOffset, int i, int j, double kappa = 1) const;
    double computePhantomEntropy(const ParametersType& p, int n, int i) const;
    void computeEnsembleEntropies(const ParametersType& p, MeasureType& ensembleCost, MatrixType& ensembleDeriv) const;
    
    // connection information between imageList
    void SetViewingDimension(int n) { m_ViewingDimension = n; }

    void SetNumberOfParams(int p) { m_nParams = p; }
    void SetNumberOfSubjects(int n) { m_nSubjects = n; }
    void SetNumberOfPoints(int p) { m_nPoints = p; }
    // number of subjects
    int GetNumberOfSubjects() const { return m_nSubjects; }
    // number of params per subject
    int GetNumberOfParams() const { return m_nParams; }
    // number of points per subject
    int GetNumberOfPoints() const { return m_nPoints; }
    
    // current parameters for optimization
    const OptimizerParametersType& GetCurrentParams() const { return m_CurrentParams; }

    // set current parameters from outside
    void SetCurrentParams(VNLVector& params);
    
    // mandatory methods to invoke before optimization
    void SetPropertyAccess(PropertyAccess props);
	void SetImageList(ImageContainer::List* list);
	void AddInitialPoints(OptimizerParametersType& points);
    void CreateRandomInitialPoints(int nPoints);
    void CreateInitialPoints(vtkPoints* pointSet);
    void CreateUniformInitialization();

    // optimization execution
	void RunOptimization();
	void ContinueOptimization();
    
    void ReportParameters(const OptimizerParametersType& params, int iterNo, double cost);
    VNLVectorArray* GetTraceVector();
    void AddTrace(VNLVector params);
    const VNLVector* GetTraceParameters(int idx);
    const VNLVector* GetDensityTraces(int idx);
    int GetNumberOfTraces() { return m_Traces.size(); }
    void SetEventCallback(EventCallback* callback) { m_EventCallback = callback; }


    void SetCurrentSliceAndView(int view, int slice) { m_ViewingDimension = view; m_Slice = slice; }
    bool IsCurrentSliceAndView(int view, int slice) { return m_ViewingDimension == view && m_Slice == slice; }

    void ProbeDerivedImage(const char* fname);

    /**
     * Particle and image association functions
     */
    void GetIndex(SliceInterpolatorType::ContinuousIndexType& idxOut, const OptimizerParametersType* params, int nsubj, int npoint) const;
    bool IsInsideBoundary(const OptimizerParametersType* params, int subj, int point) const;
    bool IsOutsideBoundary(const OptimizerParametersType* params, int subj, int point) const;

    /**
     * ODE-based particle solution
     */
    void RunODE();
    void ContinueODE();
    PropertyAccess GetProperty();

    InterpolatorList* GetAttributeInterpolators() {
        return &m_KappaMapInterpolators;
    }

    InterpolatorList* GetImageInterpolators() {
        return &m_ImageInterpolators;
    }

    GradientInterpolatorList* GetGradientInterpolators() {
        return &m_GradientInterpolators;
    }

    // return image container of n'th subject
    ImageContainer::Pointer GetImage(int n) {
        return m_ImageList->at(n);
    }

    // Apply TPS transform
    void ApplyTPSorEBSTransform(int type);
    void ApplyBSplineTransform();
    DisplacementFieldType::Pointer GetDisplacementField();

    // return bspline registration for debugging and visualization
    my::BSplineRegistration* GetBSplineRegistration() {
        return &m_BSplineRegistration;
    }

    my::ImplicitSurfaceConstraint* GetConstraint() { return &m_Constraint; }

    /**
     * UI Event handling
     */
    void OnClick(double x, double y, int imageIdx);

    
protected:
	ImageParticlesAlgorithm();
	virtual ~ImageParticlesAlgorithm();

	void PrepareOptimization();
    OptimizerType::Pointer CreateOptimizer();

private:
	ImageParticlesAlgorithm(const Self &);
	void operator=(const Self &);

    int m_iters;
    int m_nSubjects;
    int m_nPoints;
    int m_nParams;
    int m_nTotalParams;
    
    int m_Dim;
    int m_Slice;
    int m_ImageId;
    int m_ViewingDimension;

	VNLVectorArray m_InitialPoints;
    OptimizerParametersType m_CurrentParams;

    PropertyAccess m_Props;

	ImageContainer::List* m_ImageList;
    my::ImplicitSurfaceConstraint m_Constraint;
    myEnsembleEntropy::Pointer m_EnsembleEntropy;

    EventCallback* m_EventCallback;

    // ui related field
    VNLVectorArray m_Traces;
    VNLVectorArray m_DensityTraces;
    STDDoubleArray m_CostTraces;

    ImageList m_KappaMaps;
    InterpolatorList m_KappaMapInterpolators;

    LabelSliceType::Pointer m_Intersection;

    InterpolatorList m_ImageInterpolators;
    GradientInterpolatorList m_GradientInterpolators;


    my::BSplineRegistration m_BSplineRegistration;
    
    DisplacementFieldType::Pointer m_BSplineDisplacementField;
    
    /**
     * Cost function parameters
     */
    double m_Sigma;
    double m_MaxKappa;
    double m_CutoffDistance;
    double m_PhantomCutoffDistance;
    double m_EnsembleFactor;
    double m_GradientSigma;
};

#endif /* MYIMAGEPARTICLESALGORITHM_H_ */
