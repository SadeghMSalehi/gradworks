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

class vtkPoints;

class ImageParticlesAlgorithm: public itk::SingleValuedCostFunction {
public:
	typedef std::vector<OptimizerParametersType> ParametersList;
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

    typedef itk::CovariantVector<double,2> GradientType;
    typedef itk::Image<GradientType,2> GradientImageType;
	typedef itk::LinearInterpolateImageFunction<SliceType,float> InterpolatorType;
	typedef InterpolatorType::ContinuousIndexType ContinuousIndexType;
	typedef std::vector<SliceType::Pointer> ImageList;
	typedef std::vector<InterpolatorType::Pointer> InterpolatorList;
    

    
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
    
    // number of subjects
    int GetNumberOfSubjects() const { return m_nSubjects; }
    // number of params per subject
    int GetNumberOfParams() const { return m_nParams; }
    // number of points per subject
    int GetNumberOfPoints() const { return m_nPoints; }
    
    // current parameters for optimization
    const OptimizerParametersType& GetCurrentParams() const { return m_CurrentParams; }
    
    // mandatory methods to invoke before optimization
    void SetPropertyAccess(PropertyAccess props);
	void SetImageList(ImageContainer::List* list);
	void AddInitialPoints(OptimizerParametersType& points);
    void CreateRandomInitialPoints(int nPoints);
    void CreateInitialPoints(vtkPoints* pointSet);

    // optimization execution
	void RunOptimization();
	void ContinueOptimization();
    
    void ReportParameters(const OptimizerParametersType& params, int iterNo, double cost);
    const OptimizerParametersType* GetTraceParameters(int idx);
    int GetNumberOfTraces() { return m_Traces.size(); }
    void SetEventCallback(EventCallback* callback) { m_EventCallback = callback; }


    void SetCurrentSliceAndView(int view, int slice) { m_ViewingDimension = view; m_Slice = slice; }
    bool IsCurrentSliceAndView(int view, int slice) { return m_ViewingDimension == view && m_Slice == slice; }

    void ProbeDerivedImage(const char* fname);
    
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

	ParametersList m_InitialPoints;
    OptimizerParametersType m_CurrentParams;
    PropertyAccess m_Props;

	ImageContainer::List* m_ImageList;
    myImplicitSurfaceConstraint m_Constraint;
    myEnsembleEntropy::Pointer m_EnsembleEntropy;

    EventCallback* m_EventCallback;

    // ui related field
    ParametersList m_Traces;

    ImageList m_KappaMaps;
    InterpolatorList m_KappaMapInterpolators;
    
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
