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
    void SetPropertyAccess(PropertyAccess props) { m_Props = props; }
	void SetImageList(ImageContainer::List* list) { m_ImageList = list; }
	void AddInitialPoints(OptimizerParametersType& points);
    void CreateRandomInitialPoints(int nPoints);
    void CreateInitialPoints(vtkPoints* pointSet);

    // optimization execution
	void RunOptimization();
	void ContinueOptimization();
    
	bool IsRunning();
    void ReportParameters(const OptimizerParametersType& params, int iterNo, double cost);
    const OptimizerParametersType* GetTraceParameters(int idx);
    int GetNumberOfTraces() { return m_Traces.size(); }
    void SetEventCallback(EventCallback* callback) { m_EventCallback = callback; }
    
    
    // i don't like this marker but this is required to correctly render particles
    // particles are shown for each image and determined by viewing plane and slice index
    void SetSliceMarker(int slice) {
        m_Slice = slice;
        // m_ImageId = imageId;
    }
    
    bool IsMarkerValid(int dir, int slice) {
        return m_ViewingDimension == dir && m_Slice == slice;
    }

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

	ParametersList m_InitialPoints;
    OptimizerParametersType m_CurrentParams;
    PropertyAccess m_Props;

	ImageContainer::List* m_ImageList;
    myImplicitSurfaceConstraint m_Constraint;

	bool m_Running;
    EventCallback* m_EventCallback;

    // ui related field
    int m_ViewingDimension;
    ParametersList m_Traces;
    
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
