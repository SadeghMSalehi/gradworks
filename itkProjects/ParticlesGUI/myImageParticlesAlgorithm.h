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
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "PropertyAccess.h"

class vtkPoints;

class ImageParticlesAlgorithm: public itk::LightObject {
public:
    typedef itk::SingleValuedCostFunction CostFunctionType;
	typedef itk::SingleValuedNonLinearOptimizer OptimizerType;
	typedef OptimizerType::ParametersType OptimizerParametersType;
	typedef std::vector<OptimizerParametersType> ParametersList;

	typedef ImageParticlesAlgorithm Self;
	typedef itk::LightObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

    const static int Dims = 2;

	itkTypeMacro(ImageParticlesAlgorithm, itk::LightObject);
	itkNewMacro(ImageParticlesAlgorithm);

    void SetViewingDimension(int n) { m_ViewingDimension = n; }
    int GetNumberOfSubjects() const { return m_nSubjects; }
    int GetNumberOfParams() const { return m_nParams; }
    int GetNumberOfPoints() const { return m_nPoints; }
    const OptimizerParametersType& GetCurrentParams() const { return m_CurrentParams; }
    
    // mandatory methods to invoke before optimization
    void SetPropertyAccess(PropertyAccess props) { m_Props = props; }
	void SetImageList(ImageContainer::List* list) { m_ImageList = list; }
	void AddInitialPoints(OptimizerParametersType& points);
    void CreateRandomInitialPoints(int nPoints);

    // optimization execution
	void RunOptimization();
	void ContinueOptimization();
	bool IsRunning();
    void ReportParameters(const OptimizerParametersType& params);
    const OptimizerParametersType* GetTraceParameters(int idx);
    int GetNumberOfTraces() { return m_Traces.size(); }
    
protected:
	ImageParticlesAlgorithm();
	virtual ~ImageParticlesAlgorithm();

	void PrepareOptimization();
    OptimizerType::Pointer CreateOptimizer();

private:
	ImageParticlesAlgorithm(const Self &);
	void operator=(const Self &);

    int m_nSubjects;
    int m_nPoints;
    int m_nParams;
    int m_nTotalParams;

    CostFunctionType::Pointer m_CostFunc;
	ImageContainer::List* m_ImageList;
	ParametersList m_InitialPoints;
    OptimizerParametersType m_CurrentParams;
    PropertyAccess m_Props;
	bool m_Running;

    // ui related field
    int m_ViewingDimension;
    ParametersList m_Traces;
};

#endif /* MYIMAGEPARTICLESALGORITHM_H_ */
