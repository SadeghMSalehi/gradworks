//
//  myEnsembleEntropy.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/24/12.
//
//

#ifndef __ParticlesGUI__myEnsembleEntropy__
#define __ParticlesGUI__myEnsembleEntropy__

#include <iostream>

#include "myImageContainer.h"
#include "itkOptimizerCommon.h"
#include "vnlCommon.h"
#include "itkSingleValuedCostFunction.h"

class myEnsembleEntropy : public itk::SingleValuedCostFunction {
public:
    typedef myEnsembleEntropy Self;
	typedef itk::SingleValuedCostFunction Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
    
    typedef double MeasureType;
    typedef Superclass::ParametersValueType ParametersValueType;
	typedef itk::Array<ParametersValueType> DerivativeType;
    
    itkTypeMacro(myEnsembleEntropy, itk::SingleValuedCostFunction);
    itkNewMacro(myEnsembleEntropy);
    
    virtual unsigned int GetNumberOfParameters() const;
    virtual MeasureType GetValue(const ParametersType & parameters) const;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const;
	virtual void GetValueAndDerivative(const ParametersType & p,
                                       MeasureType & value, DerivativeType & derivative) const;
    
    void SetImageList(ImageContainer::List* imageList);
    void SetInitialPositions(const OptimizerParametersType& params, int nSubject, int nPoints, int nParams);
    // void GetValueAndDerivative(const OptimizerParametersType& params, double& cost, VNLMatrix& deriv) const;
    
    void SetPatchSize(int n) { m_PatchSize = n; }
    void SetTransformType(int t) { m_TransformType = t; }
    void SetTransformTypeToRigid() { m_TransformType = 1; }
    void EstimateRigidParameters(VNLMatrix& transformParams, const OptimizerParametersType& params, int target, int source) const;
    void SetVariableCounts(int s, int p, int np);
    
protected:
    myEnsembleEntropy();
    virtual ~myEnsembleEntropy();
private:
    myEnsembleEntropy(const myEnsembleEntropy&);
    void operator=(const myEnsembleEntropy&);
    
    //void EstimateRigidParameters(VNLMatrix& transformParams, const OptimizerParametersType& params, int target, int source) const;
    
    ImageContainer::List* m_ImageList;
    int m_PatchSize;
    int m_TransformType;
    int m_nSubjects;
    int m_nPoints;
    int m_nParams;
    
    // assume there's only rigid transformations
    const int m_nTransformParams = 3;

};

#endif /* defined(__ParticlesGUI__myEnsembleEntropy__) */
