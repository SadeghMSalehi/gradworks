//
//  myParticleAlgorithm.h
//  laplacePDE
//
//  Created by Joohwi Lee on 11/14/12.
//
//

#ifndef __laplacePDE__myParticleAlgorithm__
#define __laplacePDE__myParticleAlgorithm__

#include <iostream>
#include "itkLightObject.h"
#include "myImageContainer.h"
#include "PropertyAccess.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkSingleValuedNonLinearOptimizer.h"

typedef itk::SingleValuedNonLinearOptimizer OptimizerType;
typedef itk::SingleValuedNonLinearOptimizer::ParametersType OptimizerParametersType;
typedef itk::SignedDanielssonDistanceMapImageFilter<LabelType, ImageType> DistanceMapFilter;
typedef DistanceMapFilter::VectorImageType DistanceVectorImageType;
typedef std::vector<OptimizerParametersType> ParametersVectorType;

class vtkPolyData;
class vtkPointSet;

class ParticleAlgorithm: public itk::LightObject {
public:
    typedef ParticleAlgorithm Self;
    typedef itk::LightObject Superclass;
    typedef itk::SmartPointer<ParticleAlgorithm> Pointer;
    typedef itk::SmartPointer<const ParticleAlgorithm> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(ParticleAlgorithm, itk::LightObject);

    inline const int GetDimensions() const {
        return m_Dims;
    }

    inline const int GetNumberOfPoints() const {
        return m_NumberOfPoints;
    }

    ParametersVectorType GetParameterTrace() {
        return m_ParameterHistory;
    }

    /**
     * Define where particles are constrained
     */
    void SetImplicitSurface(LabelType::Pointer image);

    /**
     * Provide initial particles' coordinates
     */
    void SetInitialParticles(vtkPointSet* poly);

    void RunOptimization();
    void ContinueOptimization();
    vtkPointSet* GetResultPoints() const;
    void ReportParameters(const OptimizerParametersType& params);
    inline void SetPropertyAccess(PropertyAccess prop) { m_Props = prop; }

protected:
    ParticleAlgorithm() {
    	m_Points = NULL;
    	m_ResultPoints = NULL;
    	m_NumberOfPoints = 0;
    };
    virtual ~ParticleAlgorithm() {};


private:
    ParticleAlgorithm(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented

    LabelType::Pointer m_ImplicitSurface;
    ImageType::Pointer m_DistanceValueMap;
    DistanceVectorImageType::Pointer m_DistanceVectorMap;
    vtkPointSet* m_Points;
    vtkPointSet* m_ResultPoints;
    int m_NumberOfPoints;
    const static int m_Dims = 3;
    OptimizerParametersType m_InitialParameters;
    ParametersVectorType m_ParameterHistory;
    PropertyAccess m_Props;
};

#endif /* defined(__laplacePDE__myParticleAlgorithm__) */
