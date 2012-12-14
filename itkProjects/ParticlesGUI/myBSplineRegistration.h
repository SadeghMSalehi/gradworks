//
//  myBSplineRegistration.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 12/6/12.
//
//

#ifndef __ParticlesGUI__myBSplineRegistration__
#define __ParticlesGUI__myBSplineRegistration__

#include <iostream>
#include "vnlCommon.h"
#include "itkPointSet.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "myImageContainer.h"
#include "PropertyAccess.h"
#include "itkDisplacementFieldTransform.h"
#include "itkSingleValuedCostFunction.h"
#include "itkBSplineTransform.h"

namespace my {
    typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DisplacementFieldPointSetType, DisplacementFieldType> BSplineFilterType;
    typedef BSplineFilterType::WeightsContainerType WeightsContainerType;
    typedef itk::BSplineTransform<double, SDim, 3> BSplineTransform;
    
    // optimizer variables
    typedef double MeasureType;
    typedef itk::CostFunction::ParametersType ParametersType;
    typedef itk::SingleValuedCostFunction::DerivativeType DerivativeType;

    class BSplineRegistration;

    // cost function for deformable transformation optimizer
    class LandmarkMetric : public itk::SingleValuedCostFunction {
    public:
        itkTypeMacro(LandmarkMetric, itk::SingleValuedCostFunction);
        itkNewMacro(LandmarkMetric);
        
        virtual unsigned int GetNumberOfParameters() const;
        virtual MeasureType GetValue(const ParametersType& p) const;
        virtual void GetDerivative(const ParametersType& p, DerivativeType& d) const;
        virtual void GetValueAndDerivative(const ParametersType& p, MeasureType &b, DerivativeType& d) const;
        void SetContext(BSplineRegistration* ctx);
        
    protected:
        LandmarkMetric();
        virtual ~LandmarkMetric();
        
    private:
        LandmarkMetric(const LandmarkMetric& o);
        LandmarkMetric& operator=(const LandmarkMetric& o);

        MeasureType ComputeMSE(VNLVector& error) const;
        void ComputeDerivative(VNLVector& error, const ParametersType& p, DerivativeType& d) const;

        BSplineRegistration* m_Context;
        BSplineTransform::Pointer m_Transform;
    };

    // bspline registration
    //  1) scattered data approximation
    //  2) free-form deformation
    class BSplineRegistration {
        friend class LandmarkMetric;
    public:
        BSplineRegistration();
        ~BSplineRegistration();

        inline int GetNumberOfParams() { return m_nParams; }
        inline int GetNumberOfPoints() { return m_nPoints; }
        inline VNLVector& GetSourcePoints() { return m_Source; }
        inline VNLVector& GetTargetPoints() { return m_Target; }

        void SetUseFreeFormDeformation(bool ffd);
        void SetPropertyAccess(PropertyAccess props);
        void SetReferenceImage(SliceType::Pointer refImage);
        void SetLandmarks(int n, double* src, double* dst);
        void Update();
        void UpdateDeformation();
        void UpdateInterpolation();

        SliceType::Pointer WarpImage(SliceType::Pointer srcImage);
        SliceType::Pointer GetDeterminantOfJacobian();
        SliceType::Pointer GetDisplacementMagnitude();
        DisplacementFieldType::Pointer GetDisplacementField();
        DisplacementFieldType::Pointer GetControlPoints();
        FieldTransformType::Pointer GetTransform();
        SliceTransformType::Pointer GetRawTransform();

    private:
        bool m_UseFFD;
        VNLVector m_Params;
        SliceType::Pointer m_RefImage;
        DisplacementFieldPointSetType::Pointer m_FieldPoints;
        DisplacementFieldType::Pointer m_DisplacementField;
        DisplacementFieldType::Pointer m_PhiLattice;
        PropertyAccess m_Props;
        
        int m_nParams;
        int m_nPoints;
        VNLVector m_Source;
        VNLVector m_Target;
    };
}

#endif /* defined(__ParticlesGUI__myBSplineRegistration__) */