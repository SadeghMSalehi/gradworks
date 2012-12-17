//
//  BSplineRegistration.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 12/6/12.
//
//

#include "myBSplineRegistration.h"
#include "QElapsedTimer"
#include "itkWarpImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkImageIO.h"
#include "itkLBFGSOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkFRPROptimizer.h"
#include "iostream"

// estimation of displacement field via particle correspondence

using namespace std;

namespace my {
    typedef itk::WarpImageFilter<SliceType, SliceType, DisplacementFieldType> WarpImageFilterType;
    
    LandmarkMetric::LandmarkMetric() {
        
    }
    
    LandmarkMetric::~LandmarkMetric() {
        
    }

    MeasureType LandmarkMetric::ComputeMSE(VNLVector& error) const {
        return error.squared_magnitude();
    }

    void LandmarkMetric::ComputeDerivative(VNLVector& error, VNLVector& tX, VNLVector& Y, const ParametersType& p, DerivativeType& d) const {
        d.SetSize(m_nParams);
        d.Fill(0);

        BSplineTransform::JacobianType jac;
        BSplineTransform::InputPointType point;
        for (int i = 0; i < error.size(); i += SDim) {
            point[0] = tX[i];
            point[1] = tX[i+1];
            m_Transform->ComputeJacobianWithRespectToParameters(point, jac);
            for (int j = 0; j < m_nParams; j++) {
                d[j] += 2.0 * (jac[0][j] * error[i] + jac[1][j] * error[i+1]);
            }
        }
    }

    void LandmarkMetric::TransformPoints(VNLVector& x, VNLVector& y) const {
        const int nSize = x.size();
        y.set_size(nSize);
        for (int i = 0; i < nSize; i += SDim) {
            BSplineTransform::InputPointType p;
            p[0] = x[i];
            p[1] = x[i+1];
            BSplineTransform::OutputPointType q = m_Transform->TransformPoint(p);
            y[i] = q[0];
            y[i+1] = q[1];
        }
    }

    void LandmarkMetric::SetContext(BSplineRegistration* context) {
        m_Context = context;
    }

    void LandmarkMetric::SetTransform(BSplineTransform* transform) {
        m_Transform = transform;
        m_nParams = transform->GetNumberOfParameters();
    }

    unsigned int LandmarkMetric::GetNumberOfParameters() const {
        return m_nParams;
    }

    MeasureType LandmarkMetric::GetValue(const ParametersType& p) const {
        MeasureType v = 0;
        m_Transform->SetParameters(p);
        TransformPoints(m_Context->GetSourcePoints(), m_Tx);
        VNLVector error = m_Tx - m_Context->GetTargetPoints();
        v = ComputeMSE(error);
        
        return v;
    }
    
    void LandmarkMetric::GetDerivative(const ParametersType &p, DerivativeType &d) const {
        m_Transform->SetParameters(p);
        TransformPoints(m_Context->GetSourcePoints(), m_Tx);
        VNLVector& Y = m_Context->GetTargetPoints();
        VNLVector error = m_Tx - Y;
        ComputeDerivative(error, m_Tx, Y, p, d);
    }
    
    void LandmarkMetric::GetValueAndDerivative(const ParametersType &p, MeasureType &b, DerivativeType &d) const {
        m_Transform->SetParameters(p);
        TransformPoints(m_Context->GetSourcePoints(), m_Tx);
        VNLVector& Y = m_Context->GetTargetPoints();
        VNLVector error = m_Tx - m_Context->GetTargetPoints();
//        cout << ">> " << m_Tx << " >> " << m_Context->GetTargetPoints() << ", " << error << endl;
        b = ComputeMSE(error);
        ComputeDerivative(error, m_Tx, Y, p, d);

        // when using two control points, the result is not stable
//        cout << "Value: " << b << "; " << p << "; " << d << endl;
    }

    BSplineRegistration::BSplineRegistration() {
        m_UseFFD = false;
    }

    BSplineRegistration::~BSplineRegistration() {

    }
    
    void BSplineRegistration::SetUseFreeFormDeformation(bool ffd) {
        m_UseFFD = ffd;
    }

    void BSplineRegistration::SetPropertyAccess(PropertyAccess props) {
        m_Props = props;
    }


    void BSplineRegistration::SetReferenceImage(SliceType::Pointer refImage) {
        m_RefImage = refImage;
    }

    void BSplineRegistration::SetLandmarks(int n, double *src, double *dst) {
        m_nPoints = n;
        m_nParams = n * SDim;
        
        // store landmarks
        m_Source.set_size(m_nParams);
        m_Source.set(src);
        
        m_Target.set_size(m_nParams);
        m_Target.set(dst);
    }
    
    
    void BSplineRegistration::Update() {
        if (m_UseFFD) {
            UpdateDeformation();
        } else {
            UpdateInterpolation();
        }
        
    }

    //
    // free form deformation based on landmarks
    // the distance between landmark is optimizer via mean squared error
    //
    void BSplineRegistration::UpdateDeformation() {
        int nSize = m_Props.GetInt("numberOfControlPoints", 25);

        SliceType::SizeType refSize = m_RefImage->GetBufferedRegion().GetSize();
        SliceType::SpacingType refSpacing = m_RefImage->GetSpacing();

        BSplineTransform::PhysicalDimensionsType physicalDimensions;
        BSplineTransform::MeshSizeType meshSize;
        for (unsigned int i = 0; i < SDim; i++) {
            physicalDimensions[i] = refSpacing[i] * (refSize[i] - 1);
        }
        meshSize.Fill(nSize);


        //
        // control point domain describes locations of control points
        // they usually follows the property of the reference image
        //
        BSplineTransform::Pointer bspliner = BSplineTransform::New();
        bspliner->SetTransformDomainOrigin(m_RefImage->GetOrigin());
        bspliner->SetTransformDomainDirection(m_RefImage->GetDirection());
        bspliner->SetTransformDomainPhysicalDimensions(physicalDimensions);
        bspliner->SetTransformDomainMeshSize(meshSize);

        const int nBSplineParams = bspliner->GetNumberOfParameters();

        // debug: check if number of parameters are same as the number of control points * Dimension (meshSize)

        ParametersType params(nBSplineParams);
        params.Fill(0);
        bspliner->SetParameters(params);

        LandmarkMetric::Pointer metric = LandmarkMetric::New();
        metric->SetTransform(bspliner);
        metric->SetContext(this);

//        typedef itk::LBFGSOptimizer OptimizerType;
//        typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
        typedef itk::FRPROptimizer OptimizerType;

        OptimizerType::Pointer opti = OptimizerType::New();

        //  for LBFGS optimizer
        //        opti->SetGradientConvergenceTolerance(0.0005);
        //        opti->SetLineSearchAccuracy(0.9);
        //        opti->SetDefaultStepLength(.25);
        //        opti->SetMaximumNumberOfFunctionEvaluations(1000);
        //        opti->TraceOn();

        // for gradient descent
//        opti->SetMaximumStepLength(0.1);
//        opti->SetMinimumStepLength(0.01);
//        opti->SetNumberOfIterations(1000);
//        opti->SetGradientMagnitudeTolerance(1e-4);
//        opti->SetMinimize(true);

        // for FRPR optimizer
        opti->SetMaximumIteration(100);
        opti->SetMaximumLineIteration(100);
        opti->SetStepLength(0.1);
        opti->SetStepTolerance(0.01);
        opti->SetValueTolerance(1e-5);
        
        opti->SetCostFunction(metric);
        opti->SetInitialPosition(params);
        opti->StartOptimization();
        cout << "Stop Condition: " << opti->GetStopConditionDescription() << endl;

        // bspline transform might not have correct parameters
        bspliner->SetParameters(opti->GetCurrentPosition());
        m_FFDTransform = bspliner;
    }

    void BSplineRegistration::UpdateInterpolation() {        
        if (m_FieldPoints.IsNull()) {
            m_FieldPoints = DisplacementFieldPointSetType::New();
        }
        m_FieldPoints->Initialize();
        
        // create point structures
        PointSetType::Pointer srcPoints = PointSetType::New();
        PointSetType::Pointer dstPoints = PointSetType::New();
        
        srcPoints->Initialize();
        dstPoints->Initialize();
        
        int n = m_nPoints;
        double* pSrc = m_Source.data_block();
        double* pDst = m_Target.data_block();
        for (int i = 0; i < n; i++) {
            PointSetType::PointType iPoint;
            iPoint[0] = pSrc[0];
            iPoint[1] = pSrc[1];
            
            VectorType vector;
            vector[0] = pDst[0] - pSrc[0];
            vector[1] = pDst[1] - pSrc[1];
            m_FieldPoints->SetPoint(i, iPoint);
            m_FieldPoints->SetPointData(i, vector);
            pSrc += 2;
            pDst += 2;
        }

        int splineOrder = m_Props.GetInt("splineOrder", 3);
        int numOfLevels = m_Props.GetInt("numLevels", 1);
        int nSize = m_Props.GetInt("numberOfControlPoints", m_NumberOfControlPoints);

        BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
        BSplineFilterType::ArrayType numControlPoints;
        numControlPoints.Fill(nSize + splineOrder);

        SliceType::SizeType imageSize = m_RefImage->GetBufferedRegion().GetSize();
        SliceType::SpacingType imageSpacing = m_RefImage->GetSpacing();
        SliceType::PointType imageOrigin = m_RefImage->GetOrigin();

//        cout << "Image Size: " << imageSize << endl;
//        cout << "Image Spacing: " << imageSpacing << endl;
//        cout << "Image Origin: " << imageOrigin << endl;
//        cout << "# control points: " << numControlPoints << endl;

        try {
            // debug: reparameterized point component is outside
            QElapsedTimer timer;
            timer.start();
            bspliner->SetOrigin(imageOrigin);
            bspliner->SetSpacing(imageSpacing);
            bspliner->SetSize(imageSize);
            bspliner->SetGenerateOutputImage(true);
            bspliner->SetNumberOfLevels(numOfLevels);
            bspliner->SetSplineOrder(splineOrder);
            bspliner->SetNumberOfControlPoints(numControlPoints);
            bspliner->SetInput(m_FieldPoints);
            bspliner->Update();
            m_PhiLattice = bspliner->GetPhiLattice();
            m_DisplacementField = bspliner->GetOutput();
            cout << "BSpline Update Time: " << timer.elapsed() << endl;
        } catch (itk::ExceptionObject& e) {
            e.Print(std::cout);
        }
    }

    SliceType::Pointer BSplineRegistration::WarpImage(SliceType::Pointer srcImage) {
        if (m_DisplacementField.IsNull()) {
            return SliceType::Pointer(NULL);
        }
        WarpImageFilterType::Pointer warpFilter = WarpImageFilterType::New();
        warpFilter->SetInput(srcImage);
        warpFilter->SetDisplacementField(m_DisplacementField);
//        warpFilter->SetOutputSpacing(srcImage->GetSpacing());
//        warpFilter->SetOutputOrigin(srcImage->GetOrigin());
//        warpFilter->SetOutputSize(srcImage->GetBufferedRegion().GetSize());
        warpFilter->SetOutputParametersFromImage(srcImage);
        warpFilter->Update();
        return warpFilter->GetOutput();
    }


    // compute determinant of jacobian image
    SliceType::Pointer BSplineRegistration::GetDeterminantOfJacobian() {
        itkcmds::itkImageIO<SliceType> io;
        SliceType::Pointer detImage = io.NewImageT(m_RefImage);

        FieldTransformType::Pointer txf = GetTransform();
        SliceIteratorType iter(detImage, detImage->GetBufferedRegion());
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            FieldTransformType::JacobianType jacob;
            txf->ComputeJacobianWithRespectToParameters(iter.GetIndex(), jacob);
            double det = vnl_determinant(jacob);
            iter.Set(det);
        }
        return detImage;
    }


    DisplacementFieldType::Pointer BSplineRegistration::GetDisplacementField() {
        return m_DisplacementField;
    }

    DisplacementFieldType::Pointer BSplineRegistration::GetControlPoints() {
        return m_PhiLattice;
    }

    SliceType::Pointer BSplineRegistration::GetDisplacementMagnitude() {
        if (m_DisplacementField.IsNull()) {
            return SliceType::Pointer(NULL);
        }

        typedef itk::VectorMagnitudeImageFilter<DisplacementFieldType, SliceType> VectorFilterType;
        VectorFilterType::Pointer filter = VectorFilterType::New();
        filter->SetInput(m_DisplacementField);
        filter->Update();
        return filter->GetOutput();
    }
    
    FieldTransformType::Pointer BSplineRegistration::GetTransform() {
        FieldTransformType::Pointer txf = FieldTransformType::New();
        txf->SetDisplacementField(m_DisplacementField);
        return txf;
    }

    SliceTransformType::Pointer BSplineRegistration::GetRawTransform() {
        return SliceTransformType::Pointer(NULL);
    }

}