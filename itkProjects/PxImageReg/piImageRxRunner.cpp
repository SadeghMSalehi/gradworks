//
//  piImageRxRunner.cpp
//  PxImageReg
//
//  Created by Joohwi Lee on 11/30/13.
//
//

#include "piImageRxRunner.h"
#include "piOptions.h"
#include "piImageIO.h"
#include "piImageDef.h"
#include "piImageProc.h"
#include "piEntropyComputer.h"

#include <itkObjectToObjectOptimizerBase.h>
#include <itkCompositeTransform.h>
#include <itkEuler2DTransform.h>
#include <itkEuler3DTransform.h>
#include <itkMeanSquaresImageToImageMetricv4.h>
#include <itkCorrelationImageToImageMetricv4.h>
#include <itkMattesMutualInformationImageToImageMetricv4.h>
#include <itkRegistrationParameterScalesFromPhysicalShift.h>
#include <itkGradientDescentOptimizerv4.h>
#include <itkWindowConvergenceMonitoringFunction.h>
#include "itk/itkEntropyImageToImageMetricv4.h"

namespace pi {
    typedef itk::ObjectToObjectOptimizerBase::ScalesType ScalesType;
    typedef itk::ObjectToObjectOptimizerBase::ParametersType ParametersType;
    typedef itk::ObjectToObjectMetricBase::DerivativeType DerivativeType;

    void executeRxRunner(Options& opts, StringVector& args) {
        ImageRx runner(opts, args);
        runner.mainRigidRegistration(opts, args);
        exit(0);
    }

#pragma mark TransformParameterScalesEstimator Implementations
    class TransformParameterScalesEstimator {
    public:
        TransformParameterScalesEstimator() {
            _smallParameterVariation = 0.01;
        }

        void estimateScales(ScalesType& scales);
        float estimateStepScales(const DerivativeType& step);
        void estimateLocalStepScales(const ParametersType& step, ParametersType& localStepScales);
        float estimateMaximumStepSize();

        inline void setSmallParameterVariation(double var) { _smallParameterVariation = var; }
        inline double getSmallParameterVariation() { return _smallParameterVariation; }
        inline void setTransform(TransformType* transform) { _transform = transform; }
        inline void setImage(RealImage::Pointer image) { _baseImage = image; }


    protected:
        int getDimension();

        bool isDisplacementFieldTransform();
        bool isBsplineTransform();
        bool transformHasLocalSupportForScalesEstimation();

        void computeSquaredJacobianNorms(const RealImage::PointType& p, ParametersType& squareNorms);


        // currently sample with corners
        // there can be more sampling strategies than sampling with corners
        void sampleDomain();
        void sampleDomainWithCorners();

        // compute the shift in voxels when deltaParameters is applied onto the current parameters
        float computeMaximumVoxelShift(const ParametersType& delta);
        void computeSampleShifts(const ParametersType& delta, ScalesType& localShifts);

    private:
        RealImage::Pointer _baseImage;
        TransformType* _transform;
        std::vector<RealImage::PointType> _samplePoints;
        double _smallParameterVariation;
    };


    // get the image dimension
    int TransformParameterScalesEstimator::getDimension() {
        return RealImage::ImageDimension;
    }

    // test if the transform is a displacement field transform
    bool TransformParameterScalesEstimator::isDisplacementFieldTransform() {
        return _transform->GetTransformCategory() == TransformType::DisplacementField;
    }

    // test if the transform is a b-spline transform
    bool TransformParameterScalesEstimator::isBsplineTransform() {
        bool isBsplineTransform = false;
        if (_transform->GetTransformCategory() == TransformType::BSpline) {
            isBsplineTransform = true;
        }
        if (!isBsplineTransform) {
            typedef itk::CompositeTransform<float, RealImage::ImageDimension> CompositeTransformType;
            CompositeTransformType::Pointer compositeTransform = dynamic_cast<CompositeTransformType *>(_transform);
            if (compositeTransform) {
                isBsplineTransform = true;

                // loop over composite transforms and if they are all B-spline, it is Bspline transform
                const int nTransforms = static_cast<int>(compositeTransform->GetNumberOfTransforms());
                for (int tind = nTransforms - 1; tind >= 0; tind--) {
                    if (compositeTransform->GetNthTransformToOptimize(tind) &&
                       (compositeTransform->GetNthTransform(tind).GetPointer()->GetTransformCategory() != TransformType::BSpline))
                    {
                        isBsplineTransform = false;
                        break;
                    }
                }
            }
        }
        return isBsplineTransform;
    }

    // test if the transform has local support
    bool TransformParameterScalesEstimator::transformHasLocalSupportForScalesEstimation() {
        return isDisplacementFieldTransform() || isBsplineTransform();
    }

    void TransformParameterScalesEstimator::computeSquaredJacobianNorms(const RealImage::PointType &p, ParametersType &squareNorms) {
        TransformType::JacobianType jacobian;
        const int numPara = _transform->GetNumberOfLocalParameters();
        const int dim = getDimension();

        _transform->ComputeJacobianWithRespectToParameters(p, jacobian);
        for (int p = 0; p < numPara; p++) {
            squareNorms[p] = 0;
            for (int d = 0; d < dim; d++) {
                squareNorms[p] += jacobian[d][p] * jacobian[d][p];
            }
        }
    }

    // just do sample domain with corners
    void TransformParameterScalesEstimator::sampleDomain() {
        sampleDomainWithCorners();
    }

    // sample corners of the image
    void TransformParameterScalesEstimator::sampleDomainWithCorners() {
        RealImage::RegionType region = _baseImage->GetBufferedRegion();
        RealImage::IndexType firstCorner = region.GetIndex();
        RealImage::IndexType corner;
        RealImage::PointType point;

        RealImage::SizeType size = region.GetSize();
        const int cornerNumber = 1 << getDimension();

        _samplePoints.resize(cornerNumber);

        // radix counting
        for (int i = 0; i < cornerNumber; i++) {
            int bit;
            for (int d = 0; d < getDimension(); d++) {
                bit = (int) ((i & (1 << d)) != 0);
                corner[d] = firstCorner[d] + bit * (size[d] - 1);
            }
            _baseImage->TransformIndexToPhysicalPoint(corner, point);
            _samplePoints[i] = point;
        }
    }

    // make small perturbation in transform parameters
    // and measure its displacement
    void TransformParameterScalesEstimator::computeSampleShifts(const ParametersType &delta, ScalesType &sampleShifts) {
        const ParametersType oldParams = _transform->GetParameters();
        const int numSamples = _samplePoints.size();

        RealImage::PointType point;
        TransformType::OutputPointType newMappedVoxel;

        // temporary storage to store original points
        std::vector<TransformType::OutputPointType> oldMappedVoxels(numSamples);
        sampleShifts.SetSize(numSamples);

        // store previous sample points
        for (int c = 0; c < numSamples; c++) {
            RealImage::PointType point = _samplePoints[c];
            oldMappedVoxels[c] = _transform->TransformPoint(point);
        }

        // make perturbation
        _transform->UpdateTransformParameters(delta);

        // measure its distance
        for (int c = 0; c < numSamples; c++) {
            point = _samplePoints[c];
            newMappedVoxel = _transform->TransformPoint(point);
            sampleShifts[c] = newMappedVoxel.EuclideanDistanceTo(oldMappedVoxels[c]);
        }

        // restore the parameters in the transform
        _transform->SetParameters(oldParams);
    }

    // Estimate scales
    void TransformParameterScalesEstimator::estimateScales(ScalesType& paramScales) {
        sampleDomain();

        const int nAllParams = _transform->GetNumberOfParameters();
        const int nLocalParams = _transform->GetNumberOfLocalParameters();
        paramScales.SetSize(nLocalParams);

        float maxShift = 0;

        const float fMax = itk::NumericTraits<float>::max();
        const float fEpsilion = itk::NumericTraits<float>::epsilon();
        float minNonZeroShift = fMax;

        ParametersType deltaParameters(nAllParams);
        itk::OffsetValueType offset = 0;

        // offset should be different for a displacementfield transform

        for (int i = 0; i < nLocalParams; i++) {
            // refill deltaParameters with zeros at each loop
            // so that a variation per parameter can be measured
            deltaParameters.Fill(0);
            deltaParameters[offset + i] = _smallParameterVariation;
            maxShift = computeMaximumVoxelShift(deltaParameters);
            paramScales[i] = maxShift;
            if (maxShift > fEpsilion && maxShift < minNonZeroShift) {
                minNonZeroShift = maxShift;
            }
        }

        // if there is no variation and minNonZeroShift wasn't changed
        if (minNonZeroShift == fMax) {
            cout  << "Variation in any parameter won't change a voxel position. The default scales (1.0) are used to avoid division-by-zero." << endl;
            paramScales.Fill(1);
        } else {
            for (int i = 0; i < nLocalParams; i++) {
                if (paramScales[i] <= fEpsilion) {
                    // if scale is too small, then set it same as the minimum shift value
                    paramScales[i] = minNonZeroShift * minNonZeroShift;
                } else {
                    // otherwise, make it square
                    paramScales[i] *= paramScales[i];
                }
                // divide by the square of the perturbation delta
                paramScales[i] *= 1 / vnl_math_sqr(_smallParameterVariation);
            }
        }
    }

    /** 
     * Compute the scale for a step. For transform T(x + t * step), the scale
     * w.r.t. the step is the shift produced by step.
     */
    float TransformParameterScalesEstimator::estimateStepScales(const DerivativeType &step) {
        sampleDomain();
        if (transformHasLocalSupportForScalesEstimation()) {
            return computeMaximumVoxelShift(step);
        }

        const float fEpsilon = itk::NumericTraits<float>::epsilon();

        // For global transforms, we want a linear approximation of the function
        // of step scale w.r.t "step". This is true only when "step" is close to
        // zero. Therefore, we need to scale "step" down.
        float maxGradient = 0;

        // identify the maximum gradient step
        for (int i = 0; i < step.GetSize(); i++) {
            if (maxGradient < vcl_abs(step[i])) {
                maxGradient = vcl_abs(step[i]);
            }
        }

        // if it is too small then return 0
        if (maxGradient <= fEpsilon) {
            return 0;
        }

        // scale down the gradient so that
        // the maximum gradient step equals to the perturbation delta
        float factor = _smallParameterVariation / maxGradient;
        ParametersType smallStep(step.size());
        smallStep = step * factor;

        // restore the original step size
        return computeMaximumVoxelShift(smallStep) / factor;
    }

    // This shouldn't be used until B-spline or DisplacementField transform is used
    void TransformParameterScalesEstimator::estimateLocalStepScales(const ParametersType &step, ParametersType &localStepScales) {
        if (!transformHasLocalSupportForScalesEstimation()) {
            throw "EstimateLocalStepScales: the transform doesn't have local support (displacement field or b-spline).";
        }

        sampleDomain();

        // collect sample shifts
        ScalesType sampleShifts;
        computeSampleShifts(step, sampleShifts);

        const int nAllParams = _transform->GetNumberOfParameters();
        const int nParams = _transform->GetNumberOfLocalParameters();
        const int nLocals = nAllParams / nParams;

        localStepScales.SetSize(nLocals);
        localStepScales.Fill(0);

        const int nSamples = _samplePoints.size();
        for (int c = 0; c < nSamples; c++) {
            RealImage::PointType &point = _samplePoints[c];
            int localId = 0; // computeParameterOffsetFromPoint
            localStepScales[localId] = sampleShifts[c];
        }
    }


    float TransformParameterScalesEstimator::estimateMaximumStepSize() {
        const RealImage::SpacingType& spacing = _baseImage->GetSpacing();
        const int dim = getDimension();

        float minSpacing = itk::NumericTraits<float>::max();

        for (int d = 0; d < dim; d++) {
            if (minSpacing > spacing[d]) {
                minSpacing = spacing[d];
            }
        }
        return minSpacing;
    }

    // Compute the maximum shift when a transform is changed with deltaParameters
    float TransformParameterScalesEstimator::computeMaximumVoxelShift(const ParametersType &delta) {
        ScalesType sampleShifts;

        // make little perturbation on the transform
        // and compute sample points' shift
        computeSampleShifts(delta, sampleShifts);

        float maxShift = itk::NumericTraits<float>::Zero;
        for (unsigned int s = 0; s < sampleShifts.size(); s++) {
            // if sampleShifts[s] is larger than maxShift, swap maxShift
            if (maxShift < sampleShifts[s]) {
                maxShift = sampleShifts[s];
            }
        }
        return maxShift;
    }



#pragma mark TransformingImage Implementations
    TransformingImage::TransformingImage()
    {
        _transform = NULL;
    }

    // initialize TransformingImage with a given intensity image
    // compute the gadient image and set up interpolators
    void TransformingImage::setImage(RealImage::Pointer image) {
        _image = image;

        // use the half size of spacing as a sigma
        double sigma = _image->GetSpacing()[0] / 2.0;

        // compute the gradient image using ImageProc
        _gradientImage = ComputeGaussianGradient(_image, sigma);

        // set up interpolators
        _intensityInterpolator = itk::LinearInterpolateImageFunction<RealImage>::New();
        _intensityInterpolator->SetInputImage(_image);

        _gradientInterpolator = itk::VectorLinearInterpolateImageFunction<GradientImage>::New();
        _gradientInterpolator->SetInputImage(_gradientImage);
    }

    RealImage::PixelType TransformingImage::samplePixel(RealImage::PointType fixedPoint, bool& isValidPixel) {
        // first compute the transformed point at this moving image
        if (_transform) {
            RealImage::PointType movingPoint = _transform->TransformPoint(fixedPoint);
            if (_intensityInterpolator->IsInsideBuffer(movingPoint)) {
                isValidPixel = true;
                return _intensityInterpolator->Evaluate(movingPoint);
            } else {
                // if the moving point is outside of the image buffer
                // the return pixel is not defined
                isValidPixel = false;
                return 0;
            }
        } else {
            return _intensityInterpolator->Evaluate(fixedPoint);
        }
    }


    GradientImage::PixelType TransformingImage::sampleGradient(RealImage::PointType fixedPoint, bool& isValidPixel) {
        // first compute the transformed point at this moving image
        if (_transform) {
            RealImage::PointType movingPoint = _transform->TransformPoint(fixedPoint);
            if (_gradientInterpolator->IsInsideBuffer(movingPoint)) {
                isValidPixel = true;
                return _gradientInterpolator->Evaluate(movingPoint);
            } else {
                // if the moving point is outside of the image buffer
                // the return gradient is not defined
                isValidPixel = false;
                return GradientImage::PixelType();
            }
        } else {
            return _gradientInterpolator->Evaluate(fixedPoint);
        }
    }

    // set a transform instance
    // and assume that its fixed parameters are the center of the image
    void TransformingImage::setTransform(TransformType *transform) {
        if (_image.IsNull()) {
            cout << "TransformingImage::setTransform() must have its image first." << endl;
            exit(0);
        }

        _transform = transform;

        // compute the image center to set the fixed parameters
        RealIndex centerIdx = _image->GetBufferedRegion().GetIndex();
        RealImage::PointType centerPoint;
        for (int i = 0; i < RealImage::ImageDimension; i++) {
            centerIdx[i] += _image->GetBufferedRegion().GetSize(i) / 2.0;
        }
        _image->TransformContinuousIndexToPhysicalPoint(centerIdx, centerPoint);


        // fixed parameters
        ParametersType fixedParams;
        fixedParams.SetSize(RealImage::ImageDimension);
        for (int i = 0; i < RealImage::ImageDimension; i++) {
            fixedParams[i] = centerPoint[i];
        }
        _transform->SetFixedParameters(fixedParams);
    }

    // provide parameters to the transform instance
    void TransformingImage::setTransformParameters(TransformType::ParametersType params) {
        if (_transform) {
            _transform->SetParameters(params);
        }
    }

    // resample image with the transform
    RealImage::Pointer TransformingImage::getResampledImage(RealImage::Pointer fixedImage) {
        typedef itk::ResampleImageFilter<RealImage,RealImage> ResampleFilter;
        ResampleFilter::Pointer resampler = ResampleFilter::New();
        resampler->SetInput(_image);
        if (_transform) {
            resampler->SetTransform(dynamic_cast<ResampleFilter::TransformType*>(_transform));
        }
        resampler->SetReferenceImage(fixedImage);
        resampler->UseReferenceImageOn();
        resampler->Update();
        return resampler->GetOutput();
    }

#pragma mark EntropyCostFunction Implementations
    typedef std::vector<ParametersType> ParametersList;
    typedef std::vector<ScalesType> ScalesList;
    typedef std::vector<DerivativeType> DerivativeList;
    typedef std::vector<float> FloatVector;

    // cost function for entropy computation
    class EntropyCostFunction {
    public:
        void estimateScalesAndStepSize(ScalesType& scales, float& maxStep);
        void computeValueAndDerivatives(double& value, DerivativeList& derivs);
        void updateTransformParameters(DerivativeList& gradient, int iter);
        void rescaleGradient(DerivativeList& gradient);
        void setParametersList(ParametersList& paramsList);
        void getParametersList(ParametersList& paramsList);

        // inline functions
        inline TransformingImage::Vector& getMovingImages() { return _movingImages; }
        inline TransformingImage& getFixedImage() { return _fixedImage; }

        // temporary metric setup for validation
        void setupMetric();

    private:
        ScalesType _scales;
        TransformParameterScalesEstimator _estimator;
        TransformingImage _fixedImage;
        TransformingImage::Vector _movingImages;
        float _maximumStepSizeInPhysicalUnit;
        FloatVector _learningRate;
        ParametersList _currentParamsList;

        typedef itk::MeanSquaresImageToImageMetricv4<RealImage, RealImage> MetricType;
        MetricType::Pointer _metric;
    };

    // set up initial parameters
    void EntropyCostFunction::setParametersList(ParametersList& paramsList) {
        _currentParamsList = paramsList;
        for (int i = 0; i < _movingImages.size(); i++) {
            _movingImages[i].getTransform()->SetParameters(paramsList[i]);
        }
    }

    // get current parameters
    void EntropyCostFunction::getParametersList(ParametersList& paramsList) {
        _currentParamsList = paramsList;
        for (int i = 0; i < _movingImages.size(); i++) {
            paramsList[i] = _movingImages[i].getTransform()->GetParameters();
        }
    }

    // estimator for parameter scales and maximumStepSize
    void EntropyCostFunction::estimateScalesAndStepSize(ScalesType& scales, float& maxStep) {
        maxStep = _estimator.estimateMaximumStepSize();
        _estimator.estimateScales(scales);
    }

    // update transform parameters
    void EntropyCostFunction::updateTransformParameters(DerivativeList &gradient, int iter) {
        const int nImages = _movingImages.size();
        _estimator.estimateScales(_scales);
        _maximumStepSizeInPhysicalUnit = _estimator.estimateMaximumStepSize();

        // modify gradient so that the maximum displacement would not exceed the maximum step size
        rescaleGradient(gradient);
        for (int i = 0; i < nImages; i++) {
            _movingImages[i].getTransform()->Print(cout);
            _movingImages[i].getTransform()->UpdateTransformParameters(gradient[i]);
        }
    }

    // modify gradient so that the maximum displacement is limited
    // below the expected minimum displacement
    void EntropyCostFunction::rescaleGradient(DerivativeList &gradient) {
        const float fEpsilon = itk::NumericTraits<float>::epsilon();
        const int nImages = gradient.size();
        _learningRate.resize(nImages);
        // loop over the gradients
        for (int i = 0; i < nImages; i++) {
            // normalize unit difference between parameters
            for (int j = 0; j < gradient[i].size(); j++) {
                gradient[i][j] /= _scales[j];
            }

            // restrict the maximum displacement in according to the gradient
            // 1) normalize the gradient so that the maximum element is equal to _smallPerturbation
            // 2) compute the displacement by the maximum element and divide by _smallPerturbation
            // 3) move to the direction but with the small step size
            float stepScale = _estimator.estimateStepScales(gradient[i]);
            if (stepScale <= fEpsilon) {
                stepScale = 1;
            } else {
                _learningRate[i] = _maximumStepSizeInPhysicalUnit / stepScale;
            }
            gradient[i] *= _learningRate[i];
        }
    }

    // compute value and derivatives of the cost function
    void EntropyCostFunction::computeValueAndDerivatives(double& value, DerivativeList &derivs) {

//        _metric->GetValueAndDerivative(value, derivs[0]);

        const int nPixels = _fixedImage.getImage()->GetPixelContainer()->Size();
        EntropyComputer<> comp(_movingImages.size() + 1, nPixels, 1);

        EntropyComputer<>::Iterator iter = comp.dataIter;
        iter.FirstData();
        for (int j = 0; j < nPixels; j++) {
            iter.data[j] = _fixedImage.getImage()->GetBufferPointer()[j];
        }
        iter.NextData();

        // loop over the moving image
        for (int i = 0; i < _movingImages.size(); i++) {
            // resample the image with respect to current transform parameters
            RealImage::Pointer resampledImage = _movingImages[i].getResampledImage(_fixedImage.getImage());

            // store the resampled pixel into the entropy computer
            RealImage::PixelType* pix = resampledImage->GetBufferPointer();
            for (int j = 0; j < nPixels; j++) {
                iter.data[j] = pix[j];
            }

            // go to next data
            iter.NextData();
        }

        comp.MoveToCenter();
        comp.ComputeCovariance();
        comp.ComputeGradient();

        value = comp.ComputeEntropy();
        cout << comp.covariance << endl;

        // compute the gradient with respect to parameters
        const int nParams = _movingImages[0].getTransform()->GetNumberOfParameters();

        // allocate the derivatives list
        derivs.resize(_movingImages.size());

        for (int i = 0; i < derivs.size(); i++) {
            // for each image, the derivative has nParams elements
            derivs[i].SetSize(nParams);

            // resample gradient from the resampled moving image
            // i don't think this gonna work, but give it a try
            GradientInterpolatorType::Pointer gradIntp = GradientInterpolatorType::New();
            RealImage::Pointer resampledMovingImage = _movingImages[0].getResampledImage(_fixedImage.getImage());

            GradientImage::Pointer gradientResampledImage = ComputeGaussianGradient(resampledMovingImage, resampledMovingImage->GetSpacing()[0]/2.0);
            gradIntp->SetInputImage(gradientResampledImage);

            // for each derivative with respect to a parameter
            // j is the index for parameter
            for (int j = 0; j < derivs[i].size(); j++) {
                // initialize to zero
                derivs[i][j] = 0;
                // temporary sum
                double dSum = 0;

                // iterator to compute the jacobian
                RealImageIteratorType iter(_fixedImage.getImage(), _fixedImage.getImage()->GetBufferedRegion());
                iter.GoToBegin();

                // iterate over each pixel
                for (int k = 0; k < nPixels; k++) {
                    // compute the physical point
                    // maybe cache saves the computation time
                    RealImage::PointType kPoint;
                    _fixedImage.getImage()->TransformIndexToPhysicalPoint(iter.GetIndex(), kPoint);

                    // compute the jacobian at kPoint
                    TransformType::JacobianType jacobian;
                    _movingImages[i].getTransform()->ComputeJacobianWithRespectToParameters(kPoint, jacobian);


                    // compute the total derivative with respect to a parameter
                    // this is done by the inner product of the gradient of image
                    // and the l-th column of the jacobian matrix
                    const int dim = RealImage::ImageDimension;

                    // l is the index representing the x-y-z dimension
                    bool isValidPixel;

                    // compute the moving image's point to sample its gradient
                    RealImage::PointType movingPoint = _movingImages[0].getTransform()->TransformPoint(kPoint);
                    
                    GradientImage::PixelType imageGradient = _movingImages[0].sampleGradient(kPoint, isValidPixel);

                    double dI = 0;
                    for (int l = 0; l < dim; l++) {
                        if (isValidPixel) {
                            dI += (imageGradient[l] * jacobian[l][j]);
                        }
                    }
                    // compute the total derivative dI_k
                    dI *= comp.gradient[i+1][k];

                    // by chain-rule
                    dSum += dI;
                    ++iter;
                } // k-loop

                derivs[i][j] = dSum;
            } // j-loop
        }

    }

    // metric set up
    void EntropyCostFunction::setupMetric() {
        _metric = MetricType::New();
        _metric->SetFixedImage(_fixedImage.getImage());
        _metric->SetMovingImage(_movingImages[0].getImage());
        _metric->SetFixedInterpolator(_fixedImage.getInterpolator());
        _metric->SetMovingInterpolator(_movingImages[0].getInterpolator());
        _metric->SetMovingTransform(_movingImages[0].getTransform());
        _metric->Initialize();

        _estimator.setImage(_movingImages[0].getImage());
        _estimator.setTransform(_movingImages[0].getTransform());
    }


#pragma mark GradientDescentOptimizer
    // Gradient Descent Optimizer
    // This class updates the cost function at each iteration
    // The scaling of gradient value should be done by the cost function
    // The parameters and gradient values are also maintained by the cost function
    // This dependency pattern is most different from the original itk implementation.
    class GradientDescentOptimizer {
    public:

        GradientDescentOptimizer() {
            _maxIterations = 500;
        }

        void startOptimization();
        void stopOptimization();
        void resumeOptimization();
        void advanceOneStep();

        inline void setCostFunction(EntropyCostFunction costFunc) {
            _costFunc = costFunc;
        }

    private:
        itk::Function::WindowConvergenceMonitoringFunction<>::Pointer _convergenceFunction;
        EntropyCostFunction _costFunc;
        ParametersList _paramsList;
        DerivativeList _gradientList;
        double _currentCost;
        int _currentIteration, _maxIterations;
        bool _stop;
    };

    // initialize variables to start optimization
    void GradientDescentOptimizer::startOptimization() {
        // temporarily set to 1
        _gradientList.resize(1);

        _stop = false;
        _currentIteration = 0;

        // initialize the convergence function
        _convergenceFunction = itk::Function::WindowConvergenceMonitoringFunction<>::New();
        _convergenceFunction->SetWindowSize(50);

        while (!_stop) {
            resumeOptimization();
        }
    }

    // continue optimization by computing gradient and advance one step
    void GradientDescentOptimizer::resumeOptimization() {
        if (!_stop) {

            // begin with iteration #1
            _currentIteration++;

            _costFunc.computeValueAndDerivatives(_currentCost, _gradientList);
            _convergenceFunction->AddEnergyValue(_currentCost);

            if (_convergenceFunction->GetConvergenceValue() < 1e-8) {
                stopOptimization();
                return;
            }

            cout << "Iter: " << _currentIteration << " Cost: " << _currentCost << endl;
            advanceOneStep();

            if (_currentIteration >= _maxIterations) {
                stopOptimization();
            }
        }
    }

    // modify gradient and update transform parameters to compute the next value
    void GradientDescentOptimizer::advanceOneStep() {
        _costFunc.updateTransformParameters(_gradientList, _currentIteration);
    }

    // stop optimization
    void GradientDescentOptimizer::stopOptimization() {
        _stop = true;
    }
    



#pragma mark EntropyImageMetric Implementations
void EntropyImageMetric::setFixedImage(RealImage::Pointer image) {
    // the fixedImage will not have transform
    _fixedImage.setImage(image);
}

void EntropyImageMetric::addMovingImage(RealImage::Pointer image, TransformType* transform, TransformType::ParametersType& params) {
    // the movingImage will have its own transform and its initial parameters
    TransformingImage movingImage;
    movingImage.setImage(image);
    movingImage.setTransform(transform);
    movingImage.setTransformParameters(params);

    // add to moving image vectors so that additional images can be added
    _movingImages.push_back(movingImage);
}

double EntropyImageMetric::GetValue() const {
    // get value will return \log \det (Covariance)
    return 0;
}

void EntropyImageMetric::GetValueAndDerivative(double &value, DerivativeType &deriv) const {
    // compute value and its derivatives
    return;
}


#pragma mark OptimizerProgress Listener

    /**
     * This class is called from a registration module in order to provide internal registration information to a user
     *
     */
    template <class T>
    class OptimizerProgress: public itk::Command {
    public:
        /** Standard class typedefs. */
        typedef OptimizerProgress Self;
        typedef itk::Command Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        itkTypeMacro(OptimizerProgress, itk::Command);
        itkNewMacro(OptimizerProgress);

        virtual void Execute(Object *caller, const itk::EventObject & event) {
            T* realCaller = dynamic_cast<T*>(caller);
            if (realCaller == NULL) {
                return;
            }
            cout << realCaller->GetCurrentIteration() << ": " << realCaller->GetCurrentPosition() << endl;
            realCaller->Print(cout);
        }

        /**
         * when a caller is defined as const
         */
        virtual void Execute(const Object *caller, const itk::EventObject & event) {
        }

    protected:
        OptimizerProgress() {}
        virtual ~OptimizerProgress() {}

    private:
        // purposely not implemented
        OptimizerProgress(const Self &);
        void operator=(const Self &);
    };




#pragma mark Registration Functions

    /***
     * perform rigid registration
     *
     * command line:
     * --rx fixedImage movingImage [resampledImage] [transformParameters]
     *
     */

    void ImageRx::mainRigidRegistration(Options& opts, StringVector& args) {

        if (args.size() < 4) {
            cout << "--rx [fixed-image] [moving-image] [output-image] [output-transform]" << endl;
            return;
        }

#define USE_AFFINE_TRANSFORM

#ifdef USE_RIGID_TRANSFORM
#if DIMENSIONS == 2
        typedef itk::Euler2DTransform<> TransformType;
#else
        typedef itk::Euler3DTransform<> TransformType;
#endif
#endif

#ifdef USE_AFFINE_TRANSFORM
        typedef itk::AffineTransform<double,DIMENSIONS> TransformType;
#endif

        ImageIO<RealImage> imageIO;

        RealImage::Pointer fixedImage  = imageIO.ReadCastedImage(args[0]);
        RealImage::Pointer movingImage = imageIO.ReadCastedImage(args[1]);


        // cost function setup
        
        EntropyCostFunction myCostFunc;
        myCostFunc.getFixedImage().setImage(fixedImage);

        TransformType::Pointer movingTransform = TransformType::New();
        TransformingImage::Vector& movingImages = myCostFunc.getMovingImages();
        movingImages.resize(1);
        movingImages[0].setImage(movingImage);
        movingImages[0].setTransform(movingTransform);

        // initialize metric
        myCostFunc.setupMetric();


        // initial parameters
        // initial transform parameters
        ParametersType myInitialParams;
        myInitialParams.SetSize(movingImages[0].getTransform()->GetNumberOfParameters());
        myInitialParams.Fill(0);
        myInitialParams[0] = 1;
        myInitialParams[3] = 1;

        ParametersList paramList;
        paramList.push_back(myInitialParams);
        myCostFunc.setParametersList(paramList);

        GradientDescentOptimizer opti;
        opti.setCostFunction(myCostFunc);
        opti.startOptimization();

        RealImage::Pointer resampledImage = movingImages[0].getResampledImage(fixedImage);
        imageIO.WriteImage(args[2], resampledImage);
        imageIO.WriteSingleTransform(args[3].c_str(), movingImages[0].getTransform());

        exit(0);

        // There are three main components in the ITK registration framework.

        // 1) Optimizer computes a gradient descent direction and updates the parameters toward a minimizing direction.
        // 2) Similarity Metric provides a cost function as well as derivative functions.
        // 3) Transformation and its parameters transform the fixed image space to moving image so that a resampled image can be computed
        //

        //        typedef itk::TranslationTransform<double,3> TransformType;
        //        typedef itk::VersorRigid3DTransform<double> TransformType;
        //        typedef itk::Similarity3DTransform<double> TransformType;
        //        typedef itk::BSplineTransform<double,2,4> TransformType;
        //        typedef itk::AffineTransform<double,3> TransformType;
        TransformType::Pointer transform;
        transform = TransformType::New();

        // define a gradient descent optimizer function
        typedef itk::GradientDescentOptimizerv4 OptimizerType;
        OptimizerType::Pointer optimizer;

        // define a MeanSquares image to image metric
        typedef OptimizerType::ParametersType ParametersType;





        // set up initial parameters
        ParametersType initialParams;
        initialParams.SetSize(transform->GetNumberOfParameters());
        initialParams.Fill(0);

#ifdef USE_AFFINE_TRANSFORM
        initialParams[0] = 1;
        initialParams[3] = 1;
#endif


        // set up fixed parameters: the center of rotation
        ParametersType fixedParams;
        fixedParams.SetSize(fixedImage->GetImageDimension());

        itk::ContinuousIndex<double,RealImage::ImageDimension> centerIdx;
        RealImage::PointType centerPoint;

        for (int i = 0; i < fixedParams.GetSize(); i++) {
            centerIdx[i] = fixedImage->GetBufferedRegion().GetIndex(i) + fixedImage->GetBufferedRegion().GetSize(i) / 2.0;
        }
        fixedImage->TransformContinuousIndexToPhysicalPoint(centerIdx, centerPoint);
        for (int i = 0; i < fixedParams.size(); i++) {
            fixedParams[i] = centerPoint[i];
        }
        transform->SetFixedParameters(fixedParams);


        typedef itk::LinearInterpolateImageFunction<RealImage> ImageInterpolator;

        // set up cost functions
        // the parameters are different for each cost function

#define USE_MSQ

#ifdef USE_MSQ
        typedef itk::MeanSquaresImageToImageMetricv4<RealImage, RealImage> CostFunctionType;
#endif
#ifdef USE_CC
        typedef itk::CorrelationImageToImageMetricv4<RealImage, RealImage> CostFunctionType;
#endif


        CostFunctionType::Pointer costFunc;
        costFunc = CostFunctionType::New();
        costFunc->SetFixedImage(fixedImage);
        costFunc->SetFixedInterpolator(ImageInterpolator::New());
        costFunc->SetMovingImage(movingImage);
        costFunc->SetMovingInterpolator(ImageInterpolator::New());
        costFunc->SetMovingTransform(transform);
        costFunc->SetParameters(initialParams);
        costFunc->Initialize();


#ifdef USE_MI
        typedef itk::MattesMutualInformationImageToImageMetricv4<RealImage, RealImage> CostFunctionType;
        CostFunctionType::Pointer costFunc;
        costFunc = CostFunctionType::New();
        costFunc->SetFixedImage(fixedImage);
        costFunc->SetFixedInterpolator(ImageInterpolator::New());
        costFunc->SetMovingImage(movingImage);
        costFunc->SetMovingInterpolator(ImageInterpolator::New());
        costFunc->SetMovingTransform(transform);
        costFunc->SetParameters(initialParams);
        if (dynamic_cast<itk::MattesMutualInformationImageToImageMetricv4<RealImage, RealImage>*>(costFunc.GetPointer()) != NULL) {
            costFunc->SetNumberOfHistogramBins(32);
        }
        costFunc->Initialize();
#endif

        typedef itk::RegistrationParameterScalesFromPhysicalShift<CostFunctionType> ScaleEstimatorType;
        ScaleEstimatorType::Pointer estimator = ScaleEstimatorType::New();
        estimator->SetMetric(costFunc);

        ScalesType affineScales;
        estimator->EstimateScales(affineScales);

        cout << affineScales << endl;

        affineScales.Fill(0);

        TransformParameterScalesEstimator scalesEstimator;
        scalesEstimator.setImage(fixedImage);
        scalesEstimator.setTransform(transform);
        scalesEstimator.estimateScales(affineScales);
        cout << affineScales << endl;

        exit(0);

        // print out optimizer progress to console
        OptimizerProgress<OptimizerType>::Pointer progress = OptimizerProgress<OptimizerType>::New();

        // normalize paramter update scaling
        OptimizerType::ScalesType scales;
        scales.SetSize(transform->GetNumberOfParameters());
        scales.Fill(1);


        // set up the optimizer
        optimizer = OptimizerType::New();
        optimizer->SetScales(scales);
        optimizer->SetScalesEstimator(estimator);
        optimizer->SetMetric(costFunc);

        // number of iterations
        optimizer->SetNumberOfIterations(1000);
        RealImage::SpacingType spacing = fixedImage->GetSpacing();
        optimizer->SetMaximumStepSizeInPhysicalUnits(spacing[0]*3);

        // set the optimizer observer
        optimizer->AddObserver(itk::IterationEvent(), progress);

        try {
            ::itk::Object::GlobalWarningDisplayOn();
            optimizer->SetDebug(true);
            costFunc->SetDebug(true);
            optimizer->Print(cout);
            optimizer->StartOptimization();
        } catch (itk::ExceptionObject& ex) {
            ex.Print(cout);
        }

        cout << "Iterations: " << optimizer->GetCurrentIteration() << endl;
        cout << "Stop Reason: " << optimizer->GetStopConditionDescription() << endl;

        const ParametersType& params = optimizer->GetCurrentPosition();
        transform->SetParameters(params);

        typedef itk::ResampleImageFilter<RealImage,RealImage> ResampleFilter;
        ResampleFilter::Pointer resampler = ResampleFilter::New();
        resampler->SetInput(movingImage);
        if (transform.IsNotNull()) {
            cout << "Setting transform:" << endl;
            transform->Print(cout);
            resampler->SetTransform(dynamic_cast<ResampleFilter::TransformType*>(transform.GetPointer()));
        }
        resampler->SetReferenceImage(fixedImage);
        resampler->UseReferenceImageOn();
        resampler->Update();

        RealImage::Pointer resampled = resampler->GetOutput();
        imageIO.WriteImage(args[2], resampled);
        imageIO.WriteSingleTransform(args[3].c_str(), transform);
    }
}
