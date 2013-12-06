//
//  piImageRxRunner.h
//  PxImageReg
//
//  Created by Joohwi Lee on 11/30/13.
//
//

#ifndef __PxImageReg__piImageRxRunner__
#define __PxImageReg__piImageRxRunner__

#include <iostream>

#include "piOptions.h"
#include "piImageDef.h"

#include <itkSingleValuedCostFunctionv4.h>

namespace pi {
    /**
     * @brief The TransformingImage class
     *
     * This class encapsulates the moving image.
     * The moving image changes its coordinate space depending on the transformation.
     * Since most similarity metrics requires the computation of derivatives,
     * this class supports functions to sample not only pixel values but also gradient values.
     */
    class TransformingImage
    {
    public:
        typedef std::vector<TransformingImage> Vector;
        typedef itk::LinearInterpolateImageFunction<RealImage> IntensityInterpolator;
        typedef itk::VectorLinearInterpolateImageFunction<GradientImage> GradientInterpolator;

        TransformingImage();

        // set a real-valued intensity image
        void setImage(RealImage::Pointer image);

        // a transform is given as a pointer
        void setTransform(TransformType* transform);

        // transform parameters are copied into this class
        void setTransformParameters(TransformType::ParametersType params);

        // sample a pixel value for the given fixed physical point
        RealImage::PixelType samplePixel(RealImage::PointType fixedPoint, bool& isValidPixel);

        // sample the gradient value for the given fixed physical point
        GradientImage::PixelType sampleGradient(RealImage::PointType fixedPoint, bool& isValidPixel);
        

        // resample the moving image into the fixed image frame
        RealImage::Pointer getResampledImage(RealImage::Pointer fixedImage);


        // inline functions
        inline TransformType* getTransform() {
            return _transform;
        }

        inline RealImage::Pointer getImage() { return _image; }

        inline IntensityInterpolator::Pointer getInterpolator() {
            return _intensityInterpolator;
        }


    public:
        RealImage::Pointer _image;
        GradientImage::Pointer _gradientImage;
        IntensityInterpolator::Pointer _intensityInterpolator;
        GradientInterpolator::Pointer _gradientInterpolator;
        TransformType* _transform;
    };


    /**
     * @brief The EntropyImageMetric class
     *
     */
    class EntropyImageMetric: public itk::SingleValuedCostFunctionv4 {
    public:
        typedef itk::SingleValuedCostFunctionv4 Superclass;
        typedef Superclass::DerivativeType DerivativeType;

        void setFixedImage(RealImage::Pointer image);
        void addMovingImage(RealImage::Pointer image, TransformType* transform, TransformType::ParametersType& params);
        virtual double GetValue() const;
        virtual void GetValueAndDerivative(double& value, DerivativeType& deriv) const;

    private:
        TransformingImage _fixedImage;
        TransformingImage::Vector _movingImages;
    };



    /**
     * @brief The ImageRx class
     */
    class ImageRx {
    public:
        ImageRx(Options& opts, StringVector& args) {}

        // main entry functions
        void mainRigidRegistration(Options& opts, StringVector& args);
    };

}
#endif /* defined(__PxImageReg__piImageRxRunner__) */
