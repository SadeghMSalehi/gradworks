#ifndef PIMOVINGIMAGE_H
#define PIMOVINGIMAGE_H

#include <vector>
#include "piImageDef.h"

namespace pi {
    /**
     * @brief The MovingImage class
     *
     * This class encapsulates the moving image.
     * The moving image changes its coordinate space depending on the transformation.
     * Since most similarity metrics requires the computation of derivatives,
     * this class supports functions to sample not only pixel values but also gradient values.
     */
    class MovingImage
    {
    public:
        typedef std::vector<MovingImage> Vector;

        MovingImage();

        // set a real-valued intensity image
        void setImage(RealImage::Pointer image);

        // a transform is given as a pointer
        void setTransform(TransformType* transform);

        // transform parameters are copied into this class
        void setTransformParameters(TransformType::ParametersType params);

        // sample a pixel value for the given fixed physical point
        RealImage::PixelType samplePixel(RealImage::PointType fixedPoint);

        // sample the gradient value for the given fixed physical point
        GradientImage::PixelType sampleGradient(RealImage::PointType fixedPoint);

    public:
        GradientImage::Pointer _gradientImage;
        itk::LinearInterpolateImageFunction::Pointer _intensityInterpolator;
        itk::VectorLinearInterpolateImageFunction::Pointer _gradientInterpolator;

    };
}

#endif // PIMOVINGIMAGE_H
