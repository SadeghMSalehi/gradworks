//
//  itkMyTypes.h
//  itktools
//
//  Created by Joohwi Lee on 9/20/12.
//
//

#ifndef itktools_itkMyTypes_h
#define itktools_itkMyTypes_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMyScaleVersor3DTransformOptimizer.h"
#include "itkMySlice.h"
#include "itkResampleImageFilter.h"

typedef itk::Image<double, 3> ImageType;
typedef itk::Image<double, 2> ImageSliceType;
typedef itk::Image<unsigned short, 3> LabelType;
typedef itk::Image<unsigned short, 2> LabelSliceType;
typedef itk::Image<unsigned int, 2> BitmapType;
typedef itk::MatrixOffsetTransformBase<double, 3> MatrixTransformType;
typedef itkMySlice<ImageType, ImageSliceType> ImageSlicer;
typedef itkMySlice<LabelType, LabelSliceType> LabelSlicer;
typedef itk::ResampleImageFilter<LabelType,LabelType> ResampleFilter;

typedef enum ImageKindEnum { Source = 0, Target = 1, Label = 2, TransformedLabel = 3 } ImageKindType;

static int __colorMetallicRainbow[] = { 0x00000000, 0xc9170b, 0xc9760b, 0xbdc90b, 0x0bc917, 0x0bbdc9, 0x170bc9, 0x760bc9 };

class GrayToLabelFunctor {
private:
    int _zeroAlpha;
    int _insideAlpha;
    int _outsideAlpha;

public:
    GrayToLabelFunctor() {
        _zeroAlpha = 0x00000000;
        _insideAlpha = 0x33000000;
        _outsideAlpha = 0x88000000;
    }
    bool operator!=(const GrayToLabelFunctor &) const
    {
        return false;
    }

    bool operator==(const GrayToLabelFunctor & other) const
    {
        return !( *this != other );
    }

    void SetAlpha(int zeroAlpha, int insideAlpha, int outsideAlpha) {
        _zeroAlpha = (zeroAlpha << 24) & 0xff000000;
        SetInsideAlpha(insideAlpha);
        _outsideAlpha = (outsideAlpha << 24) & 0xff000000;
    }

    void SetInsideAlpha(int insideAlpha) {
        _insideAlpha = (insideAlpha << 24) & 0xff000000;
    }
    
    inline int operator()(const unsigned short & A) const
    {
        if (A == 0) {
            return _zeroAlpha;
        } else if (A > 0 && A < 8) {
            return _insideAlpha | __colorMetallicRainbow[A];
        }
        return _insideAlpha;
    }
};

class GrayToRGBFunctor {
private:
    int _zeroAlpha;
    int _insideAlpha;
    int _outsideAlpha;
    int _mini, _maxi, _meani, _stdi;
    double _scale;
    double _shift;
public:
    void SetAlpha(int zeroAlpha, int insideAlpha, int outsideAlpha) {
        _zeroAlpha = (zeroAlpha << 24) & 0xff000000;
        SetInsideAlpha(insideAlpha);
        _outsideAlpha = (outsideAlpha << 24) & 0xff000000;
    }

    void SetInsideAlpha(int insideAlpha) {
        _insideAlpha = (insideAlpha << 24) & 0xff000000;
    }

    void SetImageStatistics(int mini, int maxi, int meani, int stdi) {
        _mini = mini;
        _maxi = maxi;
        _meani = meani;
        _stdi = stdi;
        _shift = mini;
        _scale = 255.0 / (_maxi - _mini);
    }

    void SetShiftAndScale(double shift, double scale) {
        _shift = shift;
        _scale = scale;
    }

    bool operator!=(const GrayToRGBFunctor &) const
    {
        return false;
    }

    bool operator==(const GrayToRGBFunctor & other) const
    {
        return !( *this != other );
    }

    inline int operator()(const unsigned short & A) const
    {
        int a = (int) ((A - _shift) * _scale);
        return (a | a << 8 | a << 16 | 0xff << 24);
    }
};

typedef itk::UnaryFunctorImageFilter<LabelSliceType, BitmapType, GrayToLabelFunctor> LabelMapFilter;
typedef itk::UnaryFunctorImageFilter<ImageSliceType, BitmapType, GrayToRGBFunctor> GrayImageFilter;

typedef itk::TransformBase TransformBase;
typedef itk::TransformFileReader TransformReaderType;
typedef itk::TransformFileWriter TransformWriterType;
typedef itk::NearestNeighborInterpolateImageFunction<LabelType> InterpolatorNN;
typedef itk::LinearInterpolateImageFunction<ImageType> LinearImageInterpolator;
typedef itk::NearestNeighborInterpolateImageFunction<ImageType> NNImageInterpolator;
typedef itk::RescaleIntensityImageFilter<ImageType, LabelType> IntensityFilter;




#endif