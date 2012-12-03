//
//  itkResampler.h
//  itkcmds
//
//  Created by Joohwi Lee on 9/14/12.
//
//

#ifndef itkcmds_itkResampler_h
#define itkcmds_itkResampler_h

#include "iostream"
#include "fstream"
#include "vector"
#include "itkImageIO.h"
#include "itkExceptionObject.h"
#include "itkMath.h"
#include "itkTimer.h"
#include "itkMathCode.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVersorRigid3DTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkSimilarity3DTransform.h"
#include "itkImageRegionIteratorWithIndex.h"

template <class TImage, class TLabel, class TTransform>
class itkResampler {
public:
    typedef itkMathCode<TImage,TTransform> MathCodeType;
    typedef typename itkMathCode<TImage,TTransform>::Matrix MatrixType;
    typedef typename itkMathCode<TImage,TTransform>::Vector VectorType;
    typedef itk::ImageRegionIteratorWithIndex<TImage> IteratorType;
    typedef typename TLabel::Pointer LabelPointer;

    itkResampler() {
    }

    void LoadImage(const char* src, const char* ref) {
        itkcmds::itkImageIO<TImage> _imageIO;
        _srcImg = _imageIO.ReadImageT(src);
        _refImg = _imageIO.ReadImageT(ref);
    }

    void SetInput(typename TImage::Pointer input) {
        _srcImg = input;
        _refImg = input;
    }

    void SetReference(typename TImage::Pointer ref) {
        _refImg = ref;
    }

    void SetTransform(typename TTransform::Pointer transform) {
        _transform = transform;
    }

    void SetInputMask(LabelPointer srcMask) {
        _srcMask = srcMask;
        
        if (_srcImg->GetBufferedRegion().GetSize() != _srcMask->GetBufferedRegion().GetSize()) {
            cout << "Regions are mismatch between image and mask" << endl;
            return;
        }
    }

    typename TTransform::Pointer GetTransform() {
        return _transform;
    }

    void Project() {
        _outImg = _imageIO.NewImageT(_refImg);

        MathCodeType mathCode;
        MatrixType imageTransform;
        mathCode.createIndex2IndexTransformMatrix(_srcImg, _transform, _refImg, imageTransform);

        typename TImage::RegionType outRegion = _outImg->GetBufferedRegion();
        IteratorType iter(_srcImg, _srcImg->GetBufferedRegion());
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            typename TImage::IndexType srcIdx = iter.GetIndex();
            typename TImage::IndexType dstIdx;
            // itk::ContinuousIndex<double> dstIdx;
            mathCode.transformIndex(imageTransform, srcIdx, dstIdx);

            //cout << dstIdx << endl;
            if (outRegion.IsInside(dstIdx)) {
                _outImg->SetPixel(dstIdx, iter.Get());
            }
        }
    }

    void Resample() {
        _outImg = _imageIO.NewImageT(_refImg);

        MathCodeType mathCode;
        MatrixType imageTransform;
        mathCode.createIndex2IndexTransformMatrix(_srcImg, _transform, _refImg, imageTransform);

        typedef itk::LinearInterpolateImageFunction<TImage> InterpolatorType;
        typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
        interpolator->SetInputImage(_srcImg);

        typename TImage::RegionType outRegion = _outImg->GetBufferedRegion();
        IteratorType iter(_outImg, _outImg->GetBufferedRegion());
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            typename TImage::IndexType outIdx = iter.GetIndex();
            typename InterpolatorType::ContinuousIndexType srcIdx;
            mathCode.template transformContinousIndex<typename InterpolatorType::ContinuousIndexType>(imageTransform, outIdx, srcIdx);

            typename TImage::IndexType maskIdx;
            for (int i = 0; i < TImage::ImageDimension; i++) {
                maskIdx[i] = (int) (srcIdx[i] + 0.5);
            }

            bool insideBuffer = interpolator->IsInsideBuffer(srcIdx);
            bool insideROI = _srcMask.IsNull() || (_srcMask->GetPixel(maskIdx) > 0);

            if (insideBuffer && insideROI) {
                typename TImage::PixelType pixel = interpolator->EvaluateAtContinuousIndex(srcIdx);
                iter.Set(pixel);
            }
        }
    }


    void ResampleBackward() {
        if (_refImg.IsNull()) {
            _outImg = _imageIO.NewImageT(_srcImg);
        } else {
            _outImg = _imageIO.NewImageT(_refImg);
        }

        MathCodeType mathCode;
        MatrixType imageTransform, imageTransformInverse;
        mathCode.createIndex2IndexTransformMatrix(_srcImg, _transform, _refImg, imageTransform);
        imageTransform.inverse(imageTransformInverse);

        typedef itk::LinearInterpolateImageFunction<TImage> InterpolatorType;
        typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
        interpolator->SetInputImage(_srcImg);
        
        typename TImage::RegionType srcRegion = _srcImg->GetBufferedRegion();
        IteratorType iter(_outImg, _outImg->GetBufferedRegion());
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            typename TImage::IndexType outIdx = iter.GetIndex();
            typename InterpolatorType::ContinuousIndexType srcIdx;
            mathCode.template transformContinousIndex<typename InterpolatorType::ContinuousIndexType>(imageTransformInverse, outIdx, srcIdx);

            typename TImage::IndexType maskIdx;
            for (int i = 0; i < TImage::ImageDimension; i++) {
                maskIdx[i] = (int) (srcIdx[i] + 0.5);
            }

            bool insideBuffer = interpolator->IsInsideBuffer(srcIdx);
            bool insideROI = _srcMask.IsNull() || (_srcMask->GetPixel(maskIdx) > 0);

            if (insideBuffer && insideROI) {
                typename TImage::PixelType pixel = interpolator->EvaluateAtContinuousIndex(srcIdx);
                iter.Set(pixel);
            }
        }
    }

    typename TImage::Pointer GetOutput() {
        return _outImg;
    }

    void Run(int argc, const char* argv[]) {
        LoadImage(argv[1], argv[2]);
        ResampleBackward();
        _imageIO.WriteImageT(argv[3], _outImg);
    }

private:
    itkcmds::itkImageIO<TImage> _imageIO;
    typename TImage::Pointer _srcImg;
    typename TImage::Pointer _refImg;
    typename TImage::Pointer _outImg;
    typename TTransform::Pointer _transform;
    LabelPointer _srcMask;
    unsigned int _srcLabel;
};

#endif