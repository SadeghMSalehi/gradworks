//
//  piImageSlice.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/23/13.
//
//

#ifndef __ParticleGuidedRegistration__piImageSlice__
#define __ParticleGuidedRegistration__piImageSlice__

#include <iostream>


#include "piImageDef.h"

#include "itkRGBAPixel.h"
#include "itkScalarToARGBColormapImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkTransform.h"
#include "itkAffineTransform.h"
#include "vtkMatrix4x4.h"

#include "QPixmap"

class QGraphicsPixmapItem;


/**
 * this is a class representing a slice in QT gui
 */
namespace pi {
    typedef itk::RGBAPixel<unsigned char> RGBAPixel;
    typedef itk::Image<RGBAPixel, 2> RGBAImageType;
    typedef itk::Image<RGBAPixel, 3> RGBAVolumeType;

    class ImageStat {
    public:
        DataReal minIntensity, maxIntensity;
    };

    template <class T>
    class ImageReslice {
    public:
        typedef itk::AffineTransform<double> TransformType;

        ImageReslice() {
        }

        ~ImageReslice() {}

        void SetVTKTransform(vtkMatrix4x4* mat) {
            TransformType::Pointer affineTransform = TransformType::New();
            TransformType::ParametersType params, center;
            params.SetSize(12);
            center.SetSize(3);

            RealImage::SizeType srcSize = srcImg->GetBufferedRegion().GetSize();
            IntIndex centerIdx;
            for (int i = 0; i < 3; i++) {
                centerIdx[i] = srcSize[i] / 2.0;
            }
            RealImage::PointType centerPoint;
            srcImg->TransformIndexToPhysicalPoint(centerIdx, centerPoint);
            for (int i = 0; i < 3; i++) {
                center[i] = centerPoint[i];
            }

            params.Fill(0);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    params[4*i+j] = mat->GetElement(i, j);
                }
            }

            affineTransform->SetParameters(params);
            affineTransform->SetFixedParameters(center);

            m_transform = affineTransform;
        }
        
        // 2d slice related functions
        void SetImage(typename T::Pointer img) {
            srcImg = img;
            resampleRegion = srcImg->GetBufferedRegion();
            typedef itk::StatisticsImageFilter<RealImage> StatFilter;
            StatFilter::Pointer statFilter = StatFilter::New();
            statFilter->SetInput(srcImg);
            statFilter->Update();
            m_stat.minIntensity = statFilter->GetMinimum();
            m_stat.maxIntensity = statFilter->GetMaximum();
        }

        typename T::Pointer GetResampled() {
            if (resampledImg.IsNull()) {
                return typename T::Pointer();
            }

            typedef itk::ResampleImageFilter<RealImage,RealImage> ResampleFilter;
            typename ResampleFilter::Pointer resampleFilter = ResampleFilter::New();

            resampleFilter->SetInput(srcImg);
            resampleFilter->SetReferenceImage(resampledImg);
            resampleFilter->UseReferenceImageOn();
            if (m_transform.IsNotNull()) {
                resampleFilter->SetTransform(m_transform);
            }
            resampleFilter->Update();
            resampledImg = resampleFilter->GetOutput();

            UpdateVisualImage();
            return resampledImg;
        }

        QPixmap GetPixmap() {
            if (rgbaImage.IsNull()) {
                return QPixmap();
            }
            RGBAVolumeType::SizeType bitmapSz = rgbaImage->GetBufferedRegion().GetSize();
            QImage qImg = QImage((unsigned char*) rgbaImage->GetBufferPointer(),
                                 bitmapSz[0], bitmapSz[1], QImage::Format_ARGB32);
            return QPixmap::fromImage(qImg);
        }

        void UpdateVisualImage() {
            typedef itk::ScalarToARGBColormapImageFilter<T, RGBAVolumeType> ScalarToRGBFilter;
            typename ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
            rgbFilter->SetInput(resampledImg);
            rgbFilter->UseManualScalingOn();
            rgbFilter->SetMinimumValue(m_stat.minIntensity);
            rgbFilter->SetMaximumValue(m_stat.maxIntensity);
            rgbFilter->SetAlphaValue(255);
            rgbFilter->Update();
            rgbaImage = rgbFilter->GetOutput();
        }

        void SelectSlice(int axis, int index) {
            if (srcImg.IsNull()) {
                return;
            }
            resampleRegion = srcImg->GetBufferedRegion();
            typename T::IndexType idx1 = resampleRegion.GetIndex();
            typename T::IndexType idx2 = resampleRegion.GetUpperIndex();
            idx1[axis] = index;
            idx2[axis] = index;
            resampleRegion.SetIndex(idx1);
            resampleRegion.SetUpperIndex(idx2);
            typename ExtractFilterType::Pointer filter = ExtractFilterType::New();
            filter->SetInput(srcImg);
            filter->SetExtractionRegion(resampleRegion);
            filter->Update();
            resampledImg = filter->GetOutput();
            UpdateVisualImage();
        }

    private:
        typedef itk::ExtractImageFilter<T,T> ExtractFilterType;
        
        bool isModified;
        typename T::Pointer srcImg;
        typename T::Pointer resampledImg;
        typename T::RegionType resampleRegion;
        TransformType::Pointer m_transform;
        RGBAVolumeType::Pointer rgbaImage;
        ImageStat m_stat;
    };
    
    template <class T>
    class ImageSlice {
    public:
        int alpha;
        QGraphicsPixmapItem* pixmapCache;
        
        ImageSlice(): alpha(255), pixmapCache(NULL) {}
        
        void SetImage(typename T::Pointer label) {
            labelImage = label;
            typedef itk::ScalarToARGBColormapImageFilter<T, RGBAImageType> ScalarToRGBFilter;
            typename ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
            rgbFilter->SetInput(labelImage);
            rgbFilter->UseManualScalingOff();
            rgbFilter->SetAlphaValue(alpha);
            rgbFilter->Update();
            rgbaImage = rgbFilter->GetOutput();
        }

        typename T::Pointer GetImage() {
            return labelImage;
        }

        QPixmap GetPixmap() {
            if (rgbaImage.IsNull()) {
                return QPixmap();
            }
            RGBAImageType::SizeType bitmapSz = rgbaImage->GetBufferedRegion().GetSize();
            QImage qImg = QImage((unsigned char*) rgbaImage->GetBufferPointer(),
                                 bitmapSz[0], bitmapSz[1], QImage::Format_ARGB32);
            return QPixmap::fromImage(qImg);
        }
    private:
        typename T::Pointer labelImage;
        RGBAImageType::Pointer rgbaImage;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piImageSlice__) */
