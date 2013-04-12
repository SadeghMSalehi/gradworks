/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkEntropyImageToImageMetricv4_h
#define __itkEntropyImageToImageMetricv4_h

#include "itkImageToImageMetricv4.h"

#include "itkEntropyImageToImageMetricv4GetValueAndDerivativeThreader.h"
#include "itkEntropyImageToImageMetricv4HelperThreader.h"

namespace itk
{


    template <class TFixedImage, class TMovingImage, class TVirtualImage = TFixedImage >
    class ITK_EXPORT EntropyImageToImageMetricv4 :
    public ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage>
    {
    public:
        /** Standard class typedefs. */
        typedef EntropyImageToImageMetricv4                                Self;
        typedef ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage> Superclass;
        typedef SmartPointer<Self>                                             Pointer;
        typedef SmartPointer<const Self>                                       ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(EntropyImageToImageMetricv4, ImageToImageMetricv4);

        /** Superclass types */
        typedef typename Superclass::MeasureType             MeasureType;
        typedef typename Superclass::DerivativeType          DerivativeType;

        typedef typename Superclass::FixedImagePointType     FixedImagePointType;
        typedef typename Superclass::FixedImagePixelType     FixedImagePixelType;
        typedef typename Superclass::FixedImageGradientType  FixedImageGradientType;

        typedef typename Superclass::MovingImagePointType    MovingImagePointType;
        typedef typename Superclass::MovingImagePixelType    MovingImagePixelType;
        typedef typename Superclass::MovingImageGradientType MovingImageGradientType;

        typedef typename Superclass::MovingTransformType        MovingTransformType;
        typedef typename Superclass::JacobianType               JacobianType;
        typedef typename Superclass::VirtualImageType           VirtualImageType;
        typedef typename Superclass::VirtualIndexType           VirtualIndexType;
        typedef typename Superclass::VirtualPointType           VirtualPointType;
        typedef typename Superclass::VirtualPointSetType        VirtualPointSetType;

        /* Image dimension accessors */
        itkStaticConstMacro(VirtualImageDimension, ImageDimensionType,
                            TVirtualImage::ImageDimension);
        itkStaticConstMacro(FixedImageDimension, ImageDimensionType,
                            TFixedImage::ImageDimension);
        itkStaticConstMacro(MovingImageDimension, ImageDimensionType,
                            TMovingImage::ImageDimension);

    protected:
        EntropyImageToImageMetricv4();
        virtual ~EntropyImageToImageMetricv4();

        /** Perform any initialization required before each evaluation of
         * \c GetValueAndDerivative. This is distinct from Initialize, which
         * is called only once before a number of iterations, e.g. before
         * a registration loop. */
        virtual void InitializeForIteration() const;

        friend class ImageToImageMetricv4GetValueAndDerivativeThreaderBase< ThreadedImageRegionPartitioner< Superclass::VirtualImageDimension >, Self >;
        friend class ImageToImageMetricv4GetValueAndDerivativeThreaderBase< ThreadedIndexedContainerPartitioner, Self >;
        friend class ImageToImageMetricv4GetValueAndDerivativeThreader< ThreadedImageRegionPartitioner< Superclass::VirtualImageDimension >, Self >;
        friend class ImageToImageMetricv4GetValueAndDerivativeThreader< ThreadedIndexedContainerPartitioner, Self >;

        friend class EntropyImageToImageMetricv4GetValueAndDerivativeThreader< ThreadedImageRegionPartitioner< Superclass::VirtualImageDimension >, Superclass, Self >;
        friend class EntropyImageToImageMetricv4GetValueAndDerivativeThreader< ThreadedIndexedContainerPartitioner, Superclass, Self >;
        typedef EntropyImageToImageMetricv4GetValueAndDerivativeThreader< ThreadedImageRegionPartitioner< Superclass::VirtualImageDimension >, Superclass, Self >
        EntropyDenseGetValueAndDerivativeThreaderType;
        typedef EntropyImageToImageMetricv4GetValueAndDerivativeThreader< ThreadedIndexedContainerPartitioner, Superclass, Self >
        EntropySparseGetValueAndDerivativeThreaderType;

        friend class EntropyImageToImageMetricv4HelperThreader< ThreadedImageRegionPartitioner< Superclass::VirtualImageDimension >, Superclass, Self >;
        friend class EntropyImageToImageMetricv4HelperThreader< ThreadedIndexedContainerPartitioner, Superclass, Self >;

        typedef EntropyImageToImageMetricv4HelperThreader< ThreadedImageRegionPartitioner< Superclass::VirtualImageDimension >, Superclass, Self >
        EntropyHelperDenseThreaderType;
        typedef EntropyImageToImageMetricv4HelperThreader< ThreadedIndexedContainerPartitioner, Superclass, Self >
        EntropyHelperSparseThreaderType;

        typename EntropyHelperDenseThreaderType::Pointer  m_HelperDenseThreader;
        typename EntropyHelperSparseThreaderType::Pointer m_HelperSparseThreader;
        
        /* These values are computed during InitializeForIteration(),
         * using the helper class
         * */
        mutable MeasureType m_AverageFix;
        mutable MeasureType m_AverageMov;
        
        void PrintSelf(std::ostream& os, Indent indent) const;
        
    private:
        EntropyImageToImageMetricv4(const Self &); //purposely not implemented
        void operator = (const Self &); //purposely not implemented
    };
    
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEntropyImageToImageMetricv4.hxx"
#endif

#endif
