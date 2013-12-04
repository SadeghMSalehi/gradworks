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
#ifndef __itkEntropyImageToImageMetricv4_hxx
#define __itkEntropyImageToImageMetricv4_hxx

#include "itkEntropyImageToImageMetricv4.h"
#include "piImageIO.h"

/**
 * This class has two separated threaders.
 * One is used in the initialization step, and another is used in the computation step.
 * The initialization step computes the covariance matrix and average of intensities.
 * Once the computation of the inverse of covariance is done, the derivatives are computed from those.
 */
namespace itk
{

    template <class T>
    const int ComputeNumberOfPixels(T& region) {
        int nPixels = 1;
        for (int i = 0; i < T::ImageDimension; i++) {
            nPixels *= region.GetSize(i);
        }
        return nPixels;
    }

    template < class TFixedImage, class TMovingImage, class TVirtualImage >
    EntropyImageToImageMetricv4<TFixedImage,TMovingImage,TVirtualImage>
    ::EntropyImageToImageMetricv4()
    {
        // We have our own GetValueAndDerivativeThreader's that we want
        // ImageToImageMetricv4 to use.
        this->m_DenseGetValueAndDerivativeThreader  = EntropyDenseGetValueAndDerivativeThreaderType::New();
        this->m_SparseGetValueAndDerivativeThreader = EntropySparseGetValueAndDerivativeThreaderType::New();

        this->m_HelperDenseThreader  = EntropyHelperDenseThreaderType::New();
        this->m_HelperSparseThreader = EntropyHelperSparseThreaderType::New();

        if( this->m_MovingTransform->GetTransformCategory() == MovingTransformType::DisplacementField )
        {
            itkExceptionMacro("does not support displacement field transforms!!");
        }
    }

    template < class TFixedImage, class TMovingImage, class TVirtualImage >
    EntropyImageToImageMetricv4<TFixedImage,TMovingImage,TVirtualImage>
    ::~EntropyImageToImageMetricv4()
    {
    }

    template < class TFixedImage, class TMovingImage, class TVirtualImage  >
    void
    EntropyImageToImageMetricv4<TFixedImage,TMovingImage,TVirtualImage>
    ::PrintSelf(std::ostream& os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
    }


    template < class TFixedImage, class TMovingImage, class TVirtualImage >
    void
    EntropyImageToImageMetricv4<TFixedImage,TMovingImage,TVirtualImage>
    ::InitializeForIteration() const
    {

        Superclass::InitializeForIteration();

        this->m_AverageFix = NumericTraits<MeasureType>::Zero;
        this->m_AverageMov = NumericTraits<MeasureType>::Zero;

        // compute the average intensity of the sampled pixels
        // Invoke the pipeline in the helper threader
        // refer to DomainThreader::Execute()

        if( this->m_UseFixedSampledPointSet ) // sparse sampling
        {
            // not supported
            itkExceptionMacro("Sparse sampling is not supported");
            SizeValueType numberOfPoints = this->GetNumberOfDomainPoints();
            if( numberOfPoints < 1 )
            {
                itkExceptionMacro("FixedSampledPointSet must have 1 or more points.");
            }
            typename ImageToImageMetricv4GetValueAndDerivativeThreader< ThreadedIndexedContainerPartitioner, Self >::DomainType range;
            range[0] = 0;
            range[1] = numberOfPoints - 1;
            this->m_HelperSparseThreader->Execute( const_cast< Self* >(this), range );
        }
        else // dense sampling
        {
            AllocateData();

            // Helper thread will fill m_Data variable
            // by applying transfrom from the virtual region to the moving image region

            this->m_HelperDenseThreader->Execute( const_cast< Self* >(this), this->GetVirtualRegion() );
        }

        /*
         * the results:
         *  this->m_AverageFix
         *  this->m_AverageMov
         * will be stored during helper::AfterThreadedExecution()
         */


    }

    template < class TFixedImage, class TMovingImage, class TVirtualImage >
    void
    EntropyImageToImageMetricv4<TFixedImage,TMovingImage,TVirtualImage>
    ::AllocateData() const {
        // compute the number of pixels in the virtual region
        // Each pixel has to be tested if it is valid or not
        // because some of region is out of the moving image region

        typename VirtualImageType::RegionType region = this->GetVirtualRegion();

        pi::ImageIO<MarkImageType> markImageIO;
        m_Mark = markImageIO.template NewImageS<VirtualImageType>(this->m_VirtualImage);

        pi::ImageIO<VirtualImageType> imageIO;
        m_Data.push_back(imageIO.NewImage(this->m_VirtualImage));
        m_Data.push_back(imageIO.NewImage(this->m_VirtualImage));
    }


    template < class TFixedImage, class TMovingImage, class TVirtualImage >
    void
    EntropyImageToImageMetricv4<TFixedImage,TMovingImage,TVirtualImage>
    ::ComputeCovariance(int iter) {
        // just for two images
        const int nImages = m_Data.size();

        // a pointer to the pixel validity marker
        bool* markPtr = m_Mark->GetBufferPointer();

        // the number of elements in the virtual domain
        const int nElems = m_Mark->GetPixelContainer()->Size();


        // loop over all elements in the virtual domain
        for (int k = 0; k < nElems; k++) {
            // only compute when the pixel overlaps each other
            if (markPtr[k]) {
                // compute the sum of pixels
                double pixelSum = 0;
                for (int i = 0; i < nImages; i++) {
                    pixelSum += m_Data[i]->GetBufferPointer()[k];
                }
                double pixelAvg = pixelSum / nImages;
                for (int i = 0; i < nImages; i++) {
                    // subtract average from the data
                    m_Data[i]->GetBufferPointer()[k] -= pixelAvg;
                }
            }
        }

        // loop over the i-th row of Y' matrix
        for (int i = 0; i < nImages; i++) {
            FixedImagePixelType* iPtr = m_Data[i]->GetBufferPointer();

            // loop over the j-th column of Y matrix
            for (int j = i; j < nImages; j++) {
                MovingImagePixelType* jPtr = m_Data[j]->GetBufferPointer();
                double sum = 0;
                for (int k = 0; k < nElems; k++) {
                    // check if k-th voxel is valid
                    if (markPtr[k]) {
                        // data is not averaged
                        sum += (iPtr[k] * jPtr[k]);
                    }
                }
                m_Covariance[i][j] = sum;
            }
        }

        // m_Covariance is an upper triangular matrix
        // because elements (i > j) are empty
        // Then, make m_Covarianece a full matrix by filling in the lower triangle
        for (int i = 0; i < nImages; i++) {
            for (int j = i + 1; j < nImages; j++) {
                m_Covariance[j][i] = m_Covariance[i][j];
            }
        }

        // compute the inverse of covariance
        m_InverseCovariance = vnl_matrix_inverse<double>(m_Covariance);
    }

} // end namespace itk


#endif
