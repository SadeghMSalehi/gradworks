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

/**
 * This class has two separated threaders.
 * One is used in the initialization step, and another is used in the computation step.
 * The initialization step computes the covariance matrix and average of intensities.
 * Once the computation of the inverse of covariance is done, the derivatives are computed from those.
 */
namespace itk
{

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

        const int nPixels = 0;

        this->m_Averages.set_size(2);
        this->m_Covariances.set_size(2, nPixels);

        // compute the average intensity of the sampled pixels
        // Invoke the pipeline in the helper threader
        // refer to DomainThreader::Execute()

        if( this->m_UseFixedSampledPointSet ) // sparse sampling
        {
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
            this->m_HelperDenseThreader->Execute( const_cast< Self* >(this), this->GetVirtualRegion() );
        }

        /*
         * the results:
         *  this->m_AverageFix
         *  this->m_AverageMov
         * will be stored during helper::AfterThreadedExecution()
         */


    }

} // end namespace itk


#endif
