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
#ifndef __itkEntropyImageToImageMetricv4HelperThreader_hxx
#define __itkEntropyImageToImageMetricv4HelperThreader_hxx

#include "itkEntropyImageToImageMetricv4HelperThreader.h"

namespace itk
{
    template<class TDomainPartitioner, class TImageToImageMetric, class TEntropyMetric>
    void EntropyImageToImageMetricv4HelperThreader< TDomainPartitioner, TImageToImageMetric, TEntropyMetric>
    ::BeforeThreadedExecution()
    {
        Superclass::BeforeThreadedExecution();

        /* Store the casted pointer to avoid dynamic casting in tight loops. */
        this->m_EntropyAssociate = dynamic_cast<TEntropyMetric *>(this->m_Associate);

        this->m_FixSumPerThread.resize(this->GetNumberOfThreadsUsed());
        this->m_MovSumPerThread.resize(this->GetNumberOfThreadsUsed());

        //---------------------------------------------------------------
        // Set initial values.
        for (ThreadIdType i = 0; i < this->GetNumberOfThreadsUsed(); i++)
        {
            this->m_FixSumPerThread[i] = NumericTraits<InternalComputationValueType>::Zero;
            this->m_MovSumPerThread[i] = NumericTraits<InternalComputationValueType>::Zero;
        }

        // Allocate data matrix and fill the variable with intensity values
    }

    template<class TDomainPartitioner, class TImageToImageMetric, class TEntropyMetric>
    void
    EntropyImageToImageMetricv4HelperThreader<TDomainPartitioner,
    TImageToImageMetric, TEntropyMetric>::AfterThreadedExecution()
    {
        static int iter = 0;

        /* Store the number of valid points the enclosing class \c
         * m_NumberOfValidPoints by collecting the valid points per thread. */
        this->m_EntropyAssociate->m_NumberOfValidPoints = NumericTraits<SizeValueType>::Zero;

        for (ThreadIdType i = 0; i < this->GetNumberOfThreadsUsed(); i++)
        {
            this->m_EntropyAssociate->m_NumberOfValidPoints += this->m_NumberOfValidPointsPerThread[i];
        }

        if (this->m_EntropyAssociate->m_NumberOfValidPoints <= 0 )
        {
            itkWarningMacro("collected only zero points");
            return;
        }

        InternalComputationValueType sumF = NumericTraits<InternalComputationValueType>::Zero;
        InternalComputationValueType sumM = NumericTraits<InternalComputationValueType>::Zero;

        for (size_t i = 0; i < this->m_MeasurePerThread.size(); i++)
        {
            sumF += this->m_FixSumPerThread[i];
            sumM += this->m_MovSumPerThread[i];
        }

        this->m_EntropyAssociate->m_AverageFix = sumF / this->m_EntropyAssociate->m_NumberOfValidPoints;
        this->m_EntropyAssociate->m_AverageMov = sumM / this->m_EntropyAssociate->m_NumberOfValidPoints;

        iter ++;

        this->m_EntropyAssociate->ComputeCovariance(iter);

    }

    template<class TDomainPartitioner, class TImageToImageMetric, class TEntropyMetric>
    bool
    EntropyImageToImageMetricv4HelperThreader<TDomainPartitioner,
    TImageToImageMetric, TEntropyMetric>
    ::ProcessVirtualPoint( const VirtualIndexType & virtualIndex, const VirtualPointType & virtualPoint, const ThreadIdType threadID )
    {
        FixedOutputPointType        mappedFixedPoint;
        FixedImagePixelType         mappedFixedPixelValue;
        MovingOutputPointType       mappedMovingPoint;
        MovingImagePixelType        mappedMovingPixelValue;
        bool                        pointIsValid = false;


        /* Transform the point into fixed and moving spaces, and evaluate.
         * Different behavior with pre-warping enabled is handled transparently.
         * Do this in a try block to catch exceptions and print more useful info
         * then we otherwise get when exceptions are caught in MultiThreader. */
        try
        {
            pointIsValid = this->m_EntropyAssociate->TransformAndEvaluateFixedPoint( virtualPoint, mappedFixedPoint, mappedFixedPixelValue );
        }
        catch( ExceptionObject & exc )
        {
            //NOTE: there must be a cleaner way to do this:
            std::string msg("Caught exception: \n");
            msg += exc.what();
            ExceptionObject err(__FILE__, __LINE__, msg);
            throw err;
        }
        if( !pointIsValid )
        {
            return pointIsValid;
        }

        try
        {
            pointIsValid = this->m_EntropyAssociate->TransformAndEvaluateMovingPoint( virtualPoint, mappedMovingPoint, mappedMovingPixelValue );
        }
        catch( ExceptionObject & exc )
        {
            std::string msg("Caught exception: \n");
            msg += exc.what();
            ExceptionObject err(__FILE__, __LINE__, msg);
            throw err;
        }
        if( !pointIsValid )
        {
            return pointIsValid;
        }

        /* Do the specific calculations for values */
        try
        {
            this->m_FixSumPerThread[threadID] += mappedFixedPixelValue;
            this->m_MovSumPerThread[threadID] += mappedMovingPixelValue;
        }
        catch( ExceptionObject & exc )
        {
            std::string msg("Exception in ProcessVirtualPoint:\n");
            msg += exc.what();
            ExceptionObject err(__FILE__, __LINE__, msg);
            throw err;
        }
        if( pointIsValid )
        {
            this->m_NumberOfValidPointsPerThread[threadID]++;
        }


        // update m_Mark image
        this->m_EntropyAssociate->m_Mark->SetPixel(virtualIndex, pointIsValid);
        this->m_EntropyAssociate->m_Data[0]->SetPixel(virtualIndex, mappedFixedPixelValue);
        this->m_EntropyAssociate->m_Data[1]->SetPixel(virtualIndex, mappedMovingPixelValue);

        return pointIsValid;
    }
} // end namespace itk

#endif
