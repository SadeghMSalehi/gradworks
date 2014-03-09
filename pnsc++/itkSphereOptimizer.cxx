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
#ifndef _itkSphereOptimizer_hxx
#define _itkSphereOptimizer_hxx

#include "itkSphereOptimizer.h"

namespace itk
{
    /**
     * Advance one Step following the gradient direction
     * This method will be overrided in non-vector spaces
     */
    void
    SphereOptimizer::StepAlongGradient(double factor,
                        const DerivativeType & transformedGradient)
    {
        const ParametersType & currentPosition = this->GetCurrentPosition();

        ParametersType newParameters(SpaceDimension);

        // Now do the typical update for a Vector space.
        for ( unsigned int k = 0; k < currentPosition.GetSize(); k++ )
        {
            newParameters[k] = currentPosition[k] + transformedGradient[k] * factor;
        }

        double mag = 0;
        for (int k = 0; k < currentPosition.GetSize() - 1; k++) {
            mag += (newParameters[k] * newParameters[k]);
        }
        mag = sqrt(mag);
        for (int k = 0; k < currentPosition.GetSize() - 1; k++) {
            newParameters[k] = newParameters[k] / mag;
        }
        
        this->SetCurrentPosition(newParameters);
    }
} // end namespace itk

#endif
