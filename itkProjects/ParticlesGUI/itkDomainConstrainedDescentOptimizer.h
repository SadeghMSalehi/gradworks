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
#ifndef __itkDomainConstrainedDescentOptimizer_h
#define __itkDomainConstrainedDescentOptimizer_h

#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkVersor.h"
#include "myImageContainer.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"

class myImplicitSurfaceConstraint;

namespace itk
{


    /** \class DomainConstrainedDescentOptimizer
     * \brief Implement a gradient descent optimizer for the VersorRigid3DTransform
     * parameter space.
     *
     * DomainConstrainedDescentOptimizer is a variant of the gradient descent
     * optimizer implmented in RegularStepGradientDescentOptimizer.
     *
     * Versors are not in a vector space, for that reason, the classical gradient
     * descent algorithm has to be modified in order to be applicable to Versors
     * (unit quaternions) that form the group SO(3).
     *
     * The Versor space has only three degrees of freedom, even though Versors are
     * represented using four values.
     *
     * This optimizer assumes that the CostFunction to be optimized has an
     * itk::Versor and an itk::Vector as parameters.
     *
     * \sa RegularStepGradientDescentOptimizer
     * \sa Versor
     * \sa VersorRigid3DTransform
     *
     * \ingroup Numerics Optimizers
     * \ingroup ITKOptimizers
     */
    class ITK_EXPORT DomainConstrainedDescentOptimizer:
    public RegularStepGradientDescentBaseOptimizer
    {
    public:
        /** Standard class typedefs. */
        typedef DomainConstrainedDescentOptimizer         Self;
        typedef RegularStepGradientDescentBaseOptimizer Superclass;
        typedef SmartPointer< Self >                    Pointer;
        typedef SmartPointer< const Self >              ConstPointer;

        typedef SignedDanielssonDistanceMapImageFilter<SliceType, SliceType> DistanceFilterType;
        typedef DistanceFilterType::VectorImageType DistanceVectorImageType;


        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(DomainConstrainedDescentOptimizer,
                     RegularStepGradientDescentBaseOptimizer);

        /** Advance one step following the gradient direction. */
        virtual void StepAlongGradient(double factor,
                                       const DerivativeType & transformedGradient);

        void SetConstraint(myImplicitSurfaceConstraint* constraint);

    protected:
        DomainConstrainedDescentOptimizer() {}
        virtual ~DomainConstrainedDescentOptimizer() {}
    private:
        DomainConstrainedDescentOptimizer(const Self &); //purposely not implemented
        void operator=(const Self &);                  //purposely not implemented

        myImplicitSurfaceConstraint* m_Constraint;
    };
} // end namespace itk

#endif