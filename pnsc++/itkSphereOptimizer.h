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
#ifndef __itkSphereOptimizer_h
#define __itkSphereOptimizer_h

#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkVersor.h"

namespace itk
{
/** \class SphereOptimizer
 * \brief Implement a gradient descent optimizer for the VersorRigid3DTransform
 * parameter space.
 *
 * SphereOptimizer is a variant of the gradient descent
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
class ITK_EXPORT SphereOptimizer:
  public RegularStepGradientDescentBaseOptimizer
{
public:
  /** Standard class typedefs. */
  typedef SphereOptimizer         Self;
  typedef RegularStepGradientDescentBaseOptimizer Superclass;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SphereOptimizer,
               RegularStepGradientDescentBaseOptimizer);

  /** This class is specialized for 3D  */
  itkStaticConstMacro(SpaceDimension, unsigned int, 4);

  /**  Versor Type  */
  typedef Versor< double >       VersorType;
  typedef VersorType::VectorType VectorType;

  /** Advance one step following the gradient direction. */
  virtual void StepAlongGradient(double factor,
                                 const DerivativeType & transformedGradient);

protected:
  SphereOptimizer() {}
  virtual ~SphereOptimizer() {}
private:
  SphereOptimizer(const Self &); //purposely not implemented
  void operator=(const Self &);                  //purposely not implemented
};
} // end namespace itk

#endif
