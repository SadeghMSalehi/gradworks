/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPlaneSpatialFunction.h,v $
  Language:  C++
  Date:      $Date: 2010/09/23 21:08:25 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPlaneSpatialFunction_h
#define __itkPlaneSpatialFunction_h

#include "itkInteriorExteriorSpatialFunction.h"

namespace itk
{

/** \class PlaneSpatialFunction
 * \brief Spatial function implementation of a plane
 *
 * Implements a function that returns 0 for points inside or on the surface
 * of a plane, 1 for points outside the plane
 * 
 * \ingroup SpatialFunctions
 */
template <unsigned int VImageDimension=3,typename TInput=Point<double,VImageDimension> >
class ITK_EXPORT PlaneSpatialFunction
: public InteriorExteriorSpatialFunction<VImageDimension,TInput>
{
public:
  /** Standard class typedefs. */
  typedef PlaneSpatialFunction<VImageDimension,TInput>            Self;
  typedef InteriorExteriorSpatialFunction<VImageDimension,TInput> Superclass;
  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PlaneSpatialFunction,InteriorExteriorSpatialFunction);

  /** Input type for the function. */
  typedef typename Superclass::InputType InputType;

  /** Output type for the function. */
  typedef typename Superclass::OutputType OutputType;

  /** Evaluates the function at a given position */
  OutputType Evaluate(const InputType& position) const;

  /** Get and set the origin of the sphere. */
  itkGetConstMacro(Origin, InputType);
  itkSetMacro(Origin, InputType);
  
  /** Get and set the normal of the plane */
  double* GetNormal() {
    return m_Normal;
  }

  void SetNormal(double x, double y, double z) {
    m_Normal[0] = x;
    m_Normal[1] = y;
    m_Normal[2] = z;
  }
       
protected:
  PlaneSpatialFunction();
  virtual ~PlaneSpatialFunction();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  PlaneSpatialFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** The center of the sphere (of the same type as Input). */
  InputType m_Origin;

  /** The radius of the sphere. */ double m_Normal[3]; 
};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_PlaneSpatialFunction(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT PlaneSpatialFunction< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef PlaneSpatialFunction< ITK_TEMPLATE_2 x > \
                                           PlaneSpatialFunction##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkPlaneSpatialFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkPlaneSpatialFunction.txx"
#endif

#endif
