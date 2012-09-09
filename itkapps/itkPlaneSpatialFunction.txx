/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPlaneSpatialFunction.txx,v $
  Language:  C++
  Date:      $Date: 2010/09/23 21:08:25 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPlaneSpatialFunction_txx
#define __itkPlaneSpatialFunction_txx

#include "itkPlaneSpatialFunction.h"

namespace itk
{

template <unsigned int VImageDimension,typename TInput>
PlaneSpatialFunction<VImageDimension,TInput>
::PlaneSpatialFunction()
{
  m_Normal[0] = m_Normal[1] = m_Normal[2] = 0;
  m_Origin.Fill(0.0);
}

template <unsigned int VImageDimension,typename TInput>
PlaneSpatialFunction<VImageDimension,TInput>
::~PlaneSpatialFunction()
{

}

template <unsigned int VImageDimension,typename TInput>
typename PlaneSpatialFunction<VImageDimension,TInput>::OutputType
PlaneSpatialFunction<VImageDimension,TInput>
::Evaluate(const InputType& position) const
{
  double acc = 0;
  for(unsigned int i = 0; i < VImageDimension; i++)
  {
    acc += m_Normal[i] * (position[i] - m_Origin[i]);
  }

  if (acc <= 0) // inside the plane 
  {
    return 1;
  } else {
    return 0; // outside the plane
  }
}

template <unsigned int VImageDimension,typename TInput>
void
PlaneSpatialFunction<VImageDimension,TInput>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  unsigned int i;
  os << indent << "Origin: [";
  for (i=0; i < VImageDimension - 1; i++) {
    os << m_Origin[i] << ", ";
  }
  os << "]" << std::endl;
  
  os << indent << "Normal: [";
  for (i=0; i < VImageDimension - 1; i++) {
    os << m_Normal[i] << ", ";
  }
  os << "]" << std::endl;

}

} // end namespace itk

#endif
