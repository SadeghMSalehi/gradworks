/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorConnectedImageFunction.txx,v $
  Language:  C++
  Date:      $Date: 2010/09/23 21:08:25 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorConnectedImageFunction_txx
#define __itkTensorConnectedImageFunction_txx

#include "itkTensorConnectedImageFunction.h"
#include "itkNumericTraits.h"
#include "itkConstNeighborhoodIterator.h"
#include <iostream>

#include "itkDiffusionTensor3D.h"

typedef itk::DiffusionTensor3D<double> TensorType;

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TCoordRep>
TensorConnectedImageFunction<TInputImage,TCoordRep>
::TensorConnectedImageFunction()
{
  m_Radius.Fill(1);
}


/**
 *
 */
template <class TInputImage, class TCoordRep>
void
TensorConnectedImageFunction<TInputImage,TCoordRep>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Radius: " << m_Radius << std::endl;
}


template <class TInputImage, class TCoordRep>
void
TensorConnectedImageFunction<TInputImage, TCoordRep>
::ThresholdBetween(double lower, double upper) {
  if (m_Lower != lower || m_Upper != upper) {
    m_Lower = lower;
    m_Upper = upper;
    this->Modified();
  }
}

/**
 *
 */
template <class TInputImage, class TCoordRep>
bool
TensorConnectedImageFunction<TInputImage,TCoordRep>
::EvaluateAtIndex(const IndexType& index) const
{
  
  if( !this->GetInputImage() )
    {
    return ( false );
    }
  
  if ( !this->IsInsideBuffer( index ) )
    {
    return ( false );
    }

  // Create an N-d neighborhood kernel, using a zeroflux boundary condition
  ConstNeighborhoodIterator<InputImageType>
    it(m_Radius, this->GetInputImage(), this->GetInputImage()->GetBufferedRegion());

  // Set the iterator at the desired location
  it.SetLocation(index);
  PixelType curValue = this->GetInputImage()->GetPixel(index);
  TensorType curTensor = static_cast<TensorType>(curValue);

  TensorType::EigenValuesArrayType curEigenValues;
  TensorType::EigenVectorsMatrixType curEigenVectors;

  curTensor.ComputeEigenAnalysis(curEigenValues, curEigenVectors);


  // Walk the neighborhood
  bool allInside = true;
  double lower = this->GetLower();
  double upper = this->GetUpper();

  PixelType nbrValue;
  const unsigned int size = it.Size();
  for (unsigned int i = 0; i < size; ++i) {
    nbrValue = it.GetPixel(i);

    TensorType nbrTensor = static_cast<TensorType>(nbrValue);
    TensorType::EigenValuesArrayType nbrEigenValues;
    TensorType::EigenVectorsMatrixType nbrEigenVectors;
    nbrTensor.ComputeEigenAnalysis(nbrEigenValues, nbrEigenVectors);

    // stored in ascending order of eigen values
    double innerProduct = 0.0;
    for (int j = 0; j < 3; j++) {
      innerProduct += (curEigenVectors[j][2] * nbrEigenVectors[j][2]);
    }

    std::cout << innerProduct << std::endl;

    if (lower > innerProduct || innerProduct > upper) {
      allInside = false;
      break;
    }
  }

  return ( allInside );
}


} // end namespace itk

#endif
