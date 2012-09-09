/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorConnectedImageFunction.txx,v $
  Language:  C++
  Date:      $Date: 2011/03/20 22:53:00 $
  Version:   $Revision: 1.10 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorConnectedImageFunction_txx
#define __itkVectorConnectedImageFunction_txx

#include "itkVectorConnectedImageFunction.h"
#include "itkNumericTraits.h"
#include "itkConstNeighborhoodIterator.h"
#include <iostream>
#include <math.h>

#include "itkVector.h"

using namespace std;

typedef itk::Vector<double> VectorType;

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TCoordRep>
VectorConnectedImageFunction<TInputImage,TCoordRep>
::VectorConnectedImageFunction()
{
  m_Radius.Fill(1);
	m_Ratio = 1.0;
}


/**
 *
 */
template <class TInputImage, class TCoordRep>
void
VectorConnectedImageFunction<TInputImage,TCoordRep>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Radius: " << m_Radius << std::endl;
}


template <class TInputImage, class TCoordRep>
void
VectorConnectedImageFunction<TInputImage, TCoordRep>
::ThresholdBetween(double lower, double upper) {
  if (m_Lower != lower || m_Upper != upper) {
    m_Lower = lower;
    m_Upper = upper;

    cout << "Threshold: " << m_Lower << " - " << m_Upper << endl;
    this->Modified();
  }
}

template <class TInputImage, class TCoordRep>
bool
VectorConnectedImageFunction<TInputImage,TCoordRep>
::EvaluateAtIndex(const IndexType& index) const 
{
	std::cout << "Don't come here (evaluatAtIndex) " << std::endl;
	return false;
}

/**
 *
 */
template <class TInputImage, class TCoordRep>
bool
VectorConnectedImageFunction<TInputImage,TCoordRep>
::EvaluateNeighborAtIndex(const IndexType& nbrIndex, const IndexType& index) const
{
  
  if( !this->GetInputImage() )
    {
    return ( false );
    }
  
  if ( !this->IsInsideBuffer( index ) )
    {
    return ( false );
    }

	VectorType value = static_cast<VectorType>(this->GetInputImage()->GetPixel(index));
	VectorType nbrValue = static_cast<VectorType>(this->GetInputImage()->GetPixel(nbrIndex));

	// return false for null vectors
	if (value[0] == 0 && value[1] == 0 && value[2] == 0) {
    cout << "Null value vector" << endl;
		return false;
	}
	if (nbrValue[0] == 0 && nbrValue[1] == 0 && nbrValue[2] == 0) {
    cout << "Null neighbor value vector" << endl;
		return false;
	}

	double innerProduct = 0, magValue = 0, magNbrValue = 0;
  for (int i = 0; i < 3; i++) {
    magValue += value[i] * value[i];
    magNbrValue += nbrValue[i] * nbrValue[i];
  }

  magValue = sqrt(magValue);
  magNbrValue = sqrt(magNbrValue);

  if (isnan(magNbrValue) || isnan(magValue)) {
    cout << "Null magnitude value or neighbor vector" << endl;
    return false;
  }

	for (int i = 0; i < 3; i++) {
      value[i] = value[i] / magValue;
      nbrValue[i] = nbrValue[i] / magNbrValue;
			innerProduct += (value[i] * nbrValue[i]) ;	
	}

	if (isnan(innerProduct)) {
    cout << "Null inner product " << endl;
		return false;
	}

  innerProduct = abs(innerProduct);
  double theta = acos(innerProduct);
  /* 
  cout << innerProduct << " (" << value[0] << "," << value[1] << "," << value[2] << "), ";
  cout << " (" << nbrValue[0] << "," << nbrValue[1] << "," << nbrValue[2] << ") ";
  cout << " / " << magValue << "," << magNbrValue << endl;
  */

	if (innerProduct >= this->GetLower() && innerProduct <= this->GetUpper()) {
	//if (theta >= this->GetLower() && theta <= this->GetUpper()) {
    //cout << " confirmed" << endl;
    cout << "@, " << innerProduct << ", " << nbrValue[0] << "," << nbrValue[1] << "," << nbrValue[2] << endl;
		return true;
	} else {
  }
	return false;

	/*
  // Create an N-d neighborhood kernel, using a zeroflux boundary condition
  ConstNeighborhoodIterator<InputImageType>
    it(m_Radius, this->GetInputImage(), this->GetInputImage()->GetBufferedRegion());

  // Set the iterator at the desired location
  it.SetLocation(index);
  PixelType curValue = this->GetInputImage()->GetPixel(index);
  VectorType curVector = static_cast<VectorType>(curValue);

  for (int i = 0; i < 3; i++) {
    if (isnan(curVector[i])) {
      return false;
    }
  }

  // Walk the neighborhood
  bool allInside = true;
  double lower = this->GetLower();
  double upper = this->GetUpper();

  PixelType nbrValue;
  const unsigned int size = it.Size();
	int countSimilarDirection = 0;

  double maxInnerProduct = -1;
  for (unsigned int i = 0; i < size; ++i) {
    nbrValue = it.GetPixel(i);

    VectorType nbrVector = static_cast<VectorType>(nbrValue);

    // stored in ascending order of eigen values
    double innerProduct = 0.0;
    for (int j = 0; j < 3; j++) {
      innerProduct += (curVector[j] * nbrVector[j]);
    }

    if (isnan(innerProduct)) {
      continue;
    }

    if (maxInnerProduct < innerProduct) {
      maxInnerProduct = innerProduct;
    }

    if (lower <= innerProduct && innerProduct <= upper) {
			countSimilarDirection ++;
		}
  }

  double nbrRatio = double(countSimilarDirection) / double(size);
	allInside = (nbrRatio >= m_Ratio);

//  cout << countSimilarDirection << ", " << lower << ", " << upper << ", " << nbrRatio << endl;

  return ( allInside );
	*/
}


} // end namespace itk

#endif
