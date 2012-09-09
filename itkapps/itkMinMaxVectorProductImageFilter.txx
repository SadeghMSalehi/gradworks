/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMinMaxVectorProductImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2010/09/23 21:08:25 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMinMaxVectorProductImageFilter_txx
#define __itkMinMaxVectorProductImageFilter_txx


// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#ifdef ITK_USE_CONSOLIDATED_MORPHOLOGY
#include "itkOptMinMaxVectorProductImageFilter.txx"
#else

#include "itkMinMaxVectorProductImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"

#include "iostream"
#include "iomanip"

namespace itk
{

template <class TInputImage, class TOutputImage>
MinMaxVectorProductImageFilter<TInputImage, TOutputImage>
::MinMaxVectorProductImageFilter()
{
  m_Radius.Fill(1);
}

template <class TInputImage, class TOutputImage>
void 
MinMaxVectorProductImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // get pointers to the input and output
  typename Superclass::InputImagePointer inputPtr = 
    const_cast< TInputImage * >( this->GetInput() );
  typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
  
  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius( m_Radius );

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()) )
    {
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    return;
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    
    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
    }
}


template< class TInputImage, class TOutputImage>
void
MinMaxVectorProductImageFilter< TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId)
{
  unsigned int i;
  ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;

  ConstNeighborhoodIterator<InputImageType> bit;
  ImageRegionConstIterator<InputImageType> inIt;
  ImageRegionIterator<OutputImageType> it;
  
  // Allocate output
  typename OutputImageType::Pointer output = this->GetOutput();
  typename  InputImageType::ConstPointer input  = this->GetInput();
  
  // Find the data-set boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> bC;
  faceList = bC(input, outputRegionForThread, m_Radius);

  std::cout << std::setbase(10);

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;

  // support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  
  InputRealType sum;

  // Process each of the boundary faces.  These are N-d regions which border
  // the edge of the buffer.
  for (fit=faceList.begin(); fit != faceList.end(); ++fit)
	{ 
    bit = ConstNeighborhoodIterator<InputImageType>(m_Radius,
                                                    input, *fit);
    unsigned int neighborhoodSize = bit.Size();

    inIt = ImageRegionConstIterator<InputImageType>(input, *fit);
    it = ImageRegionIterator<OutputImageType>(output, *fit);

    bit.OverrideBoundaryCondition(&nbc);
    bit.GoToBegin();

    while ( ! bit.IsAtEnd() )
		{
			bool nanStop = false;
			InputPixelType curValue = inIt.Value();
			for (int j = 0; j < 3; j++) {
			if (curValue[j] != curValue[j]) {
					nanStop = true;
				}
			}
			if (!nanStop) {
				sum = NumericTraits<InputRealType>::Zero;

				double min = 100, max = -100;
				for (i = 0; i < neighborhoodSize; ++i)
				{
					//sum += static_cast<InputRealType>( bit.GetPixel(i) );
					InputPixelType nbrValue = bit.GetPixel(i);
					for (int j = 0; j < 3; j++) {
						if (nbrValue[j] != nbrValue[j]) {
							nanStop = true;
						}
					}
					if (!nanStop) {
						double innerProduct = 0.0;
						for (int j = 0; j < 3; j++) {
							innerProduct += curValue[j] * nbrValue[j];
						}
						min = (min > innerProduct) ? innerProduct : min;
						max = (max < innerProduct) ? innerProduct : max;
					}
				}
				it.Set( static_cast<OutputPixelType>(min * 255) );
      } else {
				// get the mean value
				it.Set( static_cast<OutputPixelType>(0) );
			}
			
			++bit;
			++it;
			++inIt;
			progress.CompletedPixel();
		}
	}
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutput>
void
MinMaxVectorProductImageFilter<TInputImage, TOutput>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Radius: " << m_Radius << std::endl;

}

} // end namespace itk

#endif

#endif
