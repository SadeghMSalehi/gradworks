/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorConnectedImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2010/09/23 21:08:25 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorConnectedImageFilter_txx
#define __itkTensorConnectedImageFilter_txx

#include "itkTensorConnectedImageFilter.h"
//##include "itkNeighborhoodBinaryThresholdImageFunction.h"
#include "itkTensorConnectedImageFunction.h"
#include "itkFloodFilledImageFunctionConditionalIterator.h"
#include "itkProgressReporter.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage>
TensorConnectedImageFilter<TInputImage, TOutputImage>
::TensorConnectedImageFilter()
{
  m_Lower = -1;
  m_Upper = 1;
  m_ReplaceValue = NumericTraits<OutputImagePixelType>::One;
  m_Radius.Fill(1);
}

template <class TInputImage, class TOutputImage>
void
TensorConnectedImageFilter<TInputImage, TOutputImage>
::ClearSeeds()
{
  if( this->m_Seeds.size() > 0 )
    {
    this->m_Seeds.clear();
    this->Modified();
    }
}

template <class TInputImage, class TOutputImage>
void
TensorConnectedImageFilter<TInputImage, TOutputImage>
::SetSeed(const IndexType & seed)
{
  this->ClearSeeds();
  this->AddSeed ( seed );
}

template <class TInputImage, class TOutputImage>
void
TensorConnectedImageFilter<TInputImage, TOutputImage>
::AddSeed ( const IndexType & seed )
{
  this->m_Seeds.push_back ( seed );
  this->Modified();
}

/**
 * Standard PrintSelf method.
 */
template <class TInputImage, class TOutputImage>
void
TensorConnectedImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Upper: "
     << static_cast<typename NumericTraits<InputImagePixelType>::PrintType>(m_Upper)
     << std::endl;
  os << indent << "Lower: "
     << static_cast<typename NumericTraits<InputImagePixelType>::PrintType>(m_Lower)
     << std::endl;
  os << indent << "ReplaceValue: "
     << static_cast<typename NumericTraits<OutputImagePixelType>::PrintType>(m_ReplaceValue)
     << std::endl;
  os << indent << "Radius: " << m_Radius << std::endl;
}

template <class TInputImage, class TOutputImage>
void 
TensorConnectedImageFilter<TInputImage,TOutputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if ( this->GetInput() )
    {
    InputImagePointer image =
      const_cast< InputImageType * >( this->GetInput() );
    image->SetRequestedRegionToLargestPossibleRegion();
    }
}

template <class TInputImage, class TOutputImage>
void 
TensorConnectedImageFilter<TInputImage,TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

template <class TInputImage, class TOutputImage>
void 
TensorConnectedImageFilter<TInputImage,TOutputImage>
::GenerateData()
{
  typename Superclass::InputImageConstPointer inputImage  = this->GetInput();
  typename Superclass::OutputImagePointer     outputImage = this->GetOutput();

  // Zero the output
  outputImage->SetBufferedRegion( outputImage->GetRequestedRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer ( NumericTraits<OutputImagePixelType>::Zero );
  
//  typedef NeighborhoodBinaryThresholdImageFunction<InputImageType> FunctionType;
  typedef TensorConnectedImageFunction<InputImageType> FunctionType;
  typedef FloodFilledImageFunctionConditionalIterator<OutputImageType, FunctionType> IteratorType;

  typename FunctionType::Pointer function = FunctionType::New();
  function->SetInputImage ( inputImage );
  function->ThresholdBetween ( m_Lower, m_Upper );
  function->SetRadius (m_Radius);
  IteratorType it = IteratorType ( outputImage, function, m_Seeds );

  ProgressReporter progress( this, 0,
                             outputImage->GetRequestedRegion().GetNumberOfPixels());
  while( !it.IsAtEnd())
    {
    it.Set(m_ReplaceValue);
    ++it;
    progress.CompletedPixel();
    }
}


} // end namespace itk

#endif
