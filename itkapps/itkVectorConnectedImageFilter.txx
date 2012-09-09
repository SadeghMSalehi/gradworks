/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorConnectedImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2011/03/20 22:53:00 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorConnectedImageFilter_txx
#define __itkVectorConnectedImageFilter_txx

#include "itkVectorConnectedImageFilter.h"
//##include "itkNeighborhoodBinaryThresholdImageFunction.h"
#include "itkVectorConnectedImageFunction.h"
#include "mrRegionGrowingIterator.h"
#include "itkProgressReporter.h"

#include <vector>

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage>
VectorConnectedImageFilter<TInputImage, TOutputImage>
::VectorConnectedImageFilter()
{
  m_LowerValues.clear();
  m_UpperValues.clear();
  m_RatioValues.clear();
  m_ReplaceValue = NumericTraits<OutputImagePixelType>::One;
  m_Radius.Fill(1);
  BeginNextSeed();
}

template <class TInputImage, class TOutputImage>
void
VectorConnectedImageFilter<TInputImage, TOutputImage>
::ClearSeeds()
{
  if( this->m_SeedCollection.size() > 0 )
    {
    this->m_SeedCollection.clear();
    this->Modified();
    }

  m_LowerValues.clear();
  m_UpperValues.clear();
  m_RatioValues.clear();

  BeginNextSeed();
}

template <class TInputImage, class TOutputImage>
void
VectorConnectedImageFilter<TInputImage, TOutputImage>
::SetSeed(const IndexType & seed)
{
  this->ClearSeeds();
  this->AddSeed ( seed );
}

template <class TInputImage, class TOutputImage>
void
VectorConnectedImageFilter<TInputImage, TOutputImage>
::AddSeed ( const IndexType & seed )
{
  m_SeedCollection.back().push_back ( seed );
  this->Modified();
}

/**
 * Standard PrintSelf method.
 */
template <class TInputImage, class TOutputImage>
void
VectorConnectedImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
}

template <class TInputImage, class TOutputImage>
void 
VectorConnectedImageFilter<TInputImage,TOutputImage>
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
VectorConnectedImageFilter<TInputImage,TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

template <class TInputImage, class TOutputImage>
void 
VectorConnectedImageFilter<TInputImage,TOutputImage>
::GenerateData()
{
  typename Superclass::InputImageConstPointer inputImage  = this->GetInput();
  typename Superclass::OutputImagePointer     outputImage = this->GetOutput();

  // Zero the output
  outputImage->SetBufferedRegion( outputImage->GetRequestedRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer ( NumericTraits<OutputImagePixelType>::Zero );
  
//  typedef NeighborhoodBinaryThresholdImageFunction<InputImageType> FunctionType;
  typedef VectorConnectedImageFunction<InputImageType> FunctionType;
  typedef RegionGrowingIterator<OutputImageType, FunctionType> RegionGrowingIteratorType;

  typename FunctionType::Pointer function = FunctionType::New();
  function->SetInputImage ( inputImage );
  function->SetRadius (m_Radius);
  RegionGrowingIteratorType rgIter = RegionGrowingIteratorType ( outputImage, function );

  ProgressReporter progress( this, 0,
                             outputImage->GetRequestedRegion().GetNumberOfPixels());

  unsigned char visitingLabel = 0;
  typename IndexVectorCollectionType::iterator seedCollectionIterator = m_SeedCollection.begin();

  int seedCollectionIndex = 0;
  while (seedCollectionIterator != m_SeedCollection.end()) {
    seedCollectionIndex ++;
    IndexVectorType seeds = *seedCollectionIterator;
    typename IndexVectorType::iterator seedIterator = seeds.begin();

    function->ThresholdBetween ( m_LowerValues[seedCollectionIndex-1], m_UpperValues[seedCollectionIndex-1] );
    function->SetRatio(m_RatioValues[seedCollectionIndex-1]);

    rgIter.ClearSeeds();
    while (seedIterator != seeds.end()) {
      IndexType seedIndex = *seedIterator;
      rgIter.AddSeed(seedIndex);
      ++seedIterator;
      cout << seedCollectionIndex << "] " << seedIndex << endl;
    }

    //it.SetVisitedLabel(++visitingLabel);
    //it.SetIncludedLabel(++visitingLabel);

    rgIter.GoToNextBegin();

    typedef itk::Vector<double> VectorType;
    while(!rgIter.IsAtEnd())
      {
      rgIter.Set(m_ReplaceValue++);
      ++rgIter;
      IndexType idx = rgIter.GetIndex();
      VectorType val = static_cast<VectorType>(inputImage->GetPixel(idx));
      cout << "#," << val[0] << "," << val[1] << "," << val[2] << endl;
      progress.CompletedPixel();
      }

    ++seedCollectionIterator;
  }

  std::string filename = "temp.nrrd";
  rgIter.WriteVisitingMap(filename);
}


} // end namespace itk

#endif
