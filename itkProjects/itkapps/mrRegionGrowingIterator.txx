/*=========================================================================
 *

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: mrRegionGrowingIterator.txx,v $
  Language:  C++
  Date:      $Date: 2011/03/20 22:53:00 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __mrRegionGrowingIterator_txx
#define __mrRegionGrowingIterator_txx

#include "mrRegionGrowingIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk
{
template<class TImage, class TFunction>
RegionGrowingIterator<TImage, TFunction>
::RegionGrowingIterator(const ImageType *imagePtr,
                                              FunctionType *fnPtr,
                                              IndexType startIndex)
{
  this->m_Image = imagePtr;
  m_Function = fnPtr;
  m_StartIndices.push_back ( startIndex );

  m_VisitedLabel = 1;
  m_IncludedLabel = 30;

  // Set up the temporary image
  this->InitializeIterator();
}

template<class TImage, class TFunction>
RegionGrowingIterator<TImage, TFunction>
::RegionGrowingIterator(const ImageType *imagePtr,
                                              FunctionType *fnPtr,
                                              std::vector<IndexType>& startIndex)
{
  this->m_Image = imagePtr;
  m_Function = fnPtr;
  unsigned int i;
  for (i = 0; i < startIndex.size(); i++ )
    {
    m_StartIndices.push_back ( startIndex[i] );
    }

  m_VisitedLabel = 1;
  m_IncludedLabel = 30;

  // Set up the temporary image
  this->InitializeIterator();
}

template<class TImage, class TFunction>
RegionGrowingIterator<TImage, TFunction>
::RegionGrowingIterator(const ImageType *imagePtr,
                                              FunctionType *fnPtr)
{
  this->m_Image = imagePtr;
  m_Function = fnPtr;

  m_VisitedLabel = 1;
  m_IncludedLabel = 30;

  // Set up the temporary image
  this->InitializeIterator();
}

template<class TImage, class TFunction>
void
RegionGrowingIterator<TImage, TFunction>
::InitializeIterator()
{
  // Get the origin and spacing from the image in simple arrays
  m_ImageOrigin  = this->m_Image->GetOrigin();
  m_ImageSpacing = this->m_Image->GetSpacing();
  m_ImageRegion  = this->m_Image->GetBufferedRegion();

  // Build a temporary image of chars for use in the flood algorithm
  m_TemporaryPointer = TTempImage::New();
  typename TTempImage::RegionType tempRegion = this->m_Image->GetBufferedRegion();

  m_TemporaryPointer->SetLargestPossibleRegion( tempRegion );
  m_TemporaryPointer->SetBufferedRegion( tempRegion );
  m_TemporaryPointer->SetRequestedRegion( tempRegion );
  m_TemporaryPointer->Allocate();
  m_TemporaryPointer->FillBuffer(NumericTraits<ITK_TYPENAME TTempImage::PixelType>::Zero);

  // Initialize the queue by adding the start index assuming one of
  // the m_StartIndices is "inside" This might not be true, in which
  // case it's up to the programmer to specify a correct starting
  // position later (using FindSeedPixel).  Must make sure that the
  // seed is inside the buffer before touching pixels.
  this->m_IsAtEnd = true;
  for ( unsigned int i = 0; i < m_StartIndices.size(); i++ )
    {
    if ( m_ImageRegion.IsInside ( m_StartIndices[i] ) )
      {
      m_IndexStack.push(m_StartIndices[i]);
      m_TemporaryPointer->SetPixel(m_StartIndices[i], markBitMap(0, m_IncludedLabel));
      this->m_IsAtEnd = false;
      }
    }

  this->m_Function->SetVisitingMap(m_TemporaryPointer);
  this->m_Function->SetIndexStack(m_IndexStack);
}

template<class TImage, class TFunction>
void
RegionGrowingIterator<TImage, TFunction>
::FindSeedPixel()
{
  // Create an iterator that will walk the input image
  typedef typename itk::ImageRegionConstIterator<TImage> IRIType;
  IRIType it = IRIType(this->m_Image, this->m_Image->GetBufferedRegion() );
  
  // Now we search the input image for the first pixel which is inside
  // the function of interest
  m_StartIndices.clear();
  for ( it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
    if( this->IsPixelIncluded( it.GetIndex() ) )
      {
      m_StartIndices.push_back ( it.GetIndex() );

      // We need to reset the "beginning" now that we have a real seed
      this->GoToBegin();

      return;
      }
    }
}

template<class TImage, class TFunction>
void
RegionGrowingIterator<TImage, TFunction>
::FindSeedPixels()
{
  // Create an iterator that will walk the input image
  typedef typename itk::ImageRegionConstIterator<TImage> IRIType;
  IRIType it = IRIType(this->m_Image, this->m_Image->GetBufferedRegion() );
  
  // Now we search the input image for the first pixel which is inside
  // the function of interest
  m_StartIndices.clear();
  bool found = false;
  for ( it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
    if( this->IsPixelIncluded( it.GetIndex() ) )
      {
      m_StartIndices.push_back ( it.GetIndex() );
      found = true;
      }
    }
  if ( found )
    {
    // We need to reset the "beginning" now that we have a real seed
    this->GoToBegin();
    }
}

template<class TImage, class TFunction>
void
RegionGrowingIterator<TImage, TFunction>
::DoFloodStep()
{
  // The index in the front of the queue should always be
  // valid and be inside since this is what the iterator
  // uses in the Set/Get methods. This is ensured by the
  // GoToBegin() method.
 

  // Take the index in the front of the queue  
  const IndexType & topIndex = m_IndexStack.front();
  
  if (NDimensions > 3) {
    cout << "Dimensions more than 3 are not supported " << endl;
    return;
  }

  typename TImage::SizeType radius;
  radius.Fill(1);
  ConstNeighborhoodIterator<TImage> it(radius, this->m_Image, m_ImageRegion);
  it.SetLocation(topIndex);
  m_TemporaryPointer->SetPixel(topIndex, markBitMap(0, m_IncludedLabel));

  const unsigned int nbrSize = it.Size();

  for (unsigned int i = 0; i < nbrSize; i++) {
    IndexType tempIndex = it.GetIndex(i);
    if (tempIndex == topIndex) {
      continue;
    }

    // neighbor check
    unsigned int nbrBitMap = m_TemporaryPointer->GetPixel(tempIndex);
    IndexType idxOffset;
    for (int i = 0; i < 3; i++) {
      idxOffset[i] = topIndex[i] - tempIndex[i] + 1;
    }

    // offset to neighbor order
    int nbrOrder = idxOffset[0] + idxOffset[1] * 3 + idxOffset[2] * 9; 
    if (!checkBitMap(nbrBitMap, nbrOrder) && !checkBitMap(nbrBitMap, m_IncludedLabel)) {
      nbrBitMap = markBitMap(nbrBitMap, nbrOrder);
      if (this->IsPixelIncluded(tempIndex, topIndex)) {
        m_IndexStack.push(tempIndex);
        nbrBitMap = markBitMap(nbrBitMap, m_IncludedLabel);
      }
      m_TemporaryPointer->SetPixel(tempIndex, nbrBitMap);
    }
  }

  /*
  // Iterate through all possible dimensions
  // NOTE: Replace this with a ShapeNeighborhoodIterator
  for(unsigned int i=0; i<NDimensions; i++)
    {
    // The j loop establishes either left or right neighbor (+-1)
    for(int j=-1; j<=1; j+=1)
      {
      IndexType tempIndex;

      // build the index of a neighbor
      for(unsigned int k=0; k<NDimensions; k++)
        {
        if( i!=k )
          {
          tempIndex.m_Index[k] = topIndex[k];
          }
        else
          {
          tempIndex.m_Index[k] = topIndex[k] + j;
          }
        } // end build the index of a neighbor


      // If this is a valid index and have not been tested,
      // then test it.
      if( m_ImageRegion.IsInside( tempIndex ) )
        {
        if( m_TemporaryPointer->GetPixel( tempIndex ) == 0 )
          {
          // if it is inside, push it into the queue  
          if(  this->IsPixelIncluded( tempIndex, topIndex ) )
            {
            m_IndexStack.push( tempIndex );
            m_TemporaryPointer->SetPixel( tempIndex, 2); 
            }
          else  // If the pixel is outside
            {
            // Mark the pixel as outside and remove it from the queue.
            m_TemporaryPointer->SetPixel( tempIndex, 1);
            }
          }
        }
      } // end left/right neighbor loop
    } // end check all neighbors
  */
  
  // Now that all the potential neighbors have been 
  // inserted we can get rid of the pixel in the front
  m_IndexStack.pop();
    
  if( m_IndexStack.empty() )
    {
    this->m_IsAtEnd = true;
    }


}

template<class TImage, class TFunction>
bool
RegionGrowingIterator<TImage, TFunction>
::IsPixelIncluded(const IndexType & nbrIndex, const IndexType & index) const
{
  //cout << nbrIndex << " <> " << index << endl;
  return this->m_Function->EvaluateNeighborAtIndex(nbrIndex, index);
}



} // end namespace itk

#endif
