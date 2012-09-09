/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorConnectedImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2010/10/24 18:44:12 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorConnectedImageFilter_h
#define __itkVectorConnectedImageFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include <vector>

namespace itk {

/** \class VectorConnectedImageFilter
 * \brief Label pixels that are connected to a seed and lie within a neighborhood
 * 
 * VectorConnectedImageFilter labels pixels with ReplaceValue that
 * are connected to an initial Seed AND whose neighbors all lie within a
 * Lower and Upper threshold range.
 *
 * \ingroup RegionGrowingSegmentation 
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT VectorConnectedImageFilter:
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef VectorConnectedImageFilter             Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods).  */
  itkTypeMacro(VectorConnectedImageFilter,
               ImageToImageFilter);

  typedef TInputImage                         InputImageType;
  typedef typename InputImageType::Pointer    InputImagePointer;
  typedef typename InputImageType::RegionType InputImageRegionType; 
  typedef typename InputImageType::PixelType  InputImagePixelType; 
  typedef typename InputImageType::IndexType  IndexType;
  typedef typename InputImageType::SizeType   InputImageSizeType;
  
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType; 
  typedef typename OutputImageType::PixelType  OutputImagePixelType; 

  typedef std::vector<IndexType> IndexVectorType;
  typedef std::vector<IndexVectorType> IndexVectorCollectionType;
  
  void PrintSelf ( std::ostream& os, Indent indent ) const;

  /** Clear the seeds */
  void ClearSeeds();

  /** Set seed point. */
  void SetSeed(const IndexType & seed);

  /** Add a seed point */
  void AddSeed ( const IndexType & seed );

  void BeginNextSeed() {
    IndexVectorType newSeedVector;
    // m_SeedCollection.insert(m_SeedCollection.begin(), newSeedVector);
    m_SeedCollection.push_back(newSeedVector);
  };

  /** Set/Get the lower threshold. The default is 0. */
  void AddLower(double lower) {
    m_LowerValues.push_back(lower);
  }

  void AddUpper(double upper) {
    m_UpperValues.push_back(upper);
  }

  void AddRatio(double ratio) {
    m_RatioValues.push_back(ratio);
  }

  /** Set/Get value to replace thresholded pixels. Pixels that lie *
   *  within Lower and Upper (inclusive) will be replaced with this
   *  value. The default is 1. */
  itkSetMacro(ReplaceValue, OutputImagePixelType);
  itkGetConstMacro(ReplaceValue, OutputImagePixelType);

  /** Set the radius of the neighborhood used for a mask. */
  itkSetMacro(Radius, InputImageSizeType);

  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Radius, InputImageSizeType);

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputEqualityComparableCheck,
    (Concept::EqualityComparable<InputImagePixelType>));
  itkConceptMacro(OutputEqualityComparableCheck,
    (Concept::EqualityComparable<OutputImagePixelType>));
  itkConceptMacro(SameDimensionCheck,
    (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
  itkConceptMacro(InputOStreamWritableCheck,
    (Concept::OStreamWritable<InputImagePixelType>));
  itkConceptMacro(OutputOStreamWritableCheck,
    (Concept::OStreamWritable<OutputImagePixelType>));
  /** End concept checking */
#endif

protected:
  VectorConnectedImageFilter();
  ~VectorConnectedImageFilter(){};
  IndexVectorCollectionType m_SeedCollection;

  std::vector<double> m_LowerValues;
  std::vector<double> m_UpperValues;
  std::vector<double> m_RatioValues;

  OutputImagePixelType   m_ReplaceValue;
  InputImageSizeType     m_Radius;

  
  // Override since the filter needs all the data for the algorithm
  void GenerateInputRequestedRegion();

  // Override since the filter produces the entire dataset
  void EnlargeOutputRequestedRegion(DataObject *output);
  void GenerateData();
  
private:
  VectorConnectedImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorConnectedImageFilter.txx"
#endif

#endif
