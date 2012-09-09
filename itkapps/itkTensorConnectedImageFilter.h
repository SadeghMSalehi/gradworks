/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorConnectedImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2010/09/23 21:08:25 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorConnectedImageFilter_h
#define __itkTensorConnectedImageFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"

namespace itk {

/** \class TensorConnectedImageFilter
 * \brief Label pixels that are connected to a seed and lie within a neighborhood
 * 
 * TensorConnectedImageFilter labels pixels with ReplaceValue that
 * are connected to an initial Seed AND whose neighbors all lie within a
 * Lower and Upper threshold range.
 *
 * \ingroup RegionGrowingSegmentation 
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT TensorConnectedImageFilter:
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef TensorConnectedImageFilter             Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods).  */
  itkTypeMacro(TensorConnectedImageFilter,
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
  
  void PrintSelf ( std::ostream& os, Indent indent ) const;

  /** Clear the seeds */
  void ClearSeeds();

  /** Set seed point. */
  void SetSeed(const IndexType & seed);

  /** Add a seed point */
  void AddSeed ( const IndexType & seed );

  /** Set/Get the lower threshold. The default is 0. */
  itkSetMacro(Lower, double);
  itkGetConstMacro(Lower, double);

  /** Set/Get the upper threshold. The default is the largest possible
   *  value for the InputPixelType. */
  itkSetMacro(Upper, double);
  itkGetConstMacro(Upper, double);
  
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
  TensorConnectedImageFilter();
  ~TensorConnectedImageFilter(){};
  std::vector<IndexType> m_Seeds;
  double    m_Lower;
  double    m_Upper;
  OutputImagePixelType   m_ReplaceValue;
  InputImageSizeType     m_Radius;

  
  // Override since the filter needs all the data for the algorithm
  void GenerateInputRequestedRegion();

  // Override since the filter produces the entire dataset
  void EnlargeOutputRequestedRegion(DataObject *output);
  void GenerateData();
  
private:
  TensorConnectedImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorConnectedImageFilter.txx"
#endif

#endif
