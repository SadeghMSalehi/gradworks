/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTensorConnectedImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2010/09/23 21:08:25 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTensorConnectedImageFunction_h
#define __itkTensorConnectedImageFunction_h

#include "itkImageFunction.h"

namespace itk
{


/**
 * \class TensorConnectedImageFunction
 * \brief Determine whether all the pixels in the specified neighborhood meet a threshold criteria
 *
 * Determine whether all the pixels in the specified neighborhood meet
 * a threshold criteria.
 *
 * If called with a ContinuousIndex or Point, the calculation is performed
 * at the nearest neighbor.
 *
 * This class is templated over the input image type and the coordinate
 * representation type (e.g. float or double).
 *
 * \ingroup ImageFunctions
 */
template <class TInputImage, class TCoordRep = float >
class ITK_EXPORT TensorConnectedImageFunction :
  public ImageFunction< TInputImage, bool, TCoordRep >
{
public:
  /** Standard class typedefs. */
  typedef TensorConnectedImageFunction                        Self;
  typedef ImageFunction<TInputImage,bool,TCoordRep>           Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorConnectedImageFunction, ImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** InputImageType typedef support. */
  typedef TInputImage InputImageType;

  /** PixelType typedef support. */
  typedef typename TInputImage::PixelType PixelType; 

  /** Dimension of the underlying image. */
  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** OutputType typdef support. */
  typedef typename Superclass::OutputType OutputType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** SizeType of the input image */
  typedef typename InputImageType::SizeType InputSizeType;

  /** Set the radius of the neighborhood used in computation. */
  itkSetMacro(Radius, InputSizeType);

  /** Get the lower threshold value */
  itkGetConstReferenceMacro(Lower, double);

  /** Get the upper threshold value */
  itkGetConstReferenceMacro(Upper, double);

  /** Get the radius of the neighborhood used in computation */
  itkGetConstReferenceMacro(Radius, InputSizeType);

  /** Values that lie between lower and upper inclusive are inside */
  void ThresholdBetween(double lower, double upper);

  /** Evalulate the function at specified index */
  virtual bool EvaluateAtIndex( const IndexType& index ) const;
  
  /** Evaluate the function at non-integer positions */
  virtual bool Evaluate( const PointType& point ) const
    { 
    IndexType index;
    this->ConvertPointToNearestIndex( point, index );
    return this->EvaluateAtIndex( index ); 
    }

  virtual bool EvaluateAtContinuousIndex( 
    const ContinuousIndexType& cindex ) const
    { 
    IndexType index;
    this->ConvertContinuousIndexToNearestIndex( cindex, index );
    return this->EvaluateAtIndex( index ); 
    }
  
protected:
  TensorConnectedImageFunction();
  ~TensorConnectedImageFunction(){};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  TensorConnectedImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  InputSizeType m_Radius;

  double m_Lower;
  double m_Upper;

};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_TensorConnectedImageFunction(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT TensorConnectedImageFunction< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef TensorConnectedImageFunction< ITK_TEMPLATE_2 x > \
                                         TensorConnectedImageFunction##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkTensorConnectedImageFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkTensorConnectedImageFunction.txx"
#endif

/*
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorConnectedImageFunction.txx"
#endif
*/

#endif
