/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorConnectedImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2011/03/20 22:53:00 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorConnectedImageFunction_h
#define __itkVectorConnectedImageFunction_h

#include "itkImageFunction.h"
#include <vector>
#include <queue>

namespace itk
{


/**
 * \class VectorConnectedImageFunction
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
class ITK_EXPORT VectorConnectedImageFunction :
  public ImageFunction< TInputImage, bool, TCoordRep >
{
public:
  /** Standard class typedefs. */
  typedef VectorConnectedImageFunction                        Self;
  typedef ImageFunction<TInputImage,bool,TCoordRep>           Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorConnectedImageFunction, ImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** InputImageType typedef support. */
  typedef TInputImage InputImageType;

  /** PixelType typedef support. */
  typedef typename TInputImage::PixelType PixelType; 

	typedef Image<unsigned int, InputImageType::ImageDimension> VisitingMapType;


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

  /** Set the ratio of the neighborhood direction used in computation. */
  itkSetMacro(Ratio, double);

  /** Get the ratio value */
  itkGetConstReferenceMacro(Ratio, double);

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
  
  /** Evalulate the function at specified index */
  virtual bool EvaluateNeighborAtIndex( const IndexType& nbrIndex, const IndexType& index ) const;

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

  virtual void SetVisitingMap(typename VisitingMapType::Pointer imageMapPointer) {
		m_VisitingMapPointer = imageMapPointer;
	}

  virtual void SetIndexStack(std::queue<typename InputImageType::IndexType>& indexStack) {
		m_IndexStack = indexStack;
	}
  
protected:
  VectorConnectedImageFunction();
  ~VectorConnectedImageFunction(){};

	std::queue<typename InputImageType::IndexType> m_IndexStack;
	typename VisitingMapType ::Pointer m_VisitingMapPointer;

  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  VectorConnectedImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  InputSizeType m_Radius;

  double m_Lower;
  double m_Upper;
	double m_Ratio;


};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_VectorConnectedImageFunction(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT VectorConnectedImageFunction< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef VectorConnectedImageFunction< ITK_TEMPLATE_2 x > \
                                         VectorConnectedImageFunction##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkVectorConnectedImageFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkVectorConnectedImageFunction.txx"
#endif

/*
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorConnectedImageFunction.txx"
#endif
*/

#endif
