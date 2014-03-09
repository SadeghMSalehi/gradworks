//
//  itkFastMeanSquaredMetric.h
//  itkcmds
//
//  Created by Joohwi Lee on 9/15/12.
//
//

#ifndef itkcmds_itkFastMeanSquaredMetric_h
#define itkcmds_itkFastMeanSquaredMetric_h

#include "itkImageIO.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkArray.h"
#include "itkSingleValuedCostFunction.h"

using namespace std;
using namespace itk;

template <class TImage, class TLabel, class TTransform, class TInterpolator>
class FastMeanSquaredMetric : public itk::SingleValuedCostFunction {
public:
    /** Standard class typedefs. */
    typedef FastMeanSquaredMetric   Self;
    typedef SingleValuedCostFunction               Superclass;
    typedef SmartPointer< Self >       Pointer;
    typedef SmartPointer< const Self > ConstPointer;

    typedef typename TImage::Pointer ImagePointer;
    typedef typename TLabel::Pointer LabelPointer;
    typedef typename TTransform::Pointer TransformPointer;
    typedef typename TInterpolator::Pointer InterpolatorPointer;
    typedef ImageRegionConstIteratorWithIndex<TLabel> IteratorType;

    typedef vector<typename TImage::IndexType> IndexArray;
    typedef vector<typename TImage::PointType> PointArray;
    typedef vector<typename TImage::PixelType> PixelArray;
    typedef itkMathCode<TImage, TTransform> MathCodeType;
    typedef typename itkMathCode<TImage, TTransform>::Matrix Mat4;
    typedef typename itkMathCode<TImage, TTransform>::Vector Vec4;


    /** Run-time type information (and related methods). */
    itkTypeMacro(FastMeanSquaredMetric, SingleValuedCostFunction);

    itkNewMacro(FastMeanSquaredMetric);

    /**  MeasureType typedef.
     *  It defines a type used to return the cost function value. */
    typedef double MeasureType;

    /**  ParametersType typedef.
     *  It defines a position in the optimization search space. */
    typedef Superclass::ParametersType      ParametersType;
    typedef Superclass::ParametersValueType ParametersValueType;

    /** DerivativeType typedef.
     *  It defines a type used to return the cost function derivative.  */
    typedef Array< ParametersValueType > DerivativeType;

    /** This method returns the value of the cost function corresponding
     * to the specified parameters.    */
    virtual MeasureType GetValue(const ParametersType & parameters) const {
        cout << "Parameters: " << parameters << endl;

        Timer timer;
        timer.start();

        TransformPointer transform = TTransform::New();
        transform->SetParameters(parameters);

        MathCodeType mathCode;
        Mat4 imageTransform;
        mathCode.createIndex2IndexTransformMatrix(m_MovingImage, transform, m_FixedImage, imageTransform);

        MeasureType sum = 0;

        InterpolatorPointer fixedInterpolator = TInterpolator::New();
        fixedInterpolator->SetInputImage(m_FixedImage);

        typename TImage::RegionType targetRegion = m_FixedImage->GetBufferedRegion();

        itkcmds::itkImageIO<TImage> imageIO;
        ImagePointer resampledImage = imageIO.NewImageT(m_FixedImage);

        for (unsigned int i = 0; i < _indexArray.size(); i++) {
            MathCode::Vec4 input(_indexArray[i][0], _indexArray[i][1], _indexArray[i][2], 1);
            MathCode::Vec4 output;
            MathCode::mult(imageTransform, input, output);
            typename TInterpolator::ContinuousIndexType targetIndex;
            for (int j = 0; j < 3; j++) {
                targetIndex[j] = output[j];
            }

            if (targetRegion.IsInside(targetIndex)) {
                double valueAtFixed = fixedInterpolator->EvaluateAtContinuousIndex(targetIndex);
                double valueAtMoving = _pixelArray[i];
                double diff = valueAtFixed - valueAtMoving;
                sum += (diff * diff);
            }
        }
        return sum;
    }

    /** This method returns the derivative of the cost function corresponding
     * to the specified parameters.   */
    virtual void GetDerivative(const ParametersType & parameters,
                               DerivativeType & derivative) const {
        return;
    }

    /** This method returns the value and derivative of the cost function corresponding
     * to the specified parameters    */
    virtual void GetValueAndDerivative(const ParametersType & parameters,
                                       MeasureType & value,
                                       DerivativeType & derivative) const
    {
        value = this->GetValue(parameters);
        this->GetDerivative(parameters, derivative);
    }

    /** Return the number of parameters required to compute
     *  this cost function.
     *  This method MUST be overloaded by derived classes. */
    virtual unsigned int GetNumberOfParameters(void) const {
        return m_Transform->GetNumberOfParameters();
    }

    itkSetMacro(FixedImage, ImagePointer);
    itkGetMacro(FixedImage, ImagePointer);

    itkSetMacro(MovingImage, ImagePointer);
    itkGetMacro(MovingImage, ImagePointer);

    itkSetMacro(FixedImageMask, LabelPointer);
    itkGetMacro(FixedImageMask, LabelPointer);

    itkSetMacro(MovingImageMask, LabelPointer);
    itkGetMacro(MovingImageMask, LabelPointer);

    itkSetMacro(Transform, TransformPointer);
    itkGetMacro(Transform, TransformPointer);

    itkSetMacro(Interpolator, InterpolatorPointer);
    itkGetMacro(Interpolator, InterpolatorPointer);

    void Initialize() {
        Timer timer;
        timer.start();
        int _label = 2;
        if (m_MovingImageMask.GetPointer() == NULL) {
            cout << "Fixed image is not initialized" << endl;
            return;
        }
        IteratorType iter(m_MovingImageMask, m_MovingImageMask->GetRequestedRegion());
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            if (iter.Get() != _label) {
                continue;
            }
            _indexArray.push_back(iter.GetIndex());
            _pixelArray.push_back(iter.Get());
        }

        timer.stop();
        cout << "Number of index: " << _indexArray.size() << " (" << timer.getElapsedTimeInMilliSec() << " ms)" << endl;
    }


protected:
    FastMeanSquaredMetric() {}
    virtual ~FastMeanSquaredMetric() {}
private:
    FastMeanSquaredMetric(const Self &); //purposely not implemented
    void operator=(const Self &);           //purposely not implemented

    typename TImage::Pointer m_FixedImage;
    typename TLabel::Pointer m_FixedImageMask;
    typename TImage::Pointer m_MovingImage;
    typename TLabel::Pointer m_MovingImageMask;
    TransformPointer m_Transform;
    InterpolatorPointer m_Interpolator;
    
    
    IndexArray _indexArray;
    PointArray _pointArray;
    PixelArray _pixelArray;
};


#endif
