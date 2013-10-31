//
//  piPatchCompare.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/13/13.
//
//

#ifndef __ParticleGuidedRegistration__piPatchCompare__
#define __ParticleGuidedRegistration__piPatchCompare__

#include <iostream>
#include "piImageDef.h"
#include "piParticleSystem.h"

#include <itkPDEDeformableRegistrationFilter.h>
#include <itkPDEDeformableRegistrationFunction.h>
#include <itkPoint.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImageFunction.h>
#include <itkCovariantVector.h>
#include <itkInterpolateImageFunction.h>
#include <itkDefaultConvertPixelTraits.h>
#include <itkEnableIf.h>
#include <itkIsSame.h>

#define PATCH_SIZE 5

namespace pi {

    typedef itk::Vector<ImageReal,PATCH_SIZE*PATCH_SIZE> LocalPatch;
    typedef itk::Image<LocalPatch,__Dim> PatchImage;
    typedef itk::ImageRegionIteratorWithIndex<PatchImage> PatchImageIterator;
    typedef std::vector<PatchImage::Pointer> PatchImageVector;

    class PatchCompare {
    public:
        PatchCompare() {}
        virtual ~PatchCompare() {}

        void setTargetRadius(int r);
        void setPatchRadius(int r);

        void setAtlasImages(PatchImageVector atlasImages);
        void setAtlasLabels(LabelImageVector atlasLabels);
        void setTargetImage(RealImage::Pointer target);
        void setTargetROI(LabelImage::Pointer target);
        void setParticleSystem(ParticleSystem* system);

        // return the estimated target label
        LabelImage::Pointer getTargetLabel();

        void estimateLabel(int searchRadius, int kNearest);

        PatchImage::Pointer buildPatchImage(RealImage::Pointer image);

        void buildGradientPatchImage(RealImage::Pointer image, const int radius, PatchImageVector& output);

        DisplacementFieldType::Pointer computeOpticalFlow(RealImage::Pointer fixed, RealImage::Pointer moving, double dt);

        DisplacementFieldType::Pointer computeDemonsFlow(RealImage::Pointer fixed, RealImage::Pointer moving, double dt);


        // run demon's algorithm on these patch image
        int performDemonsRegistration(RealImage::Pointer fixedImage, RealImage::Pointer movingImage, std::string outputImage);

        DisplacementFieldType::Pointer performDenseMapping(PatchImage::Pointer fixed, PatchImage::Pointer moving, PatchImage::RegionType activeRegion);

        /// main testing routines
        bool main(Options& parser, StringVector args);

    private:
        ParticleSystem* _system;
        RealImage::Pointer _targetImage;
        LabelImage::Pointer _targetROI;
        LabelImage::Pointer _targetLabel;
        PatchImageVector _atlasImages;
        LabelImageVector _atlasLabels;

        int _targetRadius;
        int _patchRadius;
    };


    /**
     * Compute label transfer from string arguments
     * \param args a list of string arguments
     * \param output filename for estimated label transfer
     */
    void transferLabelsWithPatch(StringVector& args, std::string output, int searchRadius = 3, int kNearest = 3);

#pragma mark -
    /**
     * \class CentralDifferenceImageFunction
     * \brief Calculate the derivative by central differencing.
     *
     * This class is templated over the input image type,
     * the coordinate representation type (e.g. float or double),
     * and the output derivative type.
     *
     * This class supports both scalar and vector pixel types
     * for the input image, including VectorImage types.
     *
     * For vector-pixel image types, the TOutputType template
     * parameter must be set to a vector of appropriate size, to
     * accomadate a result for each pixel component in each dimension.
     * The output is packed by pixel component, i.e.
     *
     *  [C0D0, C0D1, ..., C0DN, C1D0, ...]
     *
     * where C = pixel component, and D = image dimension.
     *
     * The output type can be, for example:
     *
     *     \code CovariantVector<double, numberOfPixelComponents * ImageDimension> \endcode
     *  or
     *     \code Matrix<double, numberOfPixelComponents, ImageDimension> \endcode
     *
     * Possible improvements:
     *
     * 1) speed performance:
     * The template-specialization of the Evaluate*() methods (needed
     * to support vector-pixel types) incur a performance penalty for the
     * scalar-pixel case, when compared with previous scalar-only
     * versions of the code. On MacOS (2.4GHz Core 2 Duo, gcc 4.2)
     * the penalty is 0.5-2%, depending on the method. To recover this loss,
     * the specialization of the methods would have to be done such that
     * a nested subroutine need not be called, ie the specialization is
     * performed on the Evaluate* methods directly. At the moment is seems
     * this can't be done without requiring a template parameter on the
     * methods.
     *
     * 2) the use of Neighborhood operators may improve efficiency.
     *
     * \ingroup ImageFunctions
     * \ingroup ITKImageFunction
     */
    template<
    class TInputImage,
    class TCoordRep = float,
    class TOutputType = itk::CovariantVector<double, TInputImage::ImageDimension >
    >
    class ITK_EXPORT CentralDifferenceImageFunction:
    public itk::ImageFunction< TInputImage,
    TOutputType,
    TCoordRep >
    {
    public:
        /** Dimension underlying input image. */
        itkStaticConstMacro(ImageDimension, unsigned int,
                            TInputImage::ImageDimension);

        /** Standard class typedefs. */
        typedef CentralDifferenceImageFunction   Self;
        typedef itk::ImageFunction< TInputImage,
        TOutputType,
        TCoordRep >       Superclass;
        typedef itk::SmartPointer< Self >             Pointer;
        typedef itk::SmartPointer< const Self >       ConstPointer;

        /** Run-time type information (and related methods). */
        itkTypeMacro(CentralDifferenceImageFunction, ImageFunction);

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** InputImageType typedef support. */
        typedef TInputImage InputImageType;

        /** InputPixelType typedef support */
        typedef typename InputImageType::PixelType InputPixelType;

        /** InputPixelConvert typedef support */
        typedef itk::DefaultConvertPixelTraits< InputPixelType > InputPixelConvertType;

        /** OutputType typdef support. */
        typedef typename Superclass::OutputType OutputType;

        /** Output convert typedef support */
        typedef itk::DefaultConvertPixelTraits<OutputType> OutputConvertType;

        /** Output value typedef support */
        typedef typename OutputConvertType::ComponentType OutputValueType;

        /** Scalar derivative typedef support */
        typedef itk::CovariantVector<OutputValueType, itkGetStaticConstMacro(ImageDimension) > ScalarDerivativeType;

        /** Index typedef support. */
        typedef typename Superclass::IndexType IndexType;

        /** ContinuousIndex typedef support. */
        typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

        /** Point typedef support. */
        typedef typename Superclass::PointType PointType;

        /** Spacing typedef support. */
        typedef typename TInputImage::SpacingType SpacingType;

        /** Interpolator typedef support. */
        typedef itk::InterpolateImageFunction< TInputImage, TCoordRep > InterpolatorType;
        typedef typename InterpolatorType::Pointer                 InterpolatorPointer;

        /** Set the input image.  This must be set by the user. */
        virtual void SetInputImage(const TInputImage *inputData);

        /** Set interpolator. The interpolator is used in the methods
         * \c Evaluate and \c EvaluateAtContinuousIndex. */
        virtual void SetInterpolator(InterpolatorType *interpolator);

        /** Get the interpolator. */
        itkGetConstObjectMacro( Interpolator, InterpolatorType );

        /** Evalulate the image derivative by central differencing at specified index.
         *
         *  No bounds checking is done.
         *  The point is assumed to lie within the image buffer.
         *
         *  If \c index lies on a boundary in a given dimension, 0 is returned for
         *  that dimension.
         *
         *  ImageFunction::IsInsideBuffer() can be used to check bounds before
         * calling the method. */
        virtual OutputType EvaluateAtIndex(const IndexType & index) const;

        /** Evalulate the image derivative by central differencing at non-integer
         *  point.
         *
         *  No bounds checking is done.
         *  The point is assumed to lie within the image buffer. If not, 0 is
         *  returned for the derivative without any error return, because of
         *  bounds-checking performed on the neighboring points.
         *
         *  If \c point lies on a boundary in a given dimension, 0 is returned for
         *  that dimension. Note that points are centered on the voxel.
         *
         *  ImageFunction::IsInsideBuffer() can be used to check bounds before
         * calling the method. */
        virtual OutputType Evaluate(const PointType & point) const;

        /** Evalulate the image derivative by central differencing at non-integer
         *  index.
         *
         *  No bounds checking is done.
         *  The point is assumed to lie within the image buffer.
         *
         *  If \c cindex lies on a boundary in a given dimension, 0 is returned for
         *  that dimension.
         *
         *  ImageFunction::IsInsideBuffer() can be used to check bounds before
         * calling the method. */
        virtual OutputType EvaluateAtContinuousIndex( const ContinuousIndexType & cindex) const;

        /** The UseImageDirection flag determines whether image derivatives are
         * computed with respect to the image grid or with respect to the physical
         * space. When this flag is ON the derivatives are computed with respect to
         * the coordinate system of physical space. The difference is whether we take
         * into account the image Direction or not.
         * For \c EvaluateAtIndex and \c EvaluateAtContinuousIndex, the flag ON will
         * take into account the image direction and will result in an extra matrix
         * multiplication compared to the amount of computation performed when the
         * flag is OFF.
         * For \c Evaluate, the opposite is true: the flag OFF will ignore the image
         * direction and will result in an extra matrix multiplication compared to the
         * amount of computation performed when the flag is ON.
         * The default value of this flag is On.
         */
        itkSetMacro(UseImageDirection, bool);
        itkGetConstMacro(UseImageDirection, bool);
        itkBooleanMacro(UseImageDirection);

    protected:
        CentralDifferenceImageFunction();
        ~CentralDifferenceImageFunction(){}
        void PrintSelf(std::ostream & os, itk::Indent indent) const;

    private:
        CentralDifferenceImageFunction(const Self &); //purposely not implemented
        void operator=(const Self &);                 //purposely not implemented


        /** Structure for specialization of Evaulate* methods on OutputType */
        template<typename T>
        struct OutputTypeSpecializationStructType
        {
            typedef T Type;
        };

        /** Specialized versions of EvaluteAtIndex() method to handle scalar or vector pixel types.*/
        template< class Type >
        inline void EvaluateAtIndexSpecialized( const IndexType & index, OutputType & derivative, OutputTypeSpecializationStructType<OutputType>) const;
        template< class Type >
        inline void EvaluateAtIndexSpecialized( const IndexType & index, OutputType & derivative, OutputTypeSpecializationStructType<Type>) const;

        /** Specialized versions of EvaluteAtContinuousIndex() method to handle scalar or vector pixel types.*/
        template< class Type >
        inline void EvaluateAtContinuousIndexSpecialized( const ContinuousIndexType & index, OutputType & derivative, OutputTypeSpecializationStructType<OutputType>) const;
        template< class Type >
        inline void EvaluateAtContinuousIndexSpecialized( const ContinuousIndexType & index, OutputType & derivative, OutputTypeSpecializationStructType<Type>) const;
        
        /** Specialized versions of Evalute() method to handle scalar or vector pixel types.*/
        // NOTE: for some unknown reason, making these methods inline (as those above are inlined) makes them run *slower*.
        template< class Type >
        void EvaluateSpecialized( const PointType & point, OutputType & derivative, OutputTypeSpecializationStructType<OutputType>) const;
        template< class Type >
        void EvaluateSpecialized( const PointType & point, OutputType & derivative, OutputTypeSpecializationStructType<Type>) const;
        
        // flag to take or not the image direction into account
        // when computing the derivatives.
        bool m_UseImageDirection;
        
        // interpolator
        InterpolatorPointer   m_Interpolator;
    };

#pragma mark -

    /**
     * \class DemonsRegistrationFunction
     *
     * This class encapsulate the PDE which drives the demons registration
     * algorithm. It is used by DemonsRegistrationFilter to compute the
     * output displacement field which will map a moving image onto a
     * a fixed image.
     *
     * Non-integer moving image values are obtained by using
     * interpolation. The default interpolator is of type
     * LinearInterpolateImageFunction. The user may set other
     * interpolators via method SetMovingImageInterpolator. Note that the input
     * interpolator must derive from baseclass InterpolateImageFunction.
     *
     * This class is templated over the fixed image type, moving image type,
     * and the displacement field type.
     *
     * \warning This filter assumes that the fixed image type, moving image type
     * and displacement field type all have the same number of dimensions.
     *
     * \sa DemonsRegistrationFilter
     * \ingroup FiniteDifferenceFunctions
     * \ingroup ITKPDEDeformableRegistration
     */
    template<class TFixedImage, class TMovingImage, class TDisplacementField>
    class DemonsRegistrationFunction : public itk::PDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>
    {
    public:
        /** Standard class typedefs. */
        typedef DemonsRegistrationFunction Self;
        typedef itk::PDEDeformableRegistrationFunction< TFixedImage,
        TMovingImage, TDisplacementField
        >                                      Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(DemonsRegistrationFunction,
                     itk::PDEDeformableRegistrationFunction);

        /** MovingImage image type. */
        typedef typename Superclass::MovingImageType    MovingImageType;
        typedef typename Superclass::MovingImagePointer MovingImagePointer;

        /** FixedImage image type. */
        typedef typename Superclass::FixedImageType    FixedImageType;
        typedef typename Superclass::FixedImagePointer FixedImagePointer;
        typedef typename FixedImageType::IndexType     IndexType;
        typedef typename FixedImageType::SizeType      SizeType;
        typedef typename FixedImageType::SpacingType   SpacingType;

        /** Deformation field type. */
        typedef typename Superclass::DisplacementFieldType        DisplacementFieldType;
        typedef typename Superclass::DisplacementFieldTypePointer DisplacementFieldTypePointer;

#ifdef ITKV3_COMPATIBILITY
        typedef typename Superclass::DeformationFieldType        DeformationFieldType;
        typedef typename Superclass::DeformationFieldTypePointer DeformationFieldTypePointer;
#endif

        /** Inherit some enums from the superclass. */
        itkStaticConstMacro(ImageDimension, unsigned
                            int, Superclass::ImageDimension);

        /** Inherit some enums from the superclass. */
        typedef typename Superclass::PixelType        PixelType;
        typedef typename Superclass::RadiusType       RadiusType;
        typedef typename Superclass::NeighborhoodType NeighborhoodType;
        typedef typename Superclass::FloatOffsetType  FloatOffsetType;
        typedef typename Superclass::TimeStepType     TimeStepType;

        /** Interpolator type. */
        typedef double                                                          CoordRepType;
        typedef itk::InterpolateImageFunction< MovingImageType, CoordRepType >       InterpolatorType;
        typedef typename InterpolatorType::Pointer                              InterpolatorPointer;
        typedef typename InterpolatorType::PointType                            PointType;
        typedef itk::LinearInterpolateImageFunction< MovingImageType, CoordRepType > DefaultInterpolatorType;

        /** Covariant vector type. */
        typedef itk::CovariantVector< double, itkGetStaticConstMacro(ImageDimension) > CovariantVectorType;

        /** Fixed image gradient calculator type. */
        typedef CentralDifferenceImageFunction< FixedImageType > GradientCalculatorType;
        typedef typename GradientCalculatorType::Pointer GradientCalculatorPointer;

        /** Moving image gradient calculator type. */
        typedef CentralDifferenceImageFunction<MovingImageType, CoordRepType>
        MovingImageGradientCalculatorType;
        typedef typename MovingImageGradientCalculatorType::Pointer
        MovingImageGradientCalculatorPointer;

        /** Set the moving image interpolator. */
        void SetMovingImageInterpolator(InterpolatorType *ptr)
        { m_MovingImageInterpolator = ptr; }

        /** Get the moving image interpolator. */
        InterpolatorType * GetMovingImageInterpolator(void)
        { return m_MovingImageInterpolator; }

        /** This class uses a constant timestep of 1. */
        virtual TimeStepType ComputeGlobalTimeStep( void *itkNotUsed(GlobalData) )
        const
        { return m_TimeStep; }

        /** Return a pointer to a global data structure that is passed to
         * this object from the solver at each calculation.  */
        virtual void * GetGlobalDataPointer() const
        {
            GlobalDataStruct *global = new GlobalDataStruct();

            global->m_SumOfSquaredDifference  = 0.0;
            global->m_NumberOfPixelsProcessed = 0L;
            global->m_SumOfSquaredChange      = 0;
            return global;
        }

        /** Release memory for global data structure. */
        virtual void ReleaseGlobalDataPointer(void *GlobalData) const;

        /** Set the object's state before each iteration. */
        virtual void InitializeIteration();

        /** This method is called by a finite difference solver image filter at
         * each pixel that does not lie on a data set boundary */
        virtual PixelType  ComputeUpdate( const NeighborhoodType & neighborhood,
                                         void *globalData,
                                         const FloatOffsetType & offset =
                                         FloatOffsetType(0.0) );

        /** Get the metric value. The metric value is the mean square difference
         * in intensity between the fixed image and transforming moving image
         * computed over the the overlapping region between the two images. */
        virtual double GetMetric() const
        { return m_Metric; }

        /** Get the rms change in displacement field. */
        virtual double GetRMSChange() const
        { return m_RMSChange; }

        /** Select if the fixed image or moving image gradient is used for
         * the computating the demon forces. The fixed image gradient is used
         * by default. */
        virtual void SetUseMovingImageGradient(bool flag)
        { m_UseMovingImageGradient = flag; }
        virtual bool GetUseMovingImageGradient() const
        { return m_UseMovingImageGradient; }

        /** Set/Get the threshold below which the absolute difference of
         * intensity yields a match. When the intensities match between a
         * moving and fixed image pixel, the update vector (for that
         * iteration) will be the zero vector. Default is 0.001. */
        virtual void SetIntensityDifferenceThreshold(double);

        virtual double GetIntensityDifferenceThreshold() const;

    protected:
        DemonsRegistrationFunction();
        ~DemonsRegistrationFunction() {}
        void PrintSelf(std::ostream & os, itk::Indent indent) const;

        /** FixedImage image neighborhood iterator type. */
        typedef itk::ConstNeighborhoodIterator< FixedImageType >
        FixedImageNeighborhoodIteratorType;

        /** A global data type for this class of equation. Used to store
         * information for computing the metric. */
        struct GlobalDataStruct {
            double m_SumOfSquaredDifference;
            itk::SizeValueType m_NumberOfPixelsProcessed;
            double m_SumOfSquaredChange;
        };

    private:
        DemonsRegistrationFunction(const Self &); //purposely not implemented
        void operator=(const Self &);             //purposely not implemented

        /** Cache fixed image information. */
        //SpacingType                  m_FixedImageSpacing;
        //PointType                    m_FixedImageOrigin;
        PixelType m_ZeroUpdateReturn;
        double    m_Normalizer;

        /** Function to compute derivatives of the fixed image. */
        GradientCalculatorPointer m_FixedImageGradientCalculator;

        /** Function to compute derivatives of the moving image. */
        MovingImageGradientCalculatorPointer m_MovingImageGradientCalculator;
        bool                                 m_UseMovingImageGradient;
        
        /** Function to interpolate the moving image. */
        InterpolatorPointer m_MovingImageInterpolator;
        
        /** The global timestep. */
        TimeStepType m_TimeStep;
        
        /** Threshold below which the denominator term is considered zero. */
        double m_DenominatorThreshold;
        
        /** Threshold below which two intensity value are assumed to match. */
        double m_IntensityDifferenceThreshold;
        
        /** The metric value is the mean square difference in intensity between
         * the fixed image and transforming moving image computed over the
         * the overlapping region between the two images. */
        mutable double        m_Metric;
        mutable double        m_SumOfSquaredDifference;
        mutable itk::SizeValueType m_NumberOfPixelsProcessed;
        mutable double        m_RMSChange;
        mutable double        m_SumOfSquaredChange;
        
        /** Mutex lock to protect modification to metric. */
        mutable itk::SimpleFastMutexLock m_MetricCalculationLock;
    };

#pragma mark -

    /** \class DemonsRegistrationFilter
     * \brief Deformably register two images using the demons algorithm.
     *
     * DemonsRegistrationFilter implements the demons deformable algorithm that
     * register two images by computing the displacement field which will map a
     * moving image onto a fixed image.
     *
     * A displacement field is represented as a image whose pixel type is some
     * vector type with at least N elements, where N is the dimension of
     * the fixed image. The vector type must support element access via operator
     * []. It is assumed that the vector elements behave like floating point
     * scalars.
     *
     * This class is templated over the fixed image type, moving image type
     * and the displacement field type.
     *
     * The input fixed and moving images are set via methods SetFixedImage
     * and SetMovingImage respectively. An initial displacement field maybe set via
     * SetInitialDisplacementField or SetInput. If no initial field is set,
     * a zero field is used as the initial condition.
     *
     * The algorithm has one parameters: the number of iteration to be performed.
     *
     * The output displacement field can be obtained via methods GetOutput
     * or GetDisplacementField.
     *
     * This class make use of the finite difference solver hierarchy. Update
     * for each iteration is computed in DemonsRegistrationFunction.
     *
     * \warning This filter assumes that the fixed image type, moving image type
     * and displacement field type all have the same number of dimensions.
     *
     * \sa DemonsRegistrationFunction
     * \ingroup DeformableImageRegistration MultiThreaded
     * \ingroup ITKPDEDeformableRegistration
     */
    template<class TFixedImage, class TMovingImage, class TDisplacementField>
    class DemonsRegistrationFilter:
    public itk::PDEDeformableRegistrationFilter< TFixedImage, TMovingImage,
    TDisplacementField >
    {
    public:
        /** Standard class typedefs. */
        typedef DemonsRegistrationFilter Self;
        typedef itk::PDEDeformableRegistrationFilter<TFixedImage, TMovingImage, TDisplacementField> Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(DemonsRegistrationFilter, itk::PDEDeformableRegistrationFilter);

        /** Inherit types from superclass. */
        typedef typename Superclass::TimeStepType TimeStepType;

        /** FixedImage image type. */
        typedef typename Superclass::FixedImageType    FixedImageType;
        typedef typename Superclass::FixedImagePointer FixedImagePointer;

        /** MovingImage image type. */
        typedef typename Superclass::MovingImageType    MovingImageType;
        typedef typename Superclass::MovingImagePointer MovingImagePointer;

        /** displacement field type. */
        typedef typename Superclass::DisplacementFieldType    DisplacementFieldType;
        typedef typename Superclass::DisplacementFieldPointer DisplacementFieldPointer;

#ifdef ITKV3_COMPATIBILITY
        typedef typename Superclass::DeformationFieldType    DeformationFieldType;
        typedef typename Superclass::DeformationFieldPointer DeformationFieldPointer;
#endif

        /** FiniteDifferenceFunction type. */
        typedef typename Superclass::FiniteDifferenceFunctionType
        FiniteDifferenceFunctionType;

        /** DemonsRegistrationFilterFunction type. */
        typedef DemonsRegistrationFunction<FixedImageType, MovingImageType, DisplacementFieldType> DemonsRegistrationFunctionType;

        /** Get the metric value. The metric value is the mean square difference
         * in intensity between the fixed image and transforming moving image
         * computed over the the overlapping region between the two images.
         * This is value is only available for the previous iteration and
         * NOT the current iteration. */
        virtual double GetMetric() const;

        /** Switch between using the fixed image and moving image gradient
         * for computing the displacement field updates. */
        itkSetMacro(UseMovingImageGradient, bool);
        itkGetConstMacro(UseMovingImageGradient, bool);
        itkBooleanMacro(UseMovingImageGradient);

        /** Set/Get the threshold below which the absolute difference of
         * intensity yields a match. When the intensities match between a
         * moving and fixed image pixel, the update vector (for that
         * iteration) will be the zero vector. Default is 0.001. */
        virtual void SetIntensityDifferenceThreshold(double);

        virtual double GetIntensityDifferenceThreshold() const;
        
    protected:
        DemonsRegistrationFilter();
        // ~DemonsRegistrationFilter() {} default implementation ok
        void PrintSelf(std::ostream & os, itk::Indent indent) const;
        
        /** Initialize the state of filter and equation before each iteration. */
        virtual void InitializeIteration();
        
        /** Apply update. */
        virtual void ApplyUpdate(const TimeStepType& dt);
        
        /** Override VeriyInputInformation() since this filter's inputs do
         * not need to occoupy the same physical space.
         *
         * \sa ProcessObject::VerifyInputInformation
         */
        virtual void VerifyInputInformation() {}
        
    private:
        DemonsRegistrationFilter(const Self &); //purposely not implemented
        void operator=(const Self &);           //purposely not implemented
        
        bool m_UseMovingImageGradient;
    };

}

#endif /* defined(__ParticleGuidedRegistration__piPatchCompare__) */