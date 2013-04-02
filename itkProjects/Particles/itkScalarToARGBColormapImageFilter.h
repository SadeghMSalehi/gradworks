//
//  itkScalarToARGBColormapImageFilter.h
//  laplacePDE
//
//  Created by Joohwi Lee on 11/16/12.
//
//

#ifndef laplacePDE_itkScalarToARGBColormapImageFilter_h
#define laplacePDE_itkScalarToARGBColormapImageFilter_h

#include "itkARGBColorFunction.h"
#include "itkImageToImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"
#include "itkImageSliceConstIteratorWithIndex.h"

namespace itk {


    template< class TInputImage, class TOutputImage >
    class ITK_EXPORT ScalarToARGBColormapImageFilter:
    public ImageToImageFilter< TInputImage, TOutputImage >
    {
    public:
        /** Standard class typedefs. */
        typedef ScalarToARGBColormapImageFilter                  Self;
        typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
        typedef SmartPointer< Self >                            Pointer;
        typedef SmartPointer< const Self >                      ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(ScalarToARGBColormapImageFilter, ImageToImageFilter);

        /** Some typedefs. */
        typedef TInputImage                           InputImageType;
        typedef typename InputImageType::ConstPointer InputImagePointer;
        typedef typename InputImageType::RegionType   InputImageRegionType;
        typedef typename InputImageType::PixelType    InputImagePixelType;
        typedef TOutputImage                          OutputImageType;
        typedef typename OutputImageType::Pointer     OutputImagePointer;
        typedef typename OutputImageType::RegionType  OutputImageRegionType;
        typedef typename OutputImageType::PixelType   OutputImagePixelType;

        typedef Function::ARGBColormapFunction< InputImagePixelType, OutputImagePixelType> ColormapType;

        /**
         * Set/Get the colormap object.
         */
        itkSetObjectMacro(Colormap, ColormapType);
        itkGetObjectMacro(Colormap, ColormapType);

 

        void SetColormap(ColormapEnumType);

        void SetAlphaValue(int alpha) {
            m_Colormap->SetAlpha(alpha);
        }

        itkSetMacro(MinimumValue, InputImagePixelType);
        itkGetConstMacro(MinimumValue, InputImagePixelType);
        itkSetMacro(MaximumValue, InputImagePixelType);
        itkGetConstMacro(MaximumValue, InputImagePixelType);

        /**
         * Set/Get UseInputImageExtremaForScaling.  If 'true', the colormap uses the
         * min and max values from the image to scale appropriately.  Otherwise,
         * these values can be set in the colormap manually.
         */
        itkSetMacro(UseInputImageExtremaForScaling, bool);
        itkGetConstMacro(UseInputImageExtremaForScaling, bool);
        itkBooleanMacro(UseInputImageExtremaForScaling);

        itkSetMacro(UseManualScaling, bool);
        itkGetConstMacro(UseManualScaling, bool);
        itkBooleanMacro(UseManualScaling);

        itkSetMacro(UseIntensityWindow, bool);
        itkGetConstMacro(UseIntensityWindow, bool);
        itkBooleanMacro(UseIntensityWindow);
        
    protected:
        ScalarToARGBColormapImageFilter();
        virtual ~ScalarToARGBColormapImageFilter() {}

        void PrintSelf(std::ostream & os, Indent indent) const;

        virtual void GenerateOutputInformation()
        {
            // this methods is overloaded so that if the output image is a
            // VectorImage then the correct number of components are set.

            Superclass::GenerateOutputInformation();
            OutputImageType* output = this->GetOutput();

            if ( !output )
            {
                return;
            }
            if ( output->GetNumberOfComponentsPerPixel() != 4 )
            {
                output->SetNumberOfComponentsPerPixel( 4 );
            }
        }

        /** ScalarToARGBColormapImageFilter
         * can be implemented as a multithreaded filter.
         * Therefore, this implementation provides a ThreadedGenerateData() routine
         * which is called for each processing thread. The output image data is
         * allocated automatically by the superclass prior to calling
         * ThreadedGenerateData().  ThreadedGenerateData can only write to the
         * portion of the output image specified by the parameter
         * "outputRegionForThread"
         *
         * \sa ImageToImageFilter::ThreadedGenerateData(),
         *     ImageToImageFilter::GenerateData()  */
        void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                                  ThreadIdType threadId);

        /** Process to execute before entering the multithreaded section */
        void BeforeThreadedGenerateData();


    private:
        ScalarToARGBColormapImageFilter(const Self &); //purposely not implemented
        void operator=(const Self &);                 //purposely not implemented

        typename ColormapType::Pointer m_Colormap;

        bool m_UseInputImageExtremaForScaling;
        bool m_UseManualScaling;
        bool m_UseIntensityWindow;
        InputImagePixelType m_MinimumValue;
        InputImagePixelType m_MaximumValue;
    };

    /**
     * Constructor
     */
    template< class TInputImage, class TOutputImage >
    ScalarToARGBColormapImageFilter< TInputImage, TOutputImage >
    ::ScalarToARGBColormapImageFilter()
    {
        this->SetNumberOfRequiredInputs(1);
        this->m_UseInputImageExtremaForScaling = true;
        this->m_UseIntensityWindow = false;
        typedef Function::GreyColormapFunction<InputImagePixelType, OutputImagePixelType > DefaultColormapType;
        typename DefaultColormapType::Pointer greyColormap = DefaultColormapType::New();
        this->SetColormap(greyColormap);
    }

    /**
     * BeforeThreadedGenerateData
     */
    template< class TInputImage, class TOutputImage >
    void
    ScalarToARGBColormapImageFilter< TInputImage, TOutputImage >
    ::BeforeThreadedGenerateData()
    {
        if (this->m_UseManualScaling == true) {
            this->GetColormap()->SetMinimumInputValue(m_MinimumValue);
            this->GetColormap()->SetMaximumInputValue(m_MaximumValue);
        } else if ( this->m_UseInputImageExtremaForScaling == true) {
            ImageRegionConstIterator< InputImageType > It( this->GetInput(),
                                                          this->GetInput()->GetRequestedRegion() );

            InputImagePixelType minimumValue = NumericTraits< InputImagePixelType >::max();
            InputImagePixelType maximumValue = NumericTraits< InputImagePixelType >::min();

            for (It.GoToBegin(); !It.IsAtEnd(); ++It) {
                InputImagePixelType value = It.Get();
                if (value < minimumValue) {
                    minimumValue = value;
                } else if (value > maximumValue) {
                    maximumValue = value;
                }
            }

            this->GetColormap()->SetMinimumInputValue(minimumValue);
            this->GetColormap()->SetMaximumInputValue(maximumValue);
        }
    }

    /**
     * ThreadedGenerateData performs the pixel-wise mapping
     */
    template< class TInputImage, class TOutputImage >
    void
    ScalarToARGBColormapImageFilter< TInputImage, TOutputImage >
    ::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                           ThreadIdType threadId)
    {
        InputImagePointer  inputPtr = this->GetInput();
        OutputImagePointer outputPtr = this->GetOutput();

        // Define the portion of the input to walk for this thread, using
        // the CallCopyOutputRegionToInputRegion method allows for the input
        // and output images to be different dimensions
        InputImageRegionType inputRegionForThread;

        this->CallCopyOutputRegionToInputRegion(inputRegionForThread, outputRegionForThread);

        // Define the iterators
        ImageRegionConstIterator< TInputImage > inputIt(inputPtr, inputRegionForThread);
        ImageRegionIterator< TOutputImage >     outputIt(outputPtr, outputRegionForThread);

        ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels() );

        inputIt.GoToBegin();
        outputIt.GoToBegin();

        if (m_UseIntensityWindow) {
            while ( !inputIt.IsAtEnd() )
            {
                InputImagePixelType pixel = inputIt.Get();
                if (pixel < m_MinimumValue) {
                    pixel = m_MinimumValue;
                } else if (pixel > m_MaximumValue) {
                    pixel = m_MaximumValue;
                }
                outputIt.Set(this->m_Colormap->operator()(pixel));
                ++inputIt;
                ++outputIt;
                progress.CompletedPixel();  // potential exception thrown here
            }
        } else {
            while ( !inputIt.IsAtEnd() )
            {
                outputIt.Set( this->m_Colormap->operator()( inputIt.Get() ) );
                ++inputIt;
                ++outputIt;
                progress.CompletedPixel();  // potential exception thrown here
            }
        }
    }

    template< class TInputImage, class TOutputImage >
    void
    ScalarToARGBColormapImageFilter< TInputImage, TOutputImage >
    ::SetColormap(ColormapEnumType map)
    {
        switch ( map )
        {
            case Red:
            {
                typedef Function::RedColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Green:
            {
                typedef Function::GreenColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Blue:
            {
                typedef Function::BlueColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Grey:
            default:
            {
                typedef Function::GreyColormapFunction<InputImagePixelType, OutputImagePixelType> SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Hot:
            {
                typedef Function::HotColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Cool:
            {
                typedef Function::CoolColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Spring:
            {
                typedef Function::SpringColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Summer:
            {
                typedef Function::SummerColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Autumn:
            {
                typedef Function::AutumnColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Winter:
            {
                typedef Function::WinterColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Copper:
            {
                typedef Function::CopperColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case Jet:
            {
                typedef Function::JetColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case HSV:
            {
                typedef Function::HSVColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
            case OverUnder:
            {
                typedef Function::OverUnderColormapFunction<
                InputImagePixelType, OutputImagePixelType > SpecificColormapType;
                typename SpecificColormapType::Pointer colormap = SpecificColormapType::New();
                this->SetColormap(colormap);
                break;
            }
        }
    }

    template< class TInputImage, class TOutputImage >
    void
    ScalarToARGBColormapImageFilter< TInputImage, TOutputImage >
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        this->Superclass::PrintSelf(os, indent);
        os << indent << "Class Name: " << this->GetNameOfClass() << std::endl;
        if ( this->m_Colormap.IsNotNull() )
        {
            os << indent << "Colormap " << this->m_Colormap << std::endl;
        }
        else
        {
            os << indent << "Colormap is NULL " << std::endl;
        }
        os << indent << "Use Input Image Extrema for Scaling " << this->m_UseInputImageExtremaForScaling << std::endl;
    }
}; // end namespace itk


#endif
