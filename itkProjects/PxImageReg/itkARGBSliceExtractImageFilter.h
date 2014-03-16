
//
//  itkARGBSliceExtractImageFilter.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/2/13.
//
//

#ifndef __ParticleGuidedRegistration__itkARGBSliceExtractImageFilter__
#define __ParticleGuidedRegistration__itkARGBSliceExtractImageFilter__

#include <iostream>
#include <itkObject.h>
#include <itkImage.h>
#include <itkRGBAPixel.h>
#include <itkImageSliceIteratorWithIndex.h>
#include "itkARGBColorFunction.h"
#include "piImageIO.h"

namespace itk {

    typedef itk::RGBAPixel<unsigned char> ARGBPixel;
    typedef itk::Image<ARGBPixel,2> ARGBSlice;
    typedef itk::ARGBSlice::Pointer ARGBSlicePointer;

    enum AxisFlip { NoFlip, MajorFlip, MinorFlip, MajorMinorFlip };

    template< typename TImage >
    class ITK_EXPORT ImageSliceLineIteratorWithIndex: public ImageConstIteratorWithIndex< TImage >
    {
    public:
        /** Standard class typedefs. */
        typedef ImageSliceLineIteratorWithIndex      Self;
        typedef ImageConstIteratorWithIndex< TImage > Superclass;

        /** Inherit types from the superclass */
        typedef typename Superclass::IndexType             IndexType;
        typedef typename Superclass::SizeType              SizeType;
        typedef typename Superclass::OffsetType            OffsetType;
        typedef typename Superclass::RegionType            RegionType;
        typedef typename Superclass::ImageType             ImageType;
        typedef typename Superclass::PixelContainer        PixelContainer;
        typedef typename Superclass::PixelContainerPointer PixelContainerPointer;
        typedef typename Superclass::InternalPixelType     InternalPixelType;
        typedef typename Superclass::PixelType             PixelType;
        typedef typename Superclass::AccessorType          AccessorType;

        /** Default constructor. Needed since we provide a cast constructor. */
        ImageSliceLineIteratorWithIndex():ImageConstIteratorWithIndex< TImage >() {}

        /** Constructor establishes an iterator to walk a particular image and a
         * particular region of that image. */
        ImageSliceLineIteratorWithIndex(const ImageType *ptr,
                                         const RegionType & region):
        ImageConstIteratorWithIndex< TImage >(ptr, region)
        {
            m_Direction_A = 0;
            m_Direction_B = 1;
        }

        /** Constructor that can be used to cast from an ImageIterator to an
         * ImageSliceConstIteratorWithIndex. Many routines return an ImageIterator, but for a
         * particular task, you may want an ImageSliceConstIteratorWithIndex.  Rather than
         * provide overloaded APIs that return different types of Iterators, itk
         * returns ImageIterators and uses constructors to cast from an
         * ImageIterator to a ImageSliceConstIteratorWithIndex. */
        ImageSliceLineIteratorWithIndex(const ImageConstIteratorWithIndex< TImage > & it)
        { this->ImageConstIteratorWithIndex< TImage >::operator=(it); }

        /** Go to the first pixel of the current slice */
        void GoToBeginOfSlice(void);

        /** Go to the next slice
         * \sa operator++ \sa EndOfLine \sa End */
        void NextSlice(void);

        void GoToLineBegin(void);
        void GoToLineEnd(void);
        void PreviousLine(void);
        void NextLine(void);


        /** Go to the previous slice
         */
        void PreviousSlice(void);

        /** Test if the index is at the end of line */
        bool IsAtEndOfLine(void);

        /** Test if the index is at the end of the slice */
        bool IsAtEndOfSlice(void);

        /** Test if the index is at the begin of line */
        bool IsAtReverseEndOfLine(void);

        /** Test if the index is at the begin of the slice */
        bool IsAtReverseEndOfSlice(void);

        /** Test if the index of second direction is at the first index of the second direction */
        bool IsAtFirstLine(void);

        /** Test if the index of second direction is at the last index of the second direction */
        bool IsAtLastLine(void);

        /** Set the fastest direction of movement */
        void SetFirstDirection(unsigned int direction);

        /** Set the second fastest direction of movement */
        void SetSecondDirection(unsigned int direction);
        
        /** Increment (prefix) the selected dimension.
         * No bounds checking is performed.
         * \sa operator-- \sa GetIndex */
        inline Self & operator++();
        
        /** Decrement (prefix) the selected dimension.
         * No bounds checking is performed.
         * \sa operator++ \sa GetIndex */
        inline Self & operator--();
        
    private:
        SizeValueType m_PixelJump;
        SizeValueType m_LineJump;
        unsigned int  m_Direction_A;
        unsigned int  m_Direction_B;
    };


    template <class TImage>
    bool
    ImageSliceLineIteratorWithIndex<TImage>
    ::IsAtFirstLine() {
        return this->m_PositionIndex[m_Direction_B] <= this->m_BeginIndex[m_Direction_B];
    }

    template <class TImage>
    bool
    ImageSliceLineIteratorWithIndex<TImage>
    ::IsAtLastLine() {
        return this->m_PositionIndex[m_Direction_B] >= this->m_EndIndex[m_Direction_B];
    }

    //----------------------------------------------------------------------
    //  Advance to the begin of line
    //----------------------------------------------------------------------
    template< class TImage >
    void
    ImageSliceLineIteratorWithIndex< TImage >
    ::GoToLineBegin(void)
    {
        // Move to beginning of line
        this->m_Position -= m_PixelJump
        * ( this->m_PositionIndex[m_Direction_A] - this->m_BeginIndex[m_Direction_A] );
        this->m_PositionIndex[m_Direction_A] = this->m_BeginIndex[m_Direction_A];
    }

    //----------------------------------------------------------------------
    //  Advance to Next Line
    //----------------------------------------------------------------------
    template< class TImage >
    void
    ImageSliceLineIteratorWithIndex< TImage >
    ::NextLine(void)
    {
        // Move to next line
        this->m_PositionIndex[m_Direction_B]++;
        this->m_Position += m_LineJump;
    }

    //----------------------------------------------------------------------
    //  Advance to Previous Line
    //----------------------------------------------------------------------
    template< class TImage >
    void
    ImageSliceLineIteratorWithIndex< TImage >
    ::PreviousLine(void)
    {
        // Move to previous line
        this->m_PositionIndex[m_Direction_B]--;
        this->m_Position -= m_LineJump;
    }

    
    //----------------------------------------------------------------------
    //  Advance to the end of line
    //----------------------------------------------------------------------
    template< class TImage >
    void
    ImageSliceLineIteratorWithIndex< TImage >
    ::GoToLineEnd(void)
    {
        // Move to end of line
        this->m_Position += m_PixelJump
        * ( this->m_EndIndex[m_Direction_A] - this->m_PositionIndex[m_Direction_A] );
        this->m_PositionIndex[m_Direction_A] = this->m_EndIndex[m_Direction_A] - 1;
    }

    //----------------------------------------------------------------------
    //  Go to the first pixel of the current slice
    //----------------------------------------------------------------------
    template< class TImage >
    void
    ImageSliceLineIteratorWithIndex< TImage >
    ::GoToBeginOfSlice(void)
    {
        // Move to beginning of Slice
        this->m_PositionIndex[m_Direction_B] = this->m_BeginIndex[m_Direction_B];
        this->m_Position -= m_LineJump
        * ( this->m_EndIndex[m_Direction_B] - this->m_BeginIndex[m_Direction_B] );
    }

    //----------------------------------------------------------------------
    //  Advance to next slice
    //----------------------------------------------------------------------
    template< class TImage >
    void
    ImageSliceLineIteratorWithIndex< TImage >
    ::NextSlice(void)
    {
        // Move to beginning of Slice
        this->m_Position -= m_LineJump * (this->m_PositionIndex[m_Direction_B] - this->m_BeginIndex[m_Direction_B]);
        this->m_PositionIndex[m_Direction_B] = this->m_BeginIndex[m_Direction_B];

        for (unsigned int n = 0; n < TImage::ImageDimension; n++) {
            this->m_Remaining = false;

            if ( n == m_Direction_B || n == m_Direction_A ) {
                continue;
            }

            this->m_PositionIndex[n]++;
            if ( this->m_PositionIndex[n] < this->m_EndIndex[n] )
            {
                this->m_Position += this->m_OffsetTable[n];
                this->m_Remaining = true;
                break;
            }
            else
            {
                this->m_Position -= this->m_OffsetTable[n + 1] - this->m_OffsetTable[n];
                this->m_PositionIndex[n] = this->m_BeginIndex[n];
            }
        }
    }

    //----------------------------------------------------------------------
    //  Go Back to previous slice
    //----------------------------------------------------------------------
    template< class TImage >
    void
    ImageSliceLineIteratorWithIndex< TImage >
    ::PreviousSlice(void)
    {
        // Move to end of Slice
        this->m_PositionIndex[m_Direction_B] = this->m_EndIndex[m_Direction_B] - 1;
        this->m_Position += m_LineJump * ( this->m_EndIndex[m_Direction_B] - this->m_BeginIndex[m_Direction_B] );

        for ( unsigned int n = 0; n < TImage::ImageDimension; n++ )
        {
            this->m_Remaining = false;

            if ( n == m_Direction_B || n == m_Direction_A )
            {
                continue;
            }

            this->m_PositionIndex[n]--;
            if ( this->m_PositionIndex[n] >= this->m_BeginIndex[n] )
            {
                this->m_Position -= this->m_OffsetTable[n];
                this->m_Remaining = true;
                break;
            }
            else
            {
                this->m_Position += this->m_OffsetTable[n + 1] - this->m_OffsetTable[n];
                this->m_PositionIndex[n] = this->m_EndIndex[n] - 1;
            }
        }
    }

    //----------------------------------------------------------------------
    //  Test for end of line
    //----------------------------------------------------------------------
    template< class TImage >
    bool
    ImageSliceLineIteratorWithIndex< TImage >
    ::IsAtEndOfLine(void)
    {
        return this->m_PositionIndex[m_Direction_A] >= this->m_EndIndex[m_Direction_A];
    }

    //----------------------------------------------------------------------
    //  Test for end of slice
    //----------------------------------------------------------------------
    template< class TImage >
    bool
    ImageSliceLineIteratorWithIndex< TImage >
    ::IsAtEndOfSlice(void)
    {
        return this->m_PositionIndex[m_Direction_B] >= this->m_EndIndex[m_Direction_B];
    }

    //----------------------------------------------------------------------
    //  Test for begin of line
    //----------------------------------------------------------------------
    template< class TImage >
    bool
    ImageSliceLineIteratorWithIndex< TImage >
    ::IsAtReverseEndOfLine(void)
    {
        return this->m_PositionIndex[m_Direction_A] < this->m_BeginIndex[m_Direction_A];
    }

    //----------------------------------------------------------------------
    //  Test for begin of slice
    //----------------------------------------------------------------------
    template< class TImage >
    bool
    ImageSliceLineIteratorWithIndex< TImage >
    ::IsAtReverseEndOfSlice(void)
    {
        return this->m_PositionIndex[m_Direction_B] < this->m_BeginIndex[m_Direction_B];
    }

    //----------------------------------------------------------------------
    //  Select the fastest changing direction
    //----------------------------------------------------------------------
    template< class TImage >
    void
    ImageSliceLineIteratorWithIndex< TImage >
    ::SetFirstDirection(unsigned int direction)
    {
        if ( direction >= TImage::ImageDimension )
        {
            itkGenericExceptionMacro(
                                     << "In image of dimension " << TImage::ImageDimension << " Direction " << direction << " sas selected");
        }
        m_Direction_A = direction;
        m_PixelJump = this->m_OffsetTable[m_Direction_A];
    }

    //----------------------------------------------------------------------
    //  Select the second fastest changing direction
    //----------------------------------------------------------------------
    template< class TImage >
    void
    ImageSliceLineIteratorWithIndex< TImage >
    ::SetSecondDirection(unsigned int direction)
    {
        if ( direction >= TImage::ImageDimension )
        {
            itkGenericExceptionMacro(
                                     << "In image of dimension " << TImage::ImageDimension << " Direction " << direction << " sas selected");
        }
        m_Direction_B = direction;
        m_LineJump = this->m_OffsetTable[m_Direction_B];
    }
    
    //----------------------------------------------------------------------
    //  Advance along a line
    //----------------------------------------------------------------------
    template< class TImage >
    ImageSliceLineIteratorWithIndex< TImage > &
    ImageSliceLineIteratorWithIndex< TImage >
    ::operator++()
    {
        this->m_PositionIndex[m_Direction_A]++;
        this->m_Position += m_PixelJump;
        return *this;
    }
    
    //----------------------------------------------------------------------
    //  Go back along a line
    //----------------------------------------------------------------------
    template< class TImage >
    ImageSliceLineIteratorWithIndex< TImage > &
    ImageSliceLineIteratorWithIndex< TImage >
    ::operator--()
    {
        this->m_PositionIndex[m_Direction_A]--;
        this->m_Position -= m_PixelJump;
        return *this;
    }
    
    template <class A>
    class ITK_EXPORT ARGBSliceExtractImageFilter: public Object {
    public:
        typedef ARGBSliceExtractImageFilter Self;
        typedef Object Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);
        itkTypeMacro(ARGBSliceExtractImageFilter,Object);

        typedef typename A::Pointer APointer;
        typedef typename A::PixelType APixel;
        typedef typename A::RegionType ARegionType;

        typedef Function::ARGBColormapFunction<APixel, ARGBPixel> ColormapType;

        void SetColormap(ColormapEnumType);
        itkSetMacro(MinimumValue, APixel);
        itkGetConstMacro(MinimumValue, APixel);
        itkSetMacro(MaximumValue, APixel);
        itkGetConstMacro(MaximumValue, APixel);

    public:
        void SetSliceAxis(int major, int minor);
        void SetSlice(int sliceIdx);
        void SetInput(APointer inputImage);
        void SetOutputImage(ARGBSlicePointer outputImage);
        void SetOutputFlip(AxisFlip);
        void SetIntensityRange(APixel min, APixel max);
        void Update();
        double Rescale(APixel val);
        ARGBSlicePointer GetOutput();

    protected:
        ARGBSliceExtractImageFilter();
        virtual ~ARGBSliceExtractImageFilter() {}
        void AllocateOutput();

    private:
        typedef itk::ImageSliceLineIteratorWithIndex<A> InputIterator;
        typedef itk::ImageRegionIterator<ARGBSlice> OutputIterator;

        void UpdateMajorFlip(InputIterator& iter, OutputIterator& outputIter);
        void UpdateMinorFlip(InputIterator& iter, OutputIterator& outputIter);
        void UpdateMajorMinorFlip(InputIterator& iter, OutputIterator& outputIter);
        void UpdateNoFlip(InputIterator& iter, OutputIterator& outputIter);

        ARGBSliceExtractImageFilter(const Self &); //purposely not implemented
        void operator=(const Self &);                 //purposely not implemented
        typename ColormapType::Pointer m_Colormap;

        int m_majorAxis;
        int m_minorAxis;
        int m_sliceAxis;
        int m_sliceIndex;

        APointer m_inputImage;
        ARegionType m_extractRegion;
        ARGBSlicePointer m_outputBuffer;

        APixel m_minimumValue;
        APixel m_maximumValue;
        AxisFlip m_outputFlip;

        void PrintError(const char* msg) {
            std::cout << msg << std::endl;
        }
    };

#pragma mark - 
#pragma mark Public Template Class Implementation
    template<class A> ARGBSliceExtractImageFilter<A>::ARGBSliceExtractImageFilter() {
        m_minorAxis = 0;
        m_majorAxis = 1;
        m_sliceAxis = 2;
        m_sliceIndex = 0;
        m_outputFlip = NoFlip;
        m_minimumValue = 0;
        m_maximumValue = 1;
    }

    template<class A> void ARGBSliceExtractImageFilter<A>::SetSliceAxis(int minor, int major) {
        m_minorAxis = minor;
        m_majorAxis = major;

        // identify slice normal axis
        int axes[3] = { 0, 1, 2 };
        axes[m_majorAxis] = 3;
        axes[m_minorAxis] = 3;
        for (int i = 0; i < 3; i++) {
            if (axes[i] < 3) {
                m_sliceAxis = axes[i];
                break;
            }
        }
    }

    template<class A> void ARGBSliceExtractImageFilter<A>::SetSlice(int sliceIdx) {
        m_sliceIndex = sliceIdx;
    }

    template<class A> void ARGBSliceExtractImageFilter<A>::SetInput(APointer inputImage) {
        m_inputImage = inputImage;
    }


    template<class A> double ARGBSliceExtractImageFilter<A>::Rescale(APixel val) {
     
    }


    template<class A> void ARGBSliceExtractImageFilter<A>::SetOutputImage(ARGBSlicePointer outputImage) {
        m_outputBuffer = outputImage;
    }

    template<class A> void ARGBSliceExtractImageFilter<A>::SetOutputFlip(itk::AxisFlip flip) {
        m_outputFlip = flip;
    }

    template<class A> ARGBSlicePointer ARGBSliceExtractImageFilter<A>::GetOutput() {
        return m_outputBuffer;
    }

    template<class A> void ARGBSliceExtractImageFilter<A>::AllocateOutput() {
        pi::ImageIO<ARGBSlice> io;
        ARGBSlice::SizeType outputSize;
        outputSize[0] = m_inputImage->GetBuffereRegion().GetSize(m_majorAxis);
        outputSize[1] = m_inputImage->GetBuffereRegion().GetSize(m_minorAxis);
        m_outputBuffer = io.NewImageT(outputSize);
    }

    template<class A> void ARGBSliceExtractImageFilter<A>::Update() {
        if (m_inputImage.IsNull()) {
            PrintError("Input is Null.");
            return;
        }
        if (m_outputBuffer.IsNull()) {
            AllocateOutput();
            if (m_outputBuffer.IsNull()) {
                PrintError("Output allocation failed");
                return;
            }
        }

        InputIterator iter(m_inputImage, m_inputImage->GetRequestedRegion());
        OutputIterator outputIter(m_outputBuffer, m_outputBuffer->GetBufferedRegion());

        iter.SetFirstDirection(m_majorAxis);
        iter.SetSecondDirection(m_minorAxis);
        for (int i = 0; i < m_sliceIndex - 1; i++) {
            iter.NextSlice();
        }
        switch (m_outputFlip) {
            case MajorFlip:
                UpdateMajorFlip(iter, outputIter);
                break;
            case MinorFlip:
                UpdateMinorFlip(iter, outputIter);
                break;
            case MajorMinorFlip:
                UpdateMajorMinorFlip(iter, outputIter);
                break;
            default:
                UpdateNoFlip(iter);
                break;
        }
    }

#pragma mark -
#pragma mark Private Template Class Implementation
    template<class A> void ARGBSliceExtractImageFilter<A>::UpdateMajorFlip(InputIterator &iter, OutputIterator &outputIter) {
        if (!iter.IsAtReverseEndOfSlice()) {
            iter.GoToBeginOfSlice();
        }
        while (!iter.IsAtEndOfSlice()) {
            iter.NextLine();
        }
        iter.PreviousLine();
        while (!iter.IsAtReverseEndOfSlice()) {
            iter.GoToLineBegin();
            while (!iter.IsAtEndOfLine()) {
                outputIter.Set(Rescale(iter.Get()));
                ++iter;
                ++outputIter;
            }
            iter.PreviousLine();
        }
    }

    template<class A> void ARGBSliceExtractImageFilter<A>::UpdateMinorFlip(InputIterator &iter, OutputIterator& outputIter) {
        if (!iter.IsAtReverseEndOfSlice()) {
            iter.GoToBeginOfSlice();
        }
        while (!iter.IsAtEndOfSlice()) {
            iter.GoToLineEnd();
            while (!iter.IsAtReverseEndOfLine()) {
                outputIter.Set(Rescale(iter.Get()));
                --iter;
                ++outputIter;
            }
            iter.NextLine();
        }
    }

    template<class A> void ARGBSliceExtractImageFilter<A>::UpdateMajorMinorFlip(InputIterator &iter, OutputIterator& outputIter) {
        ColormapType& colorFunc = *(m_Colormap.GetPointer());

        if (!iter.IsAtReverseEndOfSlice()) {
            iter.GoToBeginOfSlice();
        }
        while (!iter.IsAtEndOfSlice()) {
            iter.NextLine();
        }
        iter.PreviousLine();
        while (!iter.IsAtReverseEndOfSlice()) {
            iter.GoToLineEnd();
            while (!iter.IsAtReverseEndOfLine()) {
                outputIter.Set(colorFunc(Rescale(iter.Get())));
                --iter;
                ++outputIter;
            }
            iter.PreviousLine();
        }
    }

    template<class A> void ARGBSliceExtractImageFilter<A>::UpdateNoFlip(InputIterator &iter, OutputIterator& outputIter) {
        if (!iter.IsAtReverseEndOfSlice()) {
            iter.GoToBeginOfSlice();
        }
        while (!iter.IsAtEndOfSlice()) {
            while (!iter.IsAtEndOfLine()) {
                outputIter.Set(Rescale(iter.Get()));
                ++iter;
                ++outputIter;
            }
            iter.NextLine();
        }
    }
}
#endif /* defined(__ParticleGuidedRegistration__itkARGBSliceExtractImageFilter__) */