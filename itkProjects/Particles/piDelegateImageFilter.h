//
//  piInplaceImageFilter.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/28/13.
//
//

#ifndef ParticleGuidedRegistration_piInplaceImageFilter_h
#define ParticleGuidedRegistration_piInplaceImageFilter_h

#include "piImageIO.h"
#include "itkInPlaceImageFilter.h"

namespace pi {
template <class TInputImage, class TDelegate, class TOutputImage = TInputImage>
class DelegateImageFilter: public itk::InPlaceImageFilter<TInputImage, TOutputImage> {
public:
    typedef DelegateImageFilter Self;
    typedef itk::InPlaceImageFilter<TInputImage> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    
    itkTypeMacro(DelegateImageFilter, itk::InPlaceImageFilter);
    itkNewMacro(Self);
    
    typedef typename Superclass::OutputImageType       OutputImageType;
    typedef typename Superclass::OutputImagePointer    OutputImagePointer;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename Superclass::OutputImagePixelType  OutputImagePixelType;
    
    typedef TInputImage                           InputImageType;
    typedef typename InputImageType::Pointer      InputImagePointer;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename InputImageType::RegionType   InputImageRegionType;
    typedef typename InputImageType::PixelType    InputImagePixelType;
    
    itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

    void SetDelegate(TDelegate* _arg) {
        this->m_Delegate = _arg;
    }
    
    itkGetConstMacro(Delegate, TDelegate*);
    
protected:
    DelegateImageFilter() {}
    virtual ~DelegateImageFilter() {};
    
    void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId) {
        m_Delegate->Process(this->GetOutput(), outputRegionForThread, threadId);
    }
    
private:
    DelegateImageFilter(const Self&);
    void operator=(const Self&);
    
    TDelegate* m_Delegate;
};
}

#endif
