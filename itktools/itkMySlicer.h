//
//  itkMySlicer.h
//  itktools
//
//  Created by Joohwi Lee on 9/20/12.
//
//

#ifndef itktools_itkMySlicer_h
#define itktools_itkMySlicer_h

#include "itkObject.h"
#include "itkSmartPointer.h"
#include "itkExtractImageFilter.h"
#include "string"

/**
 * Extract a slice from a given view image and convert to ARGB format image
 *
 */
template <class TConverter>
class itkMySlicer : public itk::Object {
public:
    typedef itk::ExtractImageFilter<LabelType, SliceType> ExtractorType;

    typedef itkMySlicer Self;
    typedef itk::Object Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(itkMySlicer, itk::Object);
    itkGetMacro(SliceImage, BitmapType::Pointer);
    itkGetMacro(ViewImage, LabelType::Pointer);
    itkGetMacro(CurrentSliceIndex, int);
    itkGetMacro(Name, std::string);
    itkSetMacro(Name, std::string);
    itkGetMacro(Size, LabelType::SizeType);

    void SetInput(ImageType::Pointer originalImage) {
        m_SourceImage = originalImage;
        IntensityFilter::Pointer filter = IntensityFilter::New();
        filter->SetInput(m_SourceImage);
        filter->SetOutputMinimum(0);
        filter->SetOutputMaximum(255);
        filter->Update();
        SetLabel(filter->GetOutput());
    }

    void SetLabel(LabelType::Pointer labelImage) {
        m_ViewImage = labelImage;
        m_Size = m_ViewImage->GetBufferedRegion().GetSize();
    }


    bool UpdateSlice(int sliceIndex, int insideOpacity) {
        if (m_SliceImage.IsNotNull() && sliceIndex == m_CurrentSliceIndex && insideOpacity == m_LabelInsideOpacity) {
            return false;
        }
        m_CurrentSliceIndex = sliceIndex;
        m_LabelInsideOpacity = insideOpacity;
        GenerateSlice();
        return true;
    }

    int ComputeSliceAtCenter() {
        if (m_ViewImage.IsNull()) {
            return -1;
        }
        m_CurrentSliceIndex = m_ViewImage->GetBufferedRegion().GetSize()[2] / 2;
        return m_CurrentSliceIndex;
    }

    inline int* GetBitmapBuffer() {
        if (m_SliceImage.IsNull()) {
            return NULL;
        }
        return m_SliceImage->GetBufferPointer();
    }

protected:
    itkMySlicer() {
        m_CurrentSliceIndex = 0;
        m_LabelInsideOpacity = 100;
    }

    virtual ~itkMySlicer() {}

    void PrintSelf(std::ostream& os, itk::Indent indent) const {
        os << "Name: " << m_Name << std::endl;
        os << "Image: " << m_SourceImage << std::endl;
        os << "Slice: " << m_SliceImage << std::endl;
    }

private:
    itkMySlicer(const Self&);
    void operator=(const Self&);

    std::string m_Name;
    
    ImageType::Pointer m_SourceImage;
    LabelType::Pointer m_ViewImage;
    LabelType::SizeType m_Size;
    BitmapType::Pointer m_SliceImage;

    int m_CurrentSliceIndex;
    int m_LabelInsideOpacity;


    SliceType::Pointer ExtractSlice(LabelType::Pointer slicingImage, int slice = -1) {
        typedef itk::ExtractImageFilter<LabelType, SliceType> SliceExtractor;
        SliceExtractor::Pointer slicer = SliceExtractor::New();
        LabelType::RegionType sliceRegion;
        LabelType::SizeType sliceSize;
        LabelType::IndexType sliceIndex;
        sliceSize[0] = m_Size[0];
        sliceSize[1] = m_Size[1];
        sliceSize[2] = 0;
        sliceIndex[0] = sliceIndex[1] = 0;
        sliceIndex[2] = slice;
        sliceRegion.SetSize(sliceSize);
        sliceRegion.SetIndex(sliceIndex);
        slicer->SetInput(slicingImage);
        slicer->SetExtractionRegion(sliceRegion);
        slicer->SetDirectionCollapseToSubmatrix();
        slicer->Update();
        return slicer->GetOutput();
    }

    void GenerateSlice() {
        if (m_ViewImage.IsNotNull()) {
            m_CurrentSliceIndex = m_CurrentSliceIndex;
            SliceType::Pointer slice = ExtractSlice(m_ViewImage, m_CurrentSliceIndex);
            typename TConverter::Pointer filter = TConverter::New();
            filter->SetInput(slice);
            filter->GetFunctor().SetInsideAlpha(m_LabelInsideOpacity);
            filter->Update();
            m_SliceImage = filter->GetOutput();
        }
    }
};


typedef itkMySlicer<SlicerToBitmapFilter> GraySliceType;
typedef itkMySlicer<SlicerToLabelmapFilter> LabelSliceType;

#endif