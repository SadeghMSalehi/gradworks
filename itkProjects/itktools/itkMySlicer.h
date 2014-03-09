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
private:
    bool _modified;
    
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
    itkGetMacro(SourceImage, ImageType::Pointer);
    itkGetMacro(CurrentSliceIndex, int);
    itkGetMacro(CurrentSliceIndexX, int);
    itkGetMacro(CurrentSliceIndexY, int);
    itkGetMacro(CurrentSliceIndexZ, int);

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
        if (labelImage.IsNull()) {
            return;
        }
        m_ViewImage = labelImage;
        m_Size = m_ViewImage->GetBufferedRegion().GetSize();
    }


    bool UpdateSlice(int sliceIndex, int insideOpacity) {
        if (this->GetMTime() > m_ViewImage->GetMTime() && m_SliceImage.IsNotNull() && sliceIndex == m_CurrentSliceIndex && insideOpacity == m_LabelInsideOpacity) {
            return false;
        }
        m_CurrentSliceIndex = sliceIndex;
        m_LabelInsideOpacity = insideOpacity;
        GenerateSlice();
        return true;
    }

    bool needUpdate(int x, int y, int z, int insideOpacity) {
        return this->GetMTime() > m_ViewImage->GetMTime()
        && (m_SourceImage.IsNotNull() && (this->GetMTime() < m_SourceImage->GetMTime()))
        && (m_CurrentSliceIndexX != x || m_CurrentSliceIndexY != y || m_CurrentSliceIndexZ != z);
    }


    bool ZBoundCheck(int z) {
        return m_ViewImage.IsNotNull() && 0 <= z && z < m_ViewImage->GetBufferedRegion().GetSize()[2];
    }

    bool YBoundCheck(int y) {
        return m_ViewImage.IsNotNull() && 0 <= y && y < m_ViewImage->GetBufferedRegion().GetSize()[1];
    }

    bool XBoundCheck(int x) {
        return m_ViewImage.IsNotNull() && 0 <= x && x < m_ViewImage->GetBufferedRegion().GetSize()[0];
    }

    bool UpdateSlice(int x, int y, int z, int insideOpacity) {
        if (!needUpdate(x, y, z, insideOpacity)) {
            return false;
        }
        m_CurrentSliceIndexX = x;
        m_CurrentSliceIndexY = y;
        m_CurrentSliceIndexZ = z;
        m_LabelInsideOpacity = insideOpacity;
        GenerateSlice();
        return true;
    }
    
    int ComputeSliceAtCenter() {
        if (m_ViewImage.IsNull()) {
            return -1;
        }
        m_CurrentSliceIndex = m_ViewImage->GetBufferedRegion().GetSize()[2] / 2;

        m_CurrentSliceIndexX = m_ViewImage->GetBufferedRegion().GetSize()[0] / 2;
        m_CurrentSliceIndexY = m_ViewImage->GetBufferedRegion().GetSize()[1] / 2;
        m_CurrentSliceIndexZ = m_ViewImage->GetBufferedRegion().GetSize()[2] / 2;

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
        m_CurrentSliceIndex = m_CurrentSliceIndexX = m_CurrentSliceIndexY = m_CurrentSliceIndexZ = 0;
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
    BitmapType::Pointer m_SliceImage_X;
    BitmapType::Pointer m_SliceImage_Y;
    BitmapType::Pointer m_SliceImage_Z;

    int m_CurrentSliceIndex;
    int m_CurrentSliceIndexX, m_CurrentSliceIndexY, m_CurrentSliceIndexZ;

    int m_LabelInsideOpacity;


    SliceType::Pointer ExtractSliceZ(LabelType::Pointer slicingImage, int slice) {
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

    SliceType::Pointer ExtractSliceY(LabelType::Pointer slicingImage, int slice) {
        typedef itk::ExtractImageFilter<LabelType, SliceType> SliceExtractor;
        SliceExtractor::Pointer slicer = SliceExtractor::New();
        LabelType::RegionType sliceRegion;
        LabelType::SizeType sliceSize;
        LabelType::IndexType sliceIndex;
        sliceSize[0] = m_Size[0];
        sliceSize[1] = 0;
        sliceSize[2] = m_Size[2];
        sliceIndex[0] = sliceIndex[2] = 0;
        sliceIndex[1] = slice;
        sliceRegion.SetSize(sliceSize);
        sliceRegion.SetIndex(sliceIndex);
        slicer->SetInput(slicingImage);
        slicer->SetExtractionRegion(sliceRegion);
        slicer->SetDirectionCollapseToSubmatrix();
        slicer->Update();
        return slicer->GetOutput();
    }


    SliceType::Pointer ExtractSliceX(LabelType::Pointer slicingImage, int slice) {
        typedef itk::ExtractImageFilter<LabelType, SliceType> SliceExtractor;
        SliceExtractor::Pointer slicer = SliceExtractor::New();
        LabelType::RegionType sliceRegion;
        LabelType::SizeType sliceSize;
        LabelType::IndexType sliceIndex;
        sliceSize[0] = 0;
        sliceSize[1] = m_Size[1];
        sliceSize[2] = m_Size[2];
        sliceIndex[0] = slice;
        sliceIndex[1] = sliceIndex[2] = 0;
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
            SliceType::Pointer sliceZ = ExtractSliceZ(m_ViewImage, m_CurrentSliceIndex);
            typename TConverter::Pointer filterZ = TConverter::New();
            filterZ->SetInput(sliceZ);
            filterZ->GetFunctor().SetInsideAlpha(m_LabelInsideOpacity);
            filterZ->Update();
            m_SliceImage = m_SliceImage_Z = filterZ->GetOutput();

            SliceType::Pointer sliceX = ExtractSliceX(m_ViewImage, m_CurrentSliceIndexX);
            typename TConverter::Pointer filterX = TConverter::New();
            filterX->SetInput(sliceX);
            filterX->GetFunctor().SetInsideAlpha(m_LabelInsideOpacity);
            filterX->Update();
            m_SliceImage_X = filterX->GetOutput();

            SliceType::Pointer sliceY = ExtractSliceY(m_ViewImage, m_CurrentSliceIndexY);
            typename TConverter::Pointer filterY = TConverter::New();
            filterY->SetInput(sliceY);
            filterY->GetFunctor().SetInsideAlpha(m_LabelInsideOpacity);
            filterY->Update();
            m_SliceImage_Y = filterY->GetOutput();

            this->Modified();
        }
    }
};


typedef itkMySlicer<SlicerToBitmapFilter> GraySliceType;
typedef itkMySlicer<SlicerToLabelmapFilter> LabelSliceType;

#endif