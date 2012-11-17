//
//  ImageViewManager.cpp
//  laplacePDE
//
//  Created by Joohwi Lee on 11/13/12.
//
//

#include "myImageContainer.h"
#include "itkExtractImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageIO.h"
#include "itkARGBColorFunction.h"
#include "itkScalarToARGBColormapImageFilter.h"

void ImageContainer::SetAlpha(int alpha) {
    if (alpha == m_Alpha) {
        return;
    }
    m_Alpha = alpha;
    for (int j = 0; j < VDim; j++) {
        UpdateBitmap(j);
    }
}

void ImageContainer::SetLabelAlpha(int alpha) {
    if (alpha == m_LabelAlpha) {
        return;
    }
    m_LabelAlpha = alpha;
    for (int j = 0; j < VDim; j++) {
        UpdateLabelBitmap(j);
    }
}

void ImageContainer::SetSliceIndex(int dim, int idx) {
    if (idx < 0) {
        idx = 0;
    } else if (idx >= m_MaxSliceIndexes[dim]) {
        idx = m_MaxSliceIndexes[dim] - 1;
    }
    if (m_SliceIndexes[dim] == idx) {
        return;
    }
    m_SliceIndexes[dim] = idx;
    if (HasImage()) {
        UpdateSlice(dim);
        UpdateBitmap(dim);
    }
    if (HasLabel()) {
        UpdateLabelSlice(dim);
        UpdateLabelBitmap(dim);
    }
}

void ImageContainer::LoadImage(const char* filename) {
    itkcmds::itkImageIO<ImageType> io;
    SetImage(io.ReadImageT(filename));
}

void ImageContainer::LoadLabel(const char* filename) {
    itkcmds::itkImageIO<LabelType> io;
    SetLabel(io.ReadImageT(filename));
}

void ImageContainer::SetImage(ImageType::Pointer image) {
    if (image.IsNull()) {
        return;
    }
    
    m_Image = image;
    ImageType::SizeType sz = m_Image->GetBufferedRegion().GetSize();
    for (int i = 0; i < VDim; i++) {
        m_MaxSliceIndexes[i] = sz[i];
        m_SliceIndexes[i] = m_MaxSliceIndexes[i] / 2;
    }

    typedef itk::StatisticsImageFilter<ImageType> StatisticsImageFilterType;
    StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New ();
    statisticsImageFilter->SetInput(m_Image);
    statisticsImageFilter->Update();

    m_IntensityStats[0] = statisticsImageFilter->GetMinimum();
    m_IntensityStats[1] = statisticsImageFilter->GetMaximum();
    m_IntensityStats[2] = statisticsImageFilter->GetMean();
    m_IntensityStats[3] = statisticsImageFilter->GetSigma();
}

void ImageContainer::UpdateSlice(int dim) {
    if (m_Image.IsNull()) {
        return;
    }

    typedef itk::ExtractImageFilter<ImageType, SliceType> SliceExtractor;
    SliceExtractor::Pointer slicer = SliceExtractor::New();
    ImageType::RegionType sliceRegion;
    ImageType::RegionType::SizeType sliceSize;
    ImageType::RegionType::IndexType sliceIndex;

    for (int i = 0; i < ImageType::ImageDimension; i++) {
        if (i == dim) {
            sliceSize[i] = 0;
            sliceIndex[i] = m_SliceIndexes[i];
        } else {
            sliceSize[i] = m_MaxSliceIndexes[i];
            sliceIndex[i] = 0;
        }
    }

    sliceRegion.SetSize(sliceSize);
    sliceRegion.SetIndex(sliceIndex);

    slicer->SetInput(m_Image);
    slicer->SetExtractionRegion(sliceRegion);
    slicer->SetDirectionCollapseToGuess();
    slicer->Update();
    m_Slices[dim] = slicer->GetOutput();
}

void ImageContainer::UpdateBitmap(int dim) {
    if (m_Slices[dim].IsNull()) {
        UpdateSlice(dim);
    }
    typedef itk::ScalarToARGBColormapImageFilter<SliceType, RGBAImageType> ScalarToRGBFilter;
    ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
    rgbFilter->SetInput(m_Slices[dim]);
    rgbFilter->UseManualScalingOn();
    rgbFilter->SetAlphaValue(m_Alpha);
    rgbFilter->SetMinimumValue(m_IntensityStats[0]);
    rgbFilter->SetMaximumValue(m_IntensityStats[1]);
    rgbFilter->Update();
    m_Bitmaps[dim] = rgbFilter->GetOutput();
}


QPixmap ImageContainer::GetPixmap(int dim) {
    RGBAImageType::Pointer bitmap = m_Bitmaps[dim];
    if (bitmap.IsNull()) {
        UpdateBitmap(dim);
        bitmap = m_Bitmaps[dim];
    }
    RGBAImageType::SizeType bitmapSz = bitmap->GetBufferedRegion().GetSize();
    QImage qImg = QImage((unsigned char*) bitmap->GetBufferPointer(),
                         bitmapSz[0], bitmapSz[1], QImage::Format_ARGB32);
    return QPixmap::fromImage(qImg);
}

void ImageContainer::SetLabel(LabelType::Pointer label) {
    if (label.IsNull()) {
        return;
    }

    m_Label = label;
    LabelType::SizeType sz = m_Label->GetBufferedRegion().GetSize();
    for (int i = 0; i < VDim; i++) {
        if (m_MaxSliceIndexes[i] == 0) {
            m_MaxSliceIndexes[i] = sz[i];
            m_SliceIndexes[i] = m_MaxSliceIndexes[i] / 2;
        }
    }

    typedef itk::StatisticsImageFilter<LabelType> StatisticsImageFilterType;
    StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New ();
    statisticsImageFilter->SetInput(m_Label);
    statisticsImageFilter->Update();

    m_LabelStats[0] = statisticsImageFilter->GetMinimum();
    m_LabelStats[1] = statisticsImageFilter->GetMaximum();
    m_LabelStats[2] = statisticsImageFilter->GetMean();
    m_LabelStats[3] = statisticsImageFilter->GetSigma();
}

void ImageContainer::UpdateLabelSlice(int dim) {
    if (m_Label.IsNull()) {
        return;
    }

    typedef itk::ExtractImageFilter<LabelType, LabelSliceType> SliceExtractor;
    SliceExtractor::Pointer slicer = SliceExtractor::New();
    LabelType::RegionType sliceRegion;
    LabelType::RegionType::SizeType sliceSize;
    LabelType::RegionType::IndexType sliceIndex;

    for (int i = 0; i < LabelType::ImageDimension; i++) {
        if (i == dim) {
            sliceSize[i] = 0;
            sliceIndex[i] = m_SliceIndexes[i];
        } else {
            sliceSize[i] = m_MaxSliceIndexes[i];
            sliceIndex[i] = 0;
        }
    }

    sliceRegion.SetSize(sliceSize);
    sliceRegion.SetIndex(sliceIndex);

    slicer->SetInput(m_Label);
    slicer->SetExtractionRegion(sliceRegion);
    slicer->SetDirectionCollapseToGuess();
    slicer->Update();
    m_LabelSlices[dim] = slicer->GetOutput();
}

void ImageContainer::UpdateLabelBitmap(int dim) {
    if (m_LabelSlices[dim].IsNull()) {
        UpdateLabelSlice(dim);
    }
    typedef itk::ScalarToARGBColormapImageFilter<LabelSliceType, RGBAImageType> ScalarToRGBFilter;
    ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
    rgbFilter->SetInput(m_LabelSlices[dim]);
    rgbFilter->UseManualScalingOn();
    rgbFilter->SetAlphaValue(m_LabelAlpha);
    rgbFilter->SetMinimumValue(m_LabelStats[0]);
    rgbFilter->SetMaximumValue(m_LabelStats[1]);
    rgbFilter->Update();
    m_LabelBitmaps[dim] = rgbFilter->GetOutput();
}


QPixmap ImageContainer::GetLabelPixmap(int dim) {
    RGBAImageType::Pointer bitmap = m_LabelBitmaps[dim];
    if (bitmap.IsNull()) {
        UpdateLabelBitmap(dim);
    }
    bitmap = m_LabelBitmaps[dim];
    RGBAImageType::SizeType bitmapSz = bitmap->GetBufferedRegion().GetSize();
    QImage qImg = QImage((unsigned char*) bitmap->GetBufferPointer(),
                         bitmapSz[0], bitmapSz[1], QImage::Format_ARGB32);
    return QPixmap::fromImage(qImg);
}
