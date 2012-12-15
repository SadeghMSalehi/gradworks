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
#include "itkCheckerBoardImageFilter.h"
#include "itkTransformToDisplacementFieldSource.h"
#include "itkWarpImageFilter.h"
#include "itkCastImageFilter.h"

ImageContainer::SliceDictionary g_DerivedSliceDictionary;

int ImageContainer::g_CurrentImage = 0;
int ImageContainer::g_CurrentLabel = 0;
int ImageContainer::g_CurrentView = 0;

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
    SetName(std::string(filename));
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

    // set spacing to 1, origin to 0, and direction to identity
    // to make its physical space and index space same
    double origin[3] = { 0., 0., 0. };
    double spacing[3] = { 1., 1., 1. };
    ImageType::DirectionType identityDir;
    identityDir.SetIdentity();
    image->SetSpacing(spacing);
    image->SetOrigin(origin);
    image->SetDirection(identityDir);

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

    ComputeTransformToIndexMatrix(m_WorldToImage);
    ComputeTransformToPhysicalPointMatrix(m_ImageToWorld);
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

SliceType::Pointer ImageContainer::GetLabelSliceAsSliceType(int dim) {
    typedef itk::CastImageFilter<LabelSliceType,SliceType> CastFilterType;
    CastFilterType::Pointer f = CastFilterType::New();
    if (dim == -1) {
        f->SetInput(GetLabelSlice());
    } else {
        f->SetInput(GetLabelSlice(dim));
    }
    f->Update();
    return f->GetOutput();
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

    // set spacing to 1, origin to 0, and direction to identity
    // to make its physical space and index space same
    double origin[3] = { 0., 0., 0. };
    double spacing[3] = { 1., 1., 1. };
    LabelType::DirectionType identityDir;
    identityDir.SetIdentity();
    label->SetSpacing(spacing);
    label->SetOrigin(origin);
    label->SetDirection(identityDir);

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

void ImageContainer::AddDerivedView(std::string name, RGBAImageType::Pointer slice) {
//    cout << "Name: " << name << endl << slice << endl;

    // debug: should always erase the previous key
    g_DerivedSliceDictionary.erase(name);
    g_DerivedSliceDictionary.insert(NameSlicePair(name, slice));

    if (m_EventCallback != NULL) {
        m_EventCallback->EventRaised(0xADDCE, 0);
    }
}

RGBAImageType::Pointer ImageContainer::GetDerivedView(std::string name) {
    RGBAImageType::Pointer slice = g_DerivedSliceDictionary[name];
    // cout << "GetView: " << name << endl << slice << endl;
    return slice;
}

void ImageContainer::GetDerivedViewNames(StringList& names) {
    for (SliceDictionary::iterator iter = g_DerivedSliceDictionary.begin(); iter != g_DerivedSliceDictionary.end(); iter++) {
        names.push_back(iter->first);
    }
}

QPixmap ImageContainer::GetDerivedViewPixmap(std::string name) {
    RGBAImageType::Pointer bitmap = GetDerivedView(name);
    if (bitmap.IsNull()) {
        return QPixmap();
    }
    RGBAImageType::SizeType bitmapSz = bitmap->GetBufferedRegion().GetSize();
    QImage qImg = QImage((unsigned char*) bitmap->GetBufferPointer(),
                         bitmapSz[0], bitmapSz[1], QImage::Format_ARGB32);
    return QPixmap::fromImage(qImg);
}

void ImageContainer::ClearDerivedViews() {
    g_DerivedSliceDictionary.clear();
}



bool ImageContainer::ComputeTransformToPhysicalPointMatrix(VNLMatrix& out) {
    if (!HasImage()) {
        return false;
    }
    SliceType::Pointer sliceImg = GetSlice();
    SliceType::DirectionType dir = sliceImg->GetDirection();
    SliceType::SpacingType spacing = sliceImg->GetSpacing();
    SliceType::PointType origin = sliceImg->GetOrigin();

    // physical space point y = spacing*dir*x - origin
    out.set_size(3,3);
    out.set_identity();
    out.update(dir.GetVnlMatrix());
    for (int i = 0; i < spacing.Size(); i++) {
        out[i][i] *= spacing[i];
    }
    for (int i = 0; i < origin.Size(); i++) {
        out[i][origin.Size()] = origin[i];
    }
    return true;
}

bool ImageContainer::ComputeTransformToIndexMatrix(VNLMatrix& out) {
    VNLMatrix index2world(SDim+1,SDim+1);
    if (!ComputeTransformToPhysicalPointMatrix(index2world)) {
        return false;
    };
    VNLMatrix world2index = vnl_matrix_inverse<double>(index2world);
    out.set_size(3,3);
    out.set_identity();
    out.update(world2index);
    return true;
}

void ImageContainer::TransformToPhysicalPoints(const int n, double *pointsIn, double *pointsOut) {
    if (!HasImage()) {
        return;
    }
    VNLMatrixRef in(n, SDim, pointsIn);
    VNLMatrixRef out(n, SDim, pointsOut);
    vnl_transform_points_2d(m_ImageToWorld, in, out);
}

void ImageContainer::TransformToImagePoints(const int n, double *pointsIn, double *pointsOut) {
    if (!HasImage()) {
        return;
    }
    VNLMatrixRef in(n, SDim, pointsIn);
    VNLMatrixRef out(n, SDim, pointsOut);
    vnl_transform_points_2d(m_WorldToImage, in, out);
}

void ImageContainer::TransformToImageVector(VectorType& src, VectorType& dst) {
    VNLVec3 s;
    for (int i = 0; i < src.Size(); i++) {
        s[i] = src[i];
    }
    s[2] = 0;
    VNLVector t = m_WorldToImage * s;
    for (int i = 0; i < dst.Size(); i++) {
        dst[i] = t[i];
    }
}

// static methods
//
//
RGBAImageType::Pointer ImageContainer::CreateBitmap(SliceType::Pointer slice, int alpha) {
    typedef itk::ScalarToARGBColormapImageFilter<SliceType, RGBAImageType> ScalarToRGBFilter;
    ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
    rgbFilter->SetInput(slice);
    rgbFilter->UseManualScalingOff();
    rgbFilter->UseInputImageExtremaForScalingOn();
    rgbFilter->SetAlphaValue(alpha);
    rgbFilter->Update();
    return rgbFilter->GetOutput();
}


QPixmap ImageContainer::CreatePixmap(RGBAImageType::Pointer bitmap) {
    if (bitmap.IsNull()) {
        return QPixmap();
    }
    RGBAImageType::SizeType bitmapSz = bitmap->GetBufferedRegion().GetSize();
    QImage qImg = QImage((unsigned char*) bitmap->GetBufferPointer(),
                         bitmapSz[0], bitmapSz[1], QImage::Format_ARGB32);
    return QPixmap::fromImage(qImg);
}

SliceType::Pointer ImageContainer::TransformSlice(SliceType::Pointer src, SliceTransformType::Pointer txf, bool useNN) {

    SliceResamplerType::Pointer resampler = SliceResamplerType::New();
    resampler->SetInput(src);
    resampler->UseReferenceImageOn();
    resampler->SetReferenceImage(src);
    resampler->SetTransform(txf);
    if (useNN) {
        NNSliceInterpolatorType::Pointer nn = NNSliceInterpolatorType::New();
        nn->SetInputImage(src);
        resampler->SetInterpolator(nn);
    }
    resampler->Update();
    return resampler->GetOutput();
}

SliceType::Pointer ImageContainer::CreateCheckerBoards(SliceType::Pointer ref, VNLVector &pattern) {
    itkcmds::itkImageIO<SliceType> io;
    SliceType::Pointer blackImage = io.NewImageT(ref);
    blackImage->FillBuffer(0.0);
    SliceType::Pointer whiteImage = io.NewImageT(ref);
    whiteImage->FillBuffer(1.0);
    typedef itk::CheckerBoardImageFilter<SliceType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    FilterType::PatternArrayType checkPatterns;
    for (int k = 0; k < checkPatterns.Size(); k++) {
        checkPatterns[k] = pattern[k];
    }
    filter->SetInput1(blackImage);
    filter->SetInput2(whiteImage); 
    filter->SetCheckerPattern(checkPatterns);
    filter->Update();
    return filter->GetOutput();
}


// compute transform field
////
//void ImageContainer::ComputeTransformedField(DisplacementFieldType::Pointer inField, VNLMatrix& ox, VNLMatrix& oy) {
//    DisplacementFieldType::SizeType sz = inField->GetBufferedRegion().GetSize();
//    ox.set_size(sz[0], sz[1]);
//    oy.set_size(sz[0], sz[1]);
//
//    FieldTransformType::Pointer txf = FieldTransformType::New();
//    txf->SetDisplacementField(inField);
//
//    FieldIteratorType iter(inField, inField->GetBufferedRegion());
//    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
//        FieldIteratorType::IndexType idx = iter.GetIndex();
//        FieldTransformType::InputPointType inPoint;
//        FieldTransformType::OutputPointType outPoint;
//        inField->TransformIndexToPhysicalPoint(idx, inPoint);
//        outPoint = txf->TransformPoint(inPoint);
//        ox[idx[0]][idx[1]] = outPoint[0];
//        oy[idx[0]][idx[1]] = outPoint[1];
//    }
//}


void ImageContainer::WarpGrid(SliceTransformType::Pointer txf, VNLMatrix& sx, VNLMatrix& sy, VNLMatrix& tx, VNLMatrix& ty) {
    const int nRows = sx.rows();
    const int nCols = sx.cols();
    tx.set_size(nRows, nCols);
    ty.set_size(nRows, nCols);

    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            SliceTransformType::InputPointType inPoint;
            inPoint[0] = sx[i][j];
            inPoint[1] = sy[i][j];
            SliceTransformType::OutputPointType outPoint;
            outPoint = txf->TransformPoint(inPoint);
            tx[i][j] = outPoint[0];
            ty[i][j] = outPoint[1];
        }
    }
}


//
//SliceType::Pointer ImageContainer::WarpSlice(SliceType::Pointer src, SliceTransformType::Pointer txf) {
//    typedef itk::TransformToDisplacementFieldSource<DisplacementFieldType> TransformToFieldType;
//    TransformToFieldType::Pointer txfToField = TransformToFieldType::New();
//    txfToField->SetTransform(txf.GetPointer());
//    txfToField->SetOutputDirection(src->GetDirection());
//    txfToField->SetOutputSize(src->GetBufferedRegion().GetSize());
//    txfToField->SetOutputOrigin(src->GetOrigin());
//    txfToField->SetOutputSpacing(src->GetSpacing());
//    txfToField->Update();
//
//    FieldTransformType::Pointer fieldTransform = FieldTransformType::New();
//    fieldTransform->SetDisplacementField(txfToField->GetOutput());
//
//    return TransformSlice(src, fieldTransform.GetPointer());
//}
