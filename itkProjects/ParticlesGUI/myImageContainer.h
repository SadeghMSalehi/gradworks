//
//  ImageViewManager.h
//  laplacePDE
//
//  Created by Joohwi Lee on 11/13/12.
//
//

#ifndef __laplacePDE__ImageViewManager__
#define __laplacePDE__ImageViewManager__

#include <iostream>
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkScalarToARGBColormapImageFilter.h"
#include "itkRGBAPixel.h"
#include "vector"
#include "map"
#include "itkArray.h"
#include "QPixmap"
#include "myEventCallback.h"
#include "vnlCommon.h"
#include "itkPointSet.h"

#define copyArray(x,y) for (int kkk = 0; kkk < SDim; kkk++) x[kkk] = y[kkk]

const int VDim = 3;
const int SDim = 2;

typedef itk::Image<double, VDim> ImageType;
typedef itk::Image<double, SDim> SliceType;
typedef itk::ImageRegionIteratorWithIndex<SliceType> SliceIteratorType;
typedef itk::LinearInterpolateImageFunction<SliceType> SliceInterpolatorType;
typedef itk::RGBAPixel<unsigned char> RGBAPixel;
typedef itk::Image<RGBAPixel, SDim> RGBAImageType;
typedef itk::Image<unsigned short, VDim> LabelType;
typedef itk::ImageRegionIteratorWithIndex<SliceType> SliceIteratorType;
typedef itk::ImageRegionIteratorWithIndex<LabelType> LabelIteratorType;
typedef itk::Image<unsigned short, SDim> LabelSliceType;
typedef itk::ImageRegionIteratorWithIndex<LabelSliceType> LabelSliceIteratorType;
typedef itk::FixedArray<SliceType::Pointer,VDim> SliceTupleType;
typedef itk::FixedArray<LabelSliceType::Pointer,VDim> LabelSliceTupleType;
typedef itk::FixedArray<RGBAImageType::Pointer,VDim> RGBAImageTupleType;
typedef itk::FixedArray<double,4> StatTupleType;
typedef itk::FixedArray<int,VDim> IntTupleType;
typedef itk::CovariantVector<double,2> GradientType;
typedef itk::Image<GradientType,2> GradientImageType;
typedef SliceInterpolatorType::ContinuousIndexType ContinuousIndexType;
typedef itk::ScalarToARGBColormapImageFilter<SliceType, RGBAImageType> SliceToRGBAImageFilterType;
typedef itk::PointSet<int,2> PointSetType;
typedef itk::Vector<double,2> VectorType;
typedef itk::PointSet<VectorType,2> DisplacementFieldPointSetType;
typedef itk::Image<VectorType,2> DisplacementFieldType;


class ImageContainer: public itk::LightObject {
public:
    typedef ImageContainer Self;
    typedef itk::LightObject Superclass;
    typedef itk::SmartPointer<ImageContainer> Pointer;
    typedef itk::SmartPointer<const ImageContainer> ConstPointer;
    typedef std::vector<ImageContainer::Pointer> List;
    typedef std::map<std::string, RGBAImageType::Pointer> SliceDictionary;
    typedef std::pair<std::string, RGBAImageType::Pointer> NameSlicePair;
    typedef std::vector<std::string> StringList;

    itkNewMacro(Self);
    itkTypeMacro(ImageContainer, itk::LightObject);

    itkGetConstMacro(Name, std::string);
    itkGetConstMacro(Image, ImageType::Pointer);
    itkGetConstMacro(Label, LabelType::Pointer);

    // transformation between world and image space
    bool ComputeTransformToPhysicalPointMatrix(VNLMatrix& out);
    bool ComputeTransformToIndexMatrix(VNLMatrix& out);

    void TransformToPhysicalPoints(const int n, double* pointsIn, double* pointsOut);
    void TransformToImagePoints(const int n, double* pointsIn, double* pointsOut);
    void TransformToImageVector(VectorType& src, VectorType& dst);

    void SetAlpha(int alpha);
    void SetLabelAlpha(int alpha);
    void SetSliceIndex(int dim, int idx);

    void LoadImage(const char* filename);
    void LoadLabel(const char* filename);

    void SetImage(ImageType::Pointer image);
    void SetLabel(LabelType::Pointer label);

    bool HasImage() {
        return m_Image.IsNotNull();
    }

    bool HasLabel() {
        return m_Label.IsNotNull();
    }

    QPixmap GetPixmap(int dim);
    QPixmap GetLabelPixmap(int dim);

    QPixmap GetPixmap() {
        return GetPixmap(m_SliceDir);
    }
    QPixmap GetLabelPixmap() {
        return GetLabelPixmap(m_SliceDir);
    }

    IntTupleType GetSize() {
        return m_MaxSliceIndexes;
    }

    IntTupleType GetSliceIndex() {
        return m_SliceIndexes;
    }

    void SetName(std::string name) {
        m_Name = name;
    }

    SliceType::Pointer GetSlice(int dim) {
        if (m_Slices[dim].IsNull()) {
            UpdateSlice(dim);
        }
        return m_Slices[dim];
    }

    LabelSliceType::Pointer GetLabelSlice(int dim) {
        if (m_LabelSlices[dim].IsNull()) {
            UpdateLabelSlice(dim);
        }
        return m_LabelSlices[dim];
    }

    SliceType::Pointer GetSlice() {
        return GetSlice(m_SliceDir);
    }

    LabelSliceType::Pointer GetLabelSlice() {
        return GetLabelSlice(m_SliceDir);
    }

    void SetSliceDir(int dir) {
        m_SliceDir = dir;
    }

    void SetEventCallback(EventCallback* callback) { m_EventCallback = callback; }

    // derived image related methods
    void AddDerivedView(std::string name, RGBAImageType::Pointer slice);
    static RGBAImageType::Pointer GetDerivedView(std::string name);
    static QPixmap GetDerivedViewPixmap(std::string name);
    static void GetDerivedViewNames(StringList& nameList);
    static void ClearDerivedViews();

    // utility methods
    static RGBAImageType::Pointer CreateBitmap(SliceType::Pointer slice, int alpha = 255);

    static int g_CurrentView;
    static int g_CurrentImage;
    static int g_CurrentLabel;

    static void SetCurrentView(int v) { g_CurrentView = v; }
    static void SetCurrentImage(int i) { g_CurrentImage = i; }
    static void SetCurrentLabel(int l) { g_CurrentLabel = l; }
    
    static int GetCurrentView() { return g_CurrentView; }
    static int GetCurrentImage() { return g_CurrentImage; }
    static int GetCurrentLabel() { return g_CurrentLabel; }

protected:
    ImageContainer() : m_SliceDir(2), m_EventCallback(NULL) {
        m_Alpha = 255;
        m_LabelAlpha = 128;
        m_IntensityStats.Fill(0);
        m_LabelStats.Fill(0);
        m_MaxSliceIndexes.Fill(0);
        m_SliceIndexes.Fill(0);
    };
    virtual ~ImageContainer() {};

    void UpdateSlice(int dim);
    void UpdateBitmap(int dim);
    void UpdateLabelSlice(int dim);
    void UpdateLabelBitmap(int dim);
    
private:
    ImageContainer(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented

    int m_SliceDir;
    int m_Alpha, m_LabelAlpha;
    StatTupleType m_IntensityStats, m_LabelStats;
    ImageType::Pointer m_Image;
    LabelType::Pointer m_Label;
    SliceTupleType m_Slices;
    LabelSliceTupleType m_LabelSlices;
    IntTupleType m_SliceIndexes;
    IntTupleType m_MaxSliceIndexes;
    RGBAImageTupleType m_Bitmaps, m_LabelBitmaps;
    std::string m_Name;
    EventCallback* m_EventCallback;

    VNLMatrix m_ImageToWorld;
    VNLMatrix m_WorldToImage;
};

#endif /* defined(__laplacePDE__ImageViewManager__) */
