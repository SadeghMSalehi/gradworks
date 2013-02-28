//
//  ShapeGenerator.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#include "pwShapeGenerator.h"
#include "piImageDef.h"
#include "itkImageIO.h"

#include "itkEllipseSpatialObject.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkMesh.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkVTKPolyDataWriter.h"
#include "itkZeroCrossingImageFilter.h"
#include "itkAddImageFilter.h"

using namespace std;
using namespace pi;

class Ellipse {
public:
    double center[__Dim];
    double radius[__Dim];

    bool isIn(double x[__Dim]) {
        bool in = false;
        double r = 0;
        fordim(k) {
            r += (x[k]-center[k])*(x[k]-center[k])/(radius[k]*radius[k]);
        }
        in = r <= 1;

        return in;
    }

    void mark(LabelImage::Pointer img) {
        LabelImageIteratorType iter(img, img->GetBufferedRegion());
        iter.GoToBegin();
        while (!iter.IsAtEnd()) {
            IntIndex ix = iter.GetIndex();
            double x[__Dim] = { 0 };
            fordim (k) {
                x[k] = ix[k] + 0.5;
            }
            if (isIn(x)) {
                img->SetPixel(ix, 1);
            }
            ++iter;
        }
    }
};

int main(int argc, char* argv[]) {
    typedef itk::EllipseSpatialObject<__Dim> EllipseType;
    typedef itk::SpatialObjectToImageFilter<EllipseType, LabelImage> SpatialObjectToImageFilterType;
    typedef itk::Mesh<double> MeshType;
    typedef itk::BinaryMask3DMeshSource<LabelImage, MeshType> MeshSourceType;
    typedef itk::ZeroCrossingImageFilter<LabelImage, LabelImage> ZeroCrossingFilterType;
    typedef itk::AddImageFilter<LabelImage> AddFilterType;

    AddFilterType::Pointer addFilter = AddFilterType::New();



    EllipseType::Pointer ellipse = EllipseType::New();
    ellipse->SetDefaultInsideValue(1);
    ellipse->SetDefaultOutsideValue(0);

    EllipseType::ArrayType axes;
    axes[0] = 40;
    axes[1] = 40;
    axes[2] = 40;
    ellipse->SetRadius(axes);

    LabelImage::Pointer firstImage;
    itkcmds::itkImageIO<LabelImage> io;

    EllipseType::TransformType::Pointer transform = EllipseType::TransformType::New();
    transform->SetIdentity();
    TransformType::OutputVectorType translation;
    TransformType::OutputPointType center;
    translation[0] = 90;
    translation[1] = 50;
    translation[2] = 50;
    transform->Translate(translation, false);

    ellipse->SetObjectToParentTransform(transform);
    SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
    RealImage::SizeType size;
    size[0] = 180;
    size[1] = 100;
    size[2] = 100;
    imageFilter->SetSize(size);
    imageFilter->SetInput(ellipse);
    imageFilter->SetUseObjectValue(true);
    imageFilter->SetOutsideValue(0);
    imageFilter->Update();

    LabelImage::Pointer background = imageFilter->GetOutput();

    Ellipse c;
    c.center[0] = 30; c.center[1] = 50, c.center[2] = 50;
    c.radius[0] = c.radius[1] = c.radius[2] = 20;
    for (int i = 30; i <= 150; i++) {
        cout << "marking " << i << endl;
        c.center[0] = i;
//        c.mark(background);
    }

    io.WriteImageT(argv[1], background);
}

/*
int main(int argc, char* argv[]) {
    typedef itk::EllipseSpatialObject<__Dim> EllipseType;
    typedef itk::SpatialObjectToImageFilter<EllipseType, LabelImage> SpatialObjectToImageFilterType;
    typedef itk::Mesh<double> MeshType;
    typedef itk::BinaryMask3DMeshSource<LabelImage, MeshType> MeshSourceType;
    typedef itk::ZeroCrossingImageFilter<LabelImage, LabelImage> ZeroCrossingFilterType;

    SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
    DoubleImage::SizeType size;
    fordim (k) {
        size[k] = 80;
    }
    imageFilter->SetSize(size);

    EllipseType::Pointer ellipse = EllipseType::New();
    ellipse->SetDefaultInsideValue(255);
    ellipse->SetDefaultOutsideValue(0);

    EllipseType::ArrayType axes;
    axes[0] = 20;
    axes[1] = 10;
    axes[2] = 30;
    ellipse->SetRadius(axes);

    EllipseType::TransformType::Pointer transform = EllipseType::TransformType::New();
    transform->SetIdentity();
    TransformType::OutputVectorType translation;
    TransformType::OutputPointType center;
    fordim (k) {
        translation[k] = 40;
    }
    transform->Translate(translation, false);

    ellipse->SetObjectToParentTransform(transform);
    imageFilter->SetInput(ellipse);
    imageFilter->SetUseObjectValue(true);
    imageFilter->SetOutsideValue(0);
    imageFilter->Update();

    itkcmds::itkImageIO<LabelImage> io;
    io.WriteImageT(argv[1], imageFilter->GetOutput());

    ZeroCrossingFilterType::Pointer zeroFilter = ZeroCrossingFilterType::New();
    zeroFilter->SetInput(imageFilter->GetOutput());
    zeroFilter->Update();
    io.WriteImageT(argv[3], zeroFilter->GetOutput());


    MeshSourceType::Pointer meshSource = MeshSourceType::New();
    meshSource->SetInput(imageFilter->GetOutput());
    meshSource->SetObjectValue(255);
    meshSource->Update();
    MeshType::Pointer mesh = meshSource->GetOutput();

    typedef itk::VTKPolyDataWriter<MeshType> MeshWriterType;
    MeshWriterType::Pointer writer = MeshWriterType::New();
    writer->SetInput(mesh);
    writer->SetFileName(argv[2]);
    writer->Update();
    writer->Write();
}
*/