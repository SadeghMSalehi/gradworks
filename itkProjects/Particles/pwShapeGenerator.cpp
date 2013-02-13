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

using namespace pi;

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