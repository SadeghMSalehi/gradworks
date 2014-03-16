#include <itkImage.h>
#include <itkEllipseSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkImageFileWriter.h>

typedef itk::Image<short,3> ImageType;
typedef itk::EllipseSpatialObject<3> EllipseType;
typedef itk::SpatialObjectToImageFilter<EllipseType, ImageType> ConverterType;
typedef itk::ImageFileWriter<ImageType> WriterType;

int main(int argc, char* argv[]) {
    EllipseType::Pointer ellipse = EllipseType::New();
    ellipse->SetRadius(5);
    ellipse->SetDefaultInsideValue(1);
    ellipse->SetDefaultOutsideValue(0);
    ellipse->Update();
    ConverterType::Pointer converter = ConverterType::New();
    converter->SetInput(ellipse);
    converter->Update();
    converter->GetOutput();

    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(converter->GetOutput());
    writer->SetFileName(argv[1]);
    writer->UseCompressionOn();
    writer->Write();
}
