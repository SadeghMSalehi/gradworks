#include <itkImage.h>
#include <itkEllipseSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkImageFileWriter.h>

typedef itk::Image<short,3> ImageType;
typedef itk::EllipseSpatialObject<3> EllipseType;
typedef itk::SpatialObjectToImageFilter<EllipseType, ImageType> ConverterType;
typedef itk::ImageFileWriter<ImageType> WriterType;

int main(int argc, char* argv[]) {
    EllipseType::Pointer ellipse1 = EllipseType::New();
		EllipseType::ArrayType radius1;
		radius1[0] = atoi(argv[2]);
		radius1[1] = atoi(argv[3]);
		radius1[2] = atoi(argv[4]);

    ellipse1->SetRadius(radius1);
    ellipse1->SetDefaultInsideValue(1);
    ellipse1->SetDefaultOutsideValue(2);
    ellipse1->Update();

    ConverterType::Pointer converter = ConverterType::New();
    converter->SetInput(ellipse1);
    converter->Update();
    converter->GetOutput();

    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(converter->GetOutput());
    writer->SetFileName(argv[1]);
    writer->UseCompressionOn();
    writer->Write();
}
