#include "itkImageIO.h"
#include "itkContinuousIndex.h"
#include "itkTransformFileWriter.h"
#include "itkTranslationTransform.h"

typedef itk::Image<int,3> ImageType;

int main(int argc, const char* argv[]) {
    if (argc < 4) {
        printf("Usage: itkImageMatch input reference transform_output\n");
        return -1;
    }
    
    const char* inputImage = argv[1];
    const char* referenceImage = argv[2];
    const char* transformOutput = argv[3];


    itkcmds::itkImageIO<ImageType> io;
    ImageType::Pointer in = io.ReadImageT(inputImage);
    ImageType::Pointer ref = io.ReadImageT(referenceImage);

    ImageType::RegionType inRegion = in->GetLargestPossibleRegion();
    ImageType::RegionType refRegion = ref->GetLargestPossibleRegion();

    ImageType::SizeType inSz = inRegion.GetSize();
    ImageType::SizeType refSz = refRegion.GetSize();

    itk::ContinuousIndex<double, ImageType::ImageDimension> inCenter, refCenter;
    for (int i = 0; i < ImageType::ImageDimension; i++) {
        inCenter[i] = inSz[i] / 2.0;
        refCenter[i] = refSz[i] / 2.0;
    }

    ImageType::PointType inCenterPoint, refCenterPoint;
    in->TransformContinuousIndexToPhysicalPoint(inCenter, inCenterPoint);
    ref->TransformContinuousIndexToPhysicalPoint(refCenter, refCenterPoint);
    
    typedef itk::TranslationTransform<double, ImageType::ImageDimension> TransformType;
    TransformType::OutputVectorType translation;
    for (int i = 0; i < ImageType::ImageDimension; i++) {
        translation[i] = refCenterPoint[i] / 2.0 - inCenterPoint[i] / 2.0;
    }
    TransformType::Pointer transform = TransformType::New();
    transform->Translate(translation);

    itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
    transformWriter->SetInput(transform);
    transformWriter->SetFileName(transformOutput);
    transformWriter->Update();
}
