#include "iostream"
#include "piImageIO.h"
#include "piImageDef.h"

#include <itkResampleImageFilter.h>

using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    ImageIO<RealImage> io;
    RealImage::Pointer img = io.ReadCastedImage(string(argv[1]));

    typedef ImageIO<RealImage>::TransformType TransformType;
    TransformType::Pointer transformBase = io.ReadTransform(argv[2]);

    typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleFilter;
    ResampleFilter::TransformPointerType transform = dynamic_cast<ResampleFilter::TransformType*>(transformBase.GetPointer());
    transform->Print(cout);
    ResampleFilter::Pointer resample = ResampleFilter::New();
    resample->SetInput(img);
    resample->SetTransform(transform);
    resample->SetReferenceImage(img);
    resample->UseReferenceImageOn();
    resample->Update();
    io.WriteImage(string(argv[3]), resample->GetOutput());
}