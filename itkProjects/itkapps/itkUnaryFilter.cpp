#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkReassignLabelFilter.h"

typedef itk::Image<short,3> ImageType;
typedef itk::ImageFileReader<ImageType> Reader;
typedef itk::ImageFileWriter<ImageType> Writer;
typedef itk::ReassignLabelFilter<ImageType,ImageType> ReassignFilter;

int main(int argc, char* argv[]) {
    Reader::Pointer r = Reader::New();
    Writer::Pointer w = Writer::New();

    r->SetFileName(argv[1]);
    r->Update();

    ReassignFilter::Pointer f = ReassignFilter::New();
    static_cast<ReassignFilter::Superclass::FunctorType>(f->GetFunctor()).SetTargetLabel(atoi(argv[3]));
    static_cast<ReassignFilter::Superclass::FunctorType>(f->GetFunctor()).SetAssignLabel(atoi(argv[4]));
    f->SetInput(r->GetOutput());
    f->Update();

    w->SetFileName(argv[2]);
    w->UseCompressionOn();
    w->SetInput(f->GetOutput());
    w->Write();
}
