#include <itkImage.h>
#include <itkAddImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

typedef itk::Image<short,3> IT;
typedef itk::AddImageFilter<IT,IT> AIF;
typedef itk::SubtractImageFilter<IT,IT> SIF;
typedef itk::ImageFileReader<IT> IFR;
typedef itk::ImageFileWriter<IT> IFW;

int main(int argc, char* argv[]) {
    IFR::Pointer r = IFR::New();
    r->SetFileName(argv[1]);
    r->Update();

    IFR::Pointer r2 = IFR::New();
    r2->SetFileName(argv[2]);
    r2->Update();

    AIF::Pointer bf = AIF::New();
    bf->SetInput1(r->GetOutput());
    bf->SetInput2(r2->GetOutput());
    bf->Update();

    IFW::Pointer w = IFW::New();
    w->SetFileName(argv[3]);
    w->SetInput(bf->GetOutput());
    w->UseCompressionOn();
    w->Write();

}
