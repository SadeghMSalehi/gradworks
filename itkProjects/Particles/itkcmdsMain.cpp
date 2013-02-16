//
//  itkcmdsMain.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#include "itkImageIO.h"
#include "piImageDef.h"
#include "SimpleOpt.h"
#include "piOptions.h"
#include "itkZeroCrossingImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"

using namespace std;
using namespace pi;

Options g_option;

pi::LabelImage::Pointer zeroCrossingFilter(pi::LabelImage::Pointer image) {
    typedef itk::ZeroCrossingImageFilter<LabelImage,LabelImage> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(image);
    filter->Update();
    return filter->GetOutput();
}

void meshToBinary(std::string& in, std::string& out) {

}


void vectorMagnitude(std::string& in, std::string& out) {
    itkcmds::itkImageIO<VectorImage> io1;
    itkcmds::itkImageIO<DoubleImage> io2;
    typedef itk::VectorMagnitudeImageFilter<VectorImage, DoubleImage> F;
    F::Pointer f = F::New();
    f->SetInput(io1.ReadImageT(in.c_str()));
    f->Update();
    io2.WriteImageT(out.c_str(), f->GetOutput());
}

int main(int argc, char* argv[]) {
    enum { OPT_HELP };
    CSimpleOpt::SOption g_rgOptions[] = {
        // ID       TEXT          TYPE
        { OPT_HELP, "-?",     SO_NONE    }, // "-?"
        { OPT_HELP, "--help", SO_NONE    }, // "--help"
        SO_END_OF_OPTIONS                       // END
    };
    g_option.ParseOptions(argc, argv, g_rgOptions);
    pi::StringVector& args = g_option.GetStringVector("args");

    if (args.size() < 3) {
        cout << "usage: itkcmds command file1 file2 ... fileN [options]" << endl;
        cout << "\tcommands:" << endl;
        cout << "\t\tzeroCrossing - apply zero crossing filter" << endl;
        cout << endl;
        return 0;
    }

    LabelImage::Pointer image;
    LabelImage::Pointer imageOutput;
    itkcmds::itkImageIO<pi::LabelImage> io;

    string cmd = args[0];
    if (cmd == "zeroCrossing") {
        image = io.ReadImageT(args[1].c_str());
        imageOutput = zeroCrossingFilter(image);
    } else if (cmd == "mesh2binary") {
        meshToBinary(args[1], args[2]);
    } else if (cmd == "magnitude") {
        vectorMagnitude(args[1], args[2]);
    } else {
        cout << "unknown command: " << cmd << endl;
    }
    if (imageOutput.IsNotNull()) {
        io.WriteImageT(args[2].c_str(), imageOutput);
    }
}