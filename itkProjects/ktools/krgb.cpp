//
//  krgb.cpp
//  ktools
//
//  Created by Joohwi Lee on 3/9/13.
//
//

#include "krgb.h"

#include "piImageDef.h"
#include "piImageIO.h"
#include "piOptions.h"
#include "piImageProcessing.h"
#include "itkUnaryFunctorImageFilter.h"


typedef itk::Image<itk::Vector<float,3>, 3> DeformFieldImageType;
typedef itk::Image<itk::RGBPixel<u_char>, 3> RGBImageType;

using namespace std;
using namespace pi;

class Vector2RGB {
public:
    RGBImageType::PixelType operator()(DeformFieldImageType::PixelType v) {
        RGBImageType::PixelType rgb;
        v.Normalize();
        v *= 255;
        rgb.Set(v[0], v[1], v[2]);
        return rgb;
    }
};

int main(int argc, char* argv[]) {
    __noverbose = 0;
    CSimpleOpt::SOption specs[] = {
        { 0, "-o", SO_REQ_SEP },
        SO_END_OF_OPTIONS
    };

    Options argParser;
    StringVector args = argParser.ParseOptions(argc, argv, specs);

    if (args.size() < 1) {
        cout << argv[0] << " [input-vector] [output-rgb]" << endl;
        return 1;
    }

    ImageIO<DeformFieldImageType> io;
    DeformFieldImageType::Pointer img = io.ReadImage(args[0]);

//    itk::UnaryFunctorImageFilter<DeformFieldImageType,RGBImageType> CastFilter;
    typedef itk::UnaryFunctorImageFilter<DeformFieldImageType, RGBImageType, Vector2RGB> CastFilter;
    CastFilter::Pointer caster = CastFilter::New();
    caster->SetInput(img);
    caster->Update();
    RGBImageType::Pointer rgbImg = caster->GetOutput();

    rgbImg->Print(cout);
    io.WriteImageS<RGBImageType>(args[1], rgbImg);

    return 0;
}