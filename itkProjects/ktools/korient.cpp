//
//  korient.cpp
//  ktools
//
//  Created by Joohwi Lee on 3/6/13.
//
//

#ifndef ktools_korient_cpp
#define ktools_korient_cpp


#include "piImageIO.h"
#include "piOptions.h"

#include "itkOrientImageFilter.h"

#include "algorithm"
#include "string"

using namespace std;
using namespace pi;

#ifndef PIXEL_TYPE
#define PIXEL_TYPE short
#endif

typedef itk::Image<PIXEL_TYPE,3> ImageType;

static ImageIO<ImageType> imageIO;


#define __code(x) if(code==#x) return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_##x
static itk::SpatialOrientation::ValidCoordinateOrientationFlags GetCode(string code) {
    __code(RIP);
    __code(LIP);
    __code(RSP);
    __code(LSP);
    __code(RIA);
    __code(LIA);
    __code(RSA);
    __code(LSA);
    __code(IRP);
    __code(ILP);
    __code(SRP);
    __code(SLP);
    __code(IRA);
    __code(ILA);
    __code(SRA);
    __code(SLA);
    __code(RPI);
    __code(LPI);
    __code(RAI);
    __code(LAI);
    __code(RPS);
    __code(LPS);
    __code(RAS);
    __code(LAS);
    __code(PRI);
    __code(PLI);
    __code(ARI);
    __code(ALI);
    __code(PRS);
    __code(PLS);
    __code(ARS);
    __code(ALS);
    __code(IPR);
    __code(SPR);
    __code(IAR);
    __code(SAR);
    __code(IPL);
    __code(SPL);
    __code(IAL);
    __code(SAL);
    __code(PIR);
    __code(PSR);
    __code(AIR);
    __code(ASR);
    __code(PIL);
    __code(PSL);
    __code(AIL);
    __code(ASL);
    return itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID;
}

int main(int argc, char* argv[]) {
    // Option declaration
    CSimpleOpt::SOption specs[] = {
        { 0, "-i", SO_REQ_SEP },
        { 1, "-o", SO_REQ_SEP },
        SO_END_OF_OPTIONS
    };

    // Option parsing
    Options parser;
    StringVector args = parser.ParseOptions(argc, argv, specs);
    if (args.size() < 2) {
        cout << "korient [-i input-orientation (LPI)] [-o output-orientation (RAI|AIR|...)] [input-image] [output-image]" << endl;
        return 0;
    }

    string instrCode = parser.GetString("-i", "LPI");
    string outstrCode = parser.GetString("-o", "LPI");

    std::transform(instrCode.begin(), instrCode.end(), instrCode.begin(), ::toupper);
    std::transform(outstrCode.begin(), outstrCode.end(), outstrCode.begin(), ::toupper);
    
    cout << "converting " << instrCode << " => " << outstrCode << endl;
    
    itk::SpatialOrientation::ValidCoordinateOrientationFlags inputCode, outputCode;
    inputCode = GetCode(instrCode);
    outputCode = GetCode(outstrCode);

    if (inputCode == outputCode) {
        cout << "no orientation change..." << endl;
        return 0;
    }

    if (inputCode == itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID
        || outputCode == itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID) {
        cout << "wrong orientation code..." << endl;
        return 0;
    }

    // Read an image and create a new image
    ImageInfo info;
    ImageType::Pointer img = imageIO.ReadCastedImage(args[0], info);

    typedef itk::OrientImageFilter<ImageType, ImageType> OrientFilter;
    OrientFilter::Pointer filter = OrientFilter::New();
    filter->SetInput(img);
    filter->SetGivenCoordinateOrientation(inputCode);
    filter->SetDesiredCoordinateOrientation(outputCode);
    ImageType::Pointer output = filter->GetOutput();
    
    imageIO.WriteCastedImage(args[1], output, info.componenttype);
}

#endif
