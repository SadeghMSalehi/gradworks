#include "piOptions.h"
#include "piImageIO.h"
#include "itkImage.h"
#include "PixelMathImageFilter.h"

using namespace std;
using namespace pi;

#ifndef PIXEL_TYPE
#define PIXEL_TYPE float
#endif

typedef itk::Image<PIXEL_TYPE,3> ImageType;
typedef itk::Image<PIXEL_TYPE,2> Image2D;
ImageIO<ImageType> io;
std::vector<ImageType::Pointer> inputImages;

int main(int argc, char* argv[]) {
    CSimpleOpt::SOption specs[] = {
        { 0, "-o", SO_REQ_SEP },
        { 1, "-e", SO_REQ_SEP },
        { 2, "-2", SO_NONE },
        SO_END_OF_OPTIONS
    };

    Options argParser;
    StringVector args = argParser.ParseOptions(argc, argv, specs);
    string eq = argParser.GetString("-e");
    string outputFilename = argParser.GetString("-o");

    if (eq == "" || args.size() == 0 || outputFilename == "") {
        cout << argv[0] << " [ -e equation ] [ -o output-file ] [input-1:A] [input-2:B] ... [input-5:E]" << endl;
        cout << "       in the equation a pixel value from each input is represented by A,B,C,D,E" << endl;
        cout << "       for example, -e 'A+B' will compute the addition of A and B image," << endl;
        cout << "                    -e 'sqrt(A*A+B*B) will compute the magnitude of A and B," << endl;
        cout << "                    -e 'A==1 ? 2:(B==1 ? 3:4) will combine A and B by replacing into 3 and 4 labels, respectively." << endl;
        cout << "       -o will designate the output name" << endl;
        cout << "       all values are internally casted to a given PIXEL_TYPE(float in default) in compilation, " << endl;
        cout << "       and the output type will be casted to the pixel type of the last image." << endl << endl;
        cout << "       Please refer to http://muparser.beltoforion.de/mup_features.html for functions and operators." << endl;
        return 0;
    }

    if (argParser.GetBool("-2")) {
        cout << "Working on 2D images" << endl;
        ImageIO<Image2D> io2;
        ImageInfo lastImageInfo;
        PixelMathImageFilter<Image2D, Image2D>::Pointer pixelFilter = PixelMathImageFilter<Image2D, Image2D>::New();
        pixelFilter->SetEquation(eq);
        for (int i = 0; i < args.size(); i++) {
            pixelFilter->PushBackInput(io2.ReadCastedImage(args[i], lastImageInfo));
        }
        try {
            pixelFilter->Update();
        } catch (itk::ExceptionObject& e) {
            cout << e.what() << endl;
        }
        io2.WriteCastedImage(outputFilename, pixelFilter->GetOutput(), lastImageInfo.componenttype);
    } else {
        ImageInfo lastImageInfo;
        PixelMathImageFilter<ImageType, ImageType>::Pointer pixelFilter = PixelMathImageFilter<ImageType, ImageType>::New();
        pixelFilter->SetEquation(eq);
        for (int i = 0; i < args.size(); i++) {
            pixelFilter->PushBackInput(io.ReadCastedImage(args[i], lastImageInfo));
        }
        try {
            pixelFilter->Update();
        } catch (itk::ExceptionObject& e) {
            cout << e.what() << endl;
        }

        io.WriteCastedImage(outputFilename, pixelFilter->GetOutput(), lastImageInfo.componenttype);
    }
    return 0;
}
