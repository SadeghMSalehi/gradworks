#include "piImageIO.h"
#include "piOptions.h"

using namespace std;
using namespace pi;

#ifndef PIXEL_TYPE
#define PIXEL_TYPE short
#endif

typedef itk::Image<PIXEL_TYPE,3> ImageType;

static ImageIO<ImageType> imageIO;

int main(int argc, char* argv[]) {
    // Option declaration
    CSimpleOpt::SOption specs[] = {
        { 0, "-x", SO_NONE },
        { 1, "-y", SO_NONE },
        { 2, "-z", SO_NONE },
        SO_END_OF_OPTIONS
    };

    // Option parsing
    Options parser;
    StringVector args = parser.ParseOptions(argc, argv, specs);

    if (args.size() < 2) {
        cout << "kflip [-x -y -z] [input-image] [output-image]" << endl;
        return 0;
    }

    // Read an image and create a new image
    ImageType::Pointer img = imageIO.ReadCastedImage(args[0]);
    ImageType::Pointer newImg = imageIO.CopyImage(img);

    // Option processing
    bool xflip = parser.GetBool("-x", false);
    bool yflip = parser.GetBool("-y", false);
    bool zflip = parser.GetBool("-z", false);

    if (xflip) {
        cout << "flipping in x-direction" << endl;
    }
    if (yflip) {
        cout << "flipping in y-direction" << endl;
    }
    if (zflip) {
        cout << "flipping in z-direction" << endl;
    }

    ImageType::SizeType sz = img->GetBufferedRegion().GetSize();

#pragma omp parallel for
    for (int z = 0; z < sz[2]; z++) {
        for (int y = 0; y < sz[1]; y++) {
            for (int x = 0; x < sz[0]; x++) {
                ImageType::IndexType iidx;
                iidx[0] = xflip ? (sz[0] - x) : x;
                iidx[1] = yflip ? (sz[1] - y) : y;
                iidx[2] = zflip ? (sz[2] - z) : z;
                ImageType::IndexType oidx;
                oidx[0] = x;
                oidx[1] = y;
                oidx[2] = z;
                newImg->SetPixel(oidx, img->GetPixel(iidx));
            }
        }
    }

    imageIO.WriteImage(args[1], newImg);
    return 0;
}