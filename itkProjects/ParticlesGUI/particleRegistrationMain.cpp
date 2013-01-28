#include "myParticleCore.h"
#include "myParticleBSpline.h"
#include "itkImageIO.h"
#include "iostream"

using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    ParticleSystem sys;
    sys.LoadSystem(argv[1], 1);

    ImageContext& imageCtx = sys.GetImageContext();
    if (imageCtx.GetDoubleImageVector().size() < 2) {
        cout << "No images to register..." << endl;
        return 0;
    }

    if (argc == 3) {
        ParticleBSpline particleTransform;
        particleTransform.SetReferenceImage(imageCtx.GetLabel(0));
        particleTransform.EstimateTransform(sys[0], sys[1]);
        DoubleImage::Pointer outputImage = particleTransform.WarpImage(imageCtx.GetDoubleImage(1));

        itkcmds::itkImageIO<DoubleImage> io;
        io.WriteImageT(argv[2], outputImage);
    } else if (argc == 4) {
        int opt = atoi(argv[3]);
        if (opt == 0) {
            ParticleBSpline particleTransform;
            particleTransform.SetReferenceImage(imageCtx.GetLabel(0));
            particleTransform.EstimateTransform(sys[0], sys[1]);
            LabelImage::Pointer outputImage = particleTransform.WarpLabel(imageCtx.GetLabel(1));

            itkcmds::itkImageIO<LabelImage> io;
            io.WriteImageT(argv[2], outputImage);
        }
    }
    return 0;
}