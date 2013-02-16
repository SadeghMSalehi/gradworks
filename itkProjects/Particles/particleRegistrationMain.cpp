#include "piParticleCore.h"
#include "piParticleBSpline.h"
#include "piParticleSystemSolver.h"
#include "itkImageIO.h"
#include "iostream"

using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "usage: " << argv[0] << " particle-output.txt subj1_to_subj0.nrrd [0|1] " << endl;
        return 0;
    }
    
    ParticleSystemSolver solver;
    solver.LoadConfig(argv[1]);

    ParticleSystem& sys = solver.m_System;
    ImageContext& imageCtx = solver.m_ImageContext;

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