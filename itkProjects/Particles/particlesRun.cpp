#include "piParticleCore.h"
#include "piParticleBSpline.h"
#include "piParticleSystemSolver.h"
#include "piParticleTools.h"
#include "piParticleTrace.h"
#include "piOptions.h"
#include "piParticleBSpline.h"

using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    CSimpleOpt::SOption specs[] = {
        { 1, "-r", SO_NONE },
        { 2, "-o", SO_REQ_SEP },
        { 3, "--mean", SO_NONE },
        { 4, "--srcidx", SO_REQ_SEP },
        { 7, "--dstidx", SO_REQ_SEP },
        { 5, "--useEnsemble", SO_NONE },
        { 6, "--noTrace", SO_NONE },
        SO_END_OF_OPTIONS
    };
    Options parser;
    parser.ParseOptions(argc, argv, specs);
    StringVector& args = parser.GetStringVector("args");
    string output = parser.GetString("-o", "");

    ParticleSystemSolver solver;
    ParticleSystem& system = solver.m_System;
    Options& options = solver.m_Options;

    if (parser.GetBool("-r", false)) {
        if (args.size() < 1 || output == "") {
            cout << "registration requires [config.txt] -o [outputimage]" << endl;
            return 0;
        }

        // load data
        solver.LoadConfig(args[0].c_str());
        ImageContext& imageCtx = solver.m_ImageContext;

        // bspline resampling
        ParticleBSpline particleTransform;
        particleTransform.SetReferenceImage(imageCtx.GetLabel(0));

        int srcIdx = atoi(parser.GetString("--srcidx", "1").c_str());
        int dstIdx = atoi(parser.GetString("--dstidx", "0").c_str());

        if (parser.GetBool("--mean")) {
            system.ComputeMeanSubject();
            particleTransform.EstimateTransform(system.GetMeanSubject(), system[srcIdx]);
        } else {
            particleTransform.EstimateTransform(system[dstIdx], system[srcIdx]);
        }
        LabelImage::Pointer outputImage = particleTransform.WarpLabel(imageCtx.GetLabel(srcIdx));


        // write image
        itkcmds::itkImageIO<LabelImage> io;
        io.WriteImageT(output.c_str(), outputImage);
    } else {
        if (args.size() < 2) {
        	cout << "registration requires [config.txt] [output.txt]" << endl;
            return 0;
        }

        if (!solver.LoadConfig(args[0].c_str())) {
            return 0;
        }
        
        if (parser.GetBool("--noTrace")) {
            options.Set("PreprocessingTrace:", string(""));
            options.Set("RunTrace:", string(""));
            cout << options << endl;
            cout << "Trace disabled..." << endl;
        }

        solver.Preprocessing();
        
        if (parser.GetBool("--useEnsemble")) {
            options.Set("+ensemble", true);
            cout << "Ensemble term enabled..." << endl;
        }
        solver.Run();

        solver.m_Options = options;
        solver.SaveConfig(args[1].c_str());

        StringVector& markingImages = solver.m_Options.GetStringVector("FinalMarking:");
        if (markingImages.size() > 0) {
            itkcmds::itkImageIO<LabelImage> io;

            for (int i = 0; i < markingImages.size(); i++) {
                LabelImage::Pointer label = solver.m_ImageContext.GetLabel(i);
                LabelImage::Pointer canvas = io.NewImageT(label);
                ParticleArray& data = system[i].m_Particles;
                MarkAtImage<ParticleArray>(data, data.size(), canvas, 1);
                io.WriteImageT(markingImages[i].c_str(), canvas);
            }
        }
    }
    return 0;
}
