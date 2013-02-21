#include "piParticleCore.h"
#include "piParticleBSpline.h"
#include "piParticleSystemSolver.h"
#include "piParticleTools.h"
#include "piParticleTrace.h"
#include "piOptions.h"
#include "piParticleBSpline.h"
#include "piImageProcessing.h"

using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    CSimpleOpt::SOption specs[] = {
        { 9, "--seeConfig", SO_NONE },
        { 1, "-w", SO_NONE },
        { 8, "-d", SO_NONE },
        { 2, "-o", SO_REQ_SEP },
        { 3, "--mean", SO_NONE },
        { 4, "--srcidx", SO_REQ_SEP },
        { 7, "--dstidx", SO_REQ_SEP },
        { 6, "--noTrace", SO_NONE },
        { 10, "--markTrace", SO_NONE },
        { 11, "--srcsubj", SO_REQ_SEP },
        { 12, "--inputimage", SO_REQ_SEP },
        { 13, "--inputlabel", SO_REQ_SEP },
        { 14, "--normalize", SO_NONE },
        { 15, "--magnitude", SO_NONE },
        { 16, "--distancemap", SO_NONE },
        { 17, "--createmask", SO_NONE },
        { 18, "--rescaletoshort", SO_NONE },
        { 19, "--mask", SO_REQ_SEP },
        SO_END_OF_OPTIONS
    };
    Options parser;
    parser.ParseOptions(argc, argv, specs);
    StringVector& args = parser.GetStringVector("args");
    string output = parser.GetString("-o", "");

    ParticleSystemSolver solver;
    ParticleSystem& system = solver.m_System;
    Options& options = solver.m_Options;

    if (parser.GetBool("--seeConfig")) {
        solver.LoadConfig(args[0].c_str());
        cout << "Option Contents:\n\n" << options << endl;
        if (args.size() > 1) {
            solver.SaveConfig(args[1].c_str());
        }
    } else if (parser.GetBool("-w", false)) {
        if (args.size() < 1 || output == "") {
            cout << "warping requires [config.txt] -o [outputimage]" << endl;
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
            cout << "warping from " << srcIdx << " to " << dstIdx << endl;
            particleTransform.EstimateTransform(system[dstIdx], system[srcIdx]);
        }


        string input = parser.GetString("--inputimage", "");
        string label = parser.GetString("--inputlabel", "");

        cout << parser << endl;

        bool doingSomething = false;
        if (label != "") {
			// write image
			itkcmds::itkImageIO<LabelImage> io;
            LabelImage::Pointer outputImage = particleTransform.WarpLabel(io.ReadImageT(label.c_str()));
			io.WriteImageT(output.c_str(), outputImage);
            doingSomething = true;
        }
        if (input != "") {
			itkcmds::itkImageIO<RealImage> io;
        	RealImage::Pointer outputImage = particleTransform.WarpImage(io.ReadImageT(input.c_str()));
        	io.WriteImageT(output.c_str(), outputImage);
            doingSomething = true;
        }
        if (!doingSomething) {
            cout << "-w requires --inputimage or --inputlabel to warp" << endl;
        }
    } else if (parser.GetBool("--markTrace")) {
        if (args.size() < 2) {
        	cout << "--markTrace requires [trace.txt] [reference-image] [output-image]" << endl;
            return 0;
        }
        ifstream in(args[0].c_str());
        ParticleTrace trace;
        trace.Read(in);
        in.close();
        cout << trace << endl;

        int srcIdx = atoi(parser.GetString("--srcidx", "-1").c_str());
        int srcSubj = atoi(parser.GetString("--srcsubj", "-1").c_str());
        
        itkcmds::itkImageIO<LabelImage> io;
        LabelImage::Pointer ref = io.ReadImageT(args[1].c_str());
        LabelImage::Pointer canvas = io.NewImageT(ref);
        for (int i = 0; i < trace.system.size(); i++) {
            if (srcSubj == -1 || srcSubj == i) {
                for (int j = 0; j < trace.system[i].timeSeries.size(); j++) {
                    for (int k = 0; k <= trace.system[i].maxIdx; k++) {
                        if (srcIdx == -1 || srcIdx == k) {
                            Particle& p = trace.system[i].timeSeries[j][k];
                            IntIndex idx;
                            fordim (l) {
                                idx[l] = p.x[l] + 0.5;
                            }
                            (*canvas)[idx] = j;
                        }
                    }
                }
            }
        }
    } else if (parser.GetBool("--normalize")) {
        if (args.size() < 3) {
            cout << "--normalize requires [input-image] [mask-image] [output-image]" << endl;
            return 0;
        }
        itkcmds::itkImageIO<RealImage> iod;
        itkcmds::itkImageIO<LabelImage> iol;
        
        RealImage::Pointer input = iod.ReadImageT(args[0].c_str());
        LabelImage::Pointer label = iol.ReadImageT(args[1].c_str());
        
        ImageProcessing proc;
        RealImage::Pointer output = proc.NormalizeIntensity(input, label);
        
        iod.WriteImageT(args[2].c_str(), output);
    } else if (parser.GetBool("--magnitude")) {
        if (args.size() < 2) {
            cout << "--magnitude requires [input-vector-image] [output-float-image]" << endl;
            return 0;
        }
        itkcmds::itkImageIO<VectorImage> vio;
        itkcmds::itkImageIO<RealImage> rio;

        VectorImage::Pointer input = vio.ReadImageT(args[0].c_str());
        ImageProcessing proc;
        RealImage::Pointer output = proc.ComputeMagnitudeMap(input);
        rio.WriteImageT(args[1].c_str(), output);
    } else if (parser.GetBool("--rescaletoshort")) {
        if (args.size() < 2) {
            cout << "--rescaletoshort requires [input-real-image] [output-short-image] [output-min] [output-max] (--mask mask-input)" << endl;
            return 0;
        }
        itkcmds::itkImageIO<RealImage> rio;
        itkcmds::itkImageIO<LabelImage> lio;


        string maskInput;
        parser.GetStringTo("--mask", maskInput);
        LabelImage::Pointer mask;
        if (maskInput != "") {
            mask = lio.ReadImageT(maskInput.c_str());
        }

        LabelPixel outMax = 10000;
        LabelPixel outMin = 0;
        if (args.size() > 2) {
            outMin = atoi(args[2].c_str());
        }
        if (args.size() > 3) {
            outMax = atoi(args[3].c_str());
        }
        RealImage::Pointer input = rio.ReadImageT(args[0].c_str());
        ImageProcessing proc;
        LabelImage::Pointer output = proc.NormalizeToIntegralType(input, outMin, outMax, mask);
        lio.WriteImageT(args[1].c_str(), output);
    } else if (parser.GetBool("--distancemap")) {
        if (args.size() < 2) {
            cout << "--distancemap requires [input-label-image] [output-vector-image]" << endl;
            return 0;
        }
        itkcmds::itkImageIO<LabelImage> lio;
        itkcmds::itkImageIO<VectorImage> vio;

        LabelImage::Pointer input = lio.ReadImageT(args[0].c_str());
        ImageProcessing proc;
        VectorImage::Pointer output = proc.ComputeDistanceMap(input);
        vio.WriteImageT(args[1].c_str(), output);
    } else if (parser.GetBool("--createmask")) {
        if (args.size() < 2) {
            cout << "--createmask requires [input-label-image] [input-label-image]" << endl;
            return 0;
        }
        itkcmds::itkImageIO<LabelImage> lio;
        LabelImage::Pointer input = lio.ReadImageT(args[0].c_str());
        ImageProcessing proc;
        LabelImage::Pointer output = proc.ThresholdToBinary(input);
        lio.WriteImageT(args[1].c_str(), output);
    } else if (parser.GetBool("--transform")) {
//        string realInput, labelInput;
//        parser.GetStringTo("--inputimage", realInput);
//        parser.GetStringTo("--inputlabel", labelInput);
//        string transform;
//        ImageProcessing proc;
//        if (realInput != "") {
//            itkcmds::itkImageIO<RealImage> io;
//            RealImage::Pointer output = proc.TransformImage<RealImage>(io.ReadImageT(realInput.c_str()));
//            io.WriteImageT(args[0].c_str(), output);
//        }
//        if (labelInput != "") {
//            itkcmds::itkImageIO<LabelImage> io;
//            LabelImage::Pointer output = proc.TransformImage<LabelImage>(io.ReadImageT(labelInput.c_str()));
//            io.WriteImageT(args[0].c_str(), output);
//        }
    } else {
        if (args.size() < 2) {
        	cout << "registration requires [config.txt] [output.txt]" << endl;
            return 0;
        }

        if (!solver.LoadConfig(args[0].c_str())) {
            return 0;
        }
        
        if (parser.GetBool("--noTrace")) {
            options.SetString("PreprocessingTrace:", string(""));
            options.SetString("RunTrace:", string(""));
            cout << "Trace disabled..." << endl;
        }
        
        if (!options.GetBool("no_preprocessing")) {
            solver.Preprocessing();
        } else {
            cout << "Skipping preprocessing; continue with given particles from the input" << endl;
            cout << "the first particle of each subject: " << endl;
            for (int i = 0; i < system.GetNumberOfSubjects(); i++) {
                cout << system[i][0] << endl;
            }
            cout << "continue optimization;" << endl;
        }
        
        solver.Run();
        solver.SaveConfig(args[1].c_str());

        // final point location marking onto image
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

        StringVector& warpedLabels = solver.m_Options.GetStringVector("OutputLabelToFirst:");
        if (warpedLabels.size() > 0) {
            itkcmds::itkImageIO<LabelImage> io;
            StringVector& labelImages = solver.m_Options.GetStringVector("LabelImages:");
            if (warpedLabels.size() != labelImages.size()) {
                cout << "OutputLabelToFirst: and LabeImages: size are different" << endl;
            } else {
                for (int i = 0; i < warpedLabels.size(); i++) {
                    LabelImage::Pointer label = io.ReadImageT(labelImages[i].c_str());
                    // bspline resampling
                    ParticleBSpline particleTransform;
                    particleTransform.SetReferenceImage(label);
                    particleTransform.EstimateTransform(system[0], system[i]);
                    LabelImage::Pointer outputImage = particleTransform.WarpLabel(label);
                    io.WriteImageT(warpedLabels[i].c_str(), outputImage);
                }
            }
        }
    }
    return 0;
}
