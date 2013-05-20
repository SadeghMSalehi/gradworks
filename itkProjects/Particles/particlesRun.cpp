#include "piImageIO.h"
#include "piParticleCore.h"
#include "piParticleBSpline.h"
#include "piParticleSystemSolver.h"
#include "piParticleTools.h"
#include "piParticleTrace.h"
#include "piOptions.h"
#include "piParticleForces.h"
#include "piParticleBSpline.h"
#include "piImageProcessing.h"

using namespace std;
using namespace pi;

int srcIdx = 1;
int dstIdx = 0;

#define PRINT_IDX() cout << "using src: " << srcIdx << ", dst: " << dstIdx << endl


ImageIO<RealImage> realIO;
ImageIO<LabelImage> labelIO;

void zeroCrossing(Options& parser, StringVector& args) {
    if (args.size() < 2) {
        cout << "--zercorssing requires [srcimg] [dstimg]" << endl;
        return;
    }
    string srcImg = args[0];
    string dstImg = args[1];


    ImageProcessing proc;
    LabelImage::Pointer label = labelIO.ReadImage(args[0].c_str());
    LabelImage::Pointer zeroCross = proc.ZeroCrossing(label);
    labelIO.WriteImage(args[1].c_str(), zeroCross);
    
    return;
}

void computeHistogram(Options& parser, StringVector& args) {
    ImageProcessing proc;
    int nbin = 16;
    int rmin = 0;
    int rmax = 10000;
    parser.GetIntTo("--nbin", nbin);
    parser.GetIntTo("--rmin", rmin);
    parser.GetIntTo("--rmax", rmax);
    cout << proc.ComputeHistogramToString(realIO.ReadImage(args[0].c_str()), nbin, rmin, rmax) << endl;
    return;
}


void traceWarp(Options& parser, StringVector& args) {
    if (args.size() < 3) {
        cout << "--traceWarp requires [trace.txt] [srcimg] [dstimg] [--srcidx n] [--interval 1.0]" << endl;
        return;
    }

    string traceIn = args[0];
    string srcImg = args[1];
    string dstImg = args[2];

    ParticleTrace trace;
    ifstream is(traceIn.c_str());
    trace.Read(is);

    

    
    return;
}


LabelImage::Pointer RemoveBoundary(LabelImage::Pointer labelImage, int size) {
    LabelImage::RegionType region = labelImage->GetBufferedRegion();
    LabelImage::IndexType idx = region.GetIndex();
    LabelImage::IndexType upIdx = region.GetUpperIndex();

    cout << "Removing " << size << " border pixels ..." << endl;
    fordim (k) {
        idx[k] += size;
        upIdx[k] -= size;
    }
    region.SetIndex(idx);
    region.SetUpperIndex(upIdx);

    LabelImage::RegionType wholeRegion = labelImage->GetBufferedRegion();
    LabelImageIteratorType iter(labelImage, wholeRegion);

    iter.GoToBegin();
    while (!iter.IsAtEnd()) {
        if (!region.IsInside(iter.GetIndex())) {
            iter.Set(0);
        }
        ++iter;
    }

    return labelImage;
}

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
        { 23, "--markOutput", SO_NONE },
        { 11, "--srcsubj", SO_REQ_SEP },
        { 12, "--inputimage", SO_REQ_SEP },
        { 13, "--inputlabel", SO_REQ_SEP },
        { 14, "--normalize", SO_NONE },
        { 15, "--magnitude", SO_NONE },
        { 16, "--distancemap", SO_NONE },
        { 17, "--createmask", SO_NONE },
        { 18, "--rescaletoshort", SO_NONE },
        { 19, "--mask", SO_REQ_SEP },
        { 20, "--align", SO_NONE },
        { 21, "--warp", SO_NONE },
        { 22, "--norigidalign", SO_NONE },
        { 24, "--onlyrigidalign", SO_NONE },
        { 25, "--showpoints", SO_NONE },
        { 26, "--eval", SO_NONE },
        { 27, "--histo", SO_NONE },
        { 28, "--nbin", SO_NONE },
        { 29, "--rmin", SO_NONE },
        { 30, "--rmax", SO_NONE },
        { 31, "--removeborder", SO_NONE },
        { 32, "--size", SO_REQ_SEP },
        { 33, "--traceWarp", SO_NONE },
        { 34, "--interval", SO_REQ_SEP },
        { 35, "--zerocrossing", SO_NONE },
        { 36, "--meanwarp", SO_NONE },
        SO_END_OF_OPTIONS
    };

    cout << argv[0] << " version compiled at " << __TIMESTAMP__ << endl;
    Options parser;
    parser.ParseOptions(argc, argv, specs);
    StringVector& args = parser.GetStringVector("args");
    string output = parser.GetString("-o", "");

    ParticleSystemSolver solver;
    ParticleSystem& system = solver.m_System;
    Options& options = solver.m_Options;


    srcIdx = atoi(parser.GetString("--srcidx", "1").c_str());
    dstIdx = atoi(parser.GetString("--dstidx", "0").c_str());

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

        // bspline resampling
        ParticleBSpline particleTransform;
        particleTransform.SetReferenceImage(solver.m_System[0].GetLabel());

        int srcIdx = atoi(parser.GetString("--srcidx", "1").c_str());
        int dstIdx = atoi(parser.GetString("--dstidx", "0").c_str());

        if (parser.GetBool("--mean")) {
            system.ComputeXMeanSubject();
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
			ImageIO<LabelImage> io;
            LabelImage::Pointer outputImage = particleTransform.WarpLabel(io.ReadImage(label.c_str()));
			io.WriteImage(output.c_str(), outputImage);
            doingSomething = true;
        }
        if (input != "") {
			ImageIO<RealImage> io;
        	RealImage::Pointer outputImage = particleTransform.WarpImage(io.ReadImage(input.c_str()));
        	io.WriteImage(output.c_str(), outputImage);
            doingSomething = true;
        }
        if (!doingSomething) {
            cout << "-w requires --inputimage or --inputlabel to warp" << endl;
        }
    } else if (parser.GetBool("--warp")) {
        if (args.size() < 2) {
        	cout << "--warp requires [output.txt] [source-image] [warped-output-image]" << endl;
            return 0;
        }

        PRINT_IDX();
        string outputName = args[1];
        string inputImage, inputLabel, refImageName;
        parser.GetStringTo("--inputimage", inputImage);
        parser.GetStringTo("--inputlabel", inputLabel);
        parser.GetStringTo("--reference", refImageName);

        solver.LoadConfig(args[0].c_str());

        
        if (system.size() < 2) {
            cout << "system is not loaded successfully" << endl;
            return 0;
        }
        if (inputImage != "") {
            RealImage::Pointer refImage;
            // warp from srcidx to dstidx
            if (refImageName != "") {
                refImage = realIO.ReadImage(refImageName.c_str());
            }
            RealImage::Pointer output = warp_image<RealImage>(system[dstIdx], system[srcIdx], realIO.ReadImage(inputImage.c_str()), refImage, false, parser.GetBool("--norigidalign"), parser.GetBool("--onlyrigidalign"));
            realIO.WriteImage(outputName.c_str(), output);
        }
        if (inputLabel != "") {
            LabelImage::Pointer refImage;
            // warp from srcidx to dstidx
            if (refImageName != "") {
                refImage = labelIO.ReadImage(refImageName.c_str());
            }
            LabelImage::Pointer output = warp_image<LabelImage>(system[dstIdx], system[srcIdx], labelIO.ReadImage(inputLabel.c_str()), refImage, true, parser.GetBool("--norigidalign"), parser.GetBool("--onlyrigidalign"));
            labelIO.WriteImage(outputName.c_str(), output);
        }
    } else if (parser.GetBool("--meanwarp")) {
        if (args.size() < 2) {
        	cout << "--meanwarp requires [output.txt] --inputimage|--inputlabel [source-image] [warped-output-image]" << endl;
            return 0;
        }

        PRINT_IDX();
        string outputName = args[1];
        string inputImage, inputLabel, refImageName;
        parser.GetStringTo("--inputimage", inputImage);
        parser.GetStringTo("--inputlabel", inputLabel);
        parser.GetStringTo("--reference", refImageName);

        solver.LoadConfig(args[0].c_str());


        ParticleSubject meanSubj = system.ComputeXMeanSubject();
        
        if (system.size() < 2) {
            cout << "system is not loaded successfully" << endl;
            return 0;
        }
        if (inputImage != "") {
            RealImage::Pointer refImage;
            // warp from srcidx to dstidx
            if (refImageName != "") {
                refImage = realIO.ReadImage(refImageName.c_str());
            }
            RealImage::Pointer output = warp_image<RealImage>(meanSubj, system[srcIdx], realIO.ReadImage(inputImage.c_str()), refImage, false, parser.GetBool("--norigidalign"), parser.GetBool("--onlyrigidalign"));
            realIO.WriteImage(outputName.c_str(), output);
        }
        if (inputLabel != "") {
            LabelImage::Pointer refImage;
            // warp from srcidx to dstidx
            if (refImageName != "") {
                refImage = labelIO.ReadImage(refImageName.c_str());
            }
            LabelImage::Pointer output = warp_image<LabelImage>(meanSubj, system[srcIdx], labelIO.ReadImage(inputLabel.c_str()), refImage, true, parser.GetBool("--norigidalign"), parser.GetBool("--onlyrigidalign"));
            labelIO.WriteImage(outputName.c_str(), output);
        }

    } else if (parser.GetBool("--markTrace")) {
        if (args.size() < 2) {
        	cout << "--meanwarp requires [output.txt] [reference-image] [output-image] --srcidx [point-index] --srcsubj [subject-index]" << endl;
            return 0;
        }
        ifstream in(args[0].c_str());
        ParticleTrace trace;
        trace.Read(in);
        in.close();
        cout << trace << endl;

        int srcIdx = atoi(parser.GetString("--srcidx", "-1").c_str());
        int srcSubj = atoi(parser.GetString("--srcsubj", "-1").c_str());
        
        ImageIO<LabelImage> io;
        LabelImage::Pointer ref = io.ReadImage(args[1].c_str());
        LabelImage::Pointer canvas = io.NewImage(ref);
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
    } else if (parser.GetBool("--markOutput")) {
        if (args.size() < 1) {
            cout << "--markOutput requires [output.txt] [output-image]" << endl;
        }
        solver.LoadConfig(args[0].c_str());
        LabelImage::Pointer label = labelIO.NewImage(solver.m_System[srcIdx].GetLabel());
        int nPoints = system[srcIdx].GetNumberOfPoints();
        cout << "Marking " << nPoints << " points ..." << endl;
        MarkAtImage(system[srcIdx], system[srcIdx].GetNumberOfPoints(), label, 1);
        labelIO.WriteImage(args[1].c_str(), label);
    } else if (parser.GetBool("--showpoints")) {
        if (args.size() < 1) {
            cout << "--showpoints requires [output.txt]" << endl;
            return 0;
        }
        solver.LoadConfig(args[0].c_str());
        for (int i = 0; i < system.size(); i++) {
            cout << "Subject " << i << endl;
            for (int j = 0; j < system[i].GetNumberOfPoints(); j++) {
                cout << system[i][j] << endl;
            }
        }

    } else if (parser.GetBool("--normalize")) {
        if (args.size() < 3) {
            cout << "--normalize requires [input-image] [mask-image] [output-image]" << endl;
            return 0;
        }
        ImageIO<RealImage> iod;
        ImageIO<LabelImage> iol;
        
        RealImage::Pointer input = iod.ReadImage(args[0].c_str());
        LabelImage::Pointer label = iol.ReadImage(args[1].c_str());
        
        ImageProcessing proc;
        RealImage::Pointer output = proc.NormalizeIntensity(input, label);
        
        iod.WriteImage(args[2].c_str(), output);
    } else if (parser.GetBool("--magnitude")) {
        if (args.size() < 2) {
            cout << "--magnitude requires [input-vector-image] [output-float-image]" << endl;
            return 0;
        }
        ImageIO<VectorImage> vio;
        ImageIO<RealImage> rio;

        VectorImage::Pointer input = vio.ReadImage(args[0].c_str());
        ImageProcessing proc;
        RealImage::Pointer output = proc.ComputeMagnitudeMap(input);
        rio.WriteImage(args[1].c_str(), output);
    } else if (parser.GetBool("--rescaletoshort")) {
        if (args.size() < 2) {
            cout << "--rescaletoshort requires [input-real-image] [output-short-image] [output-min] [output-max] (--mask mask-input)" << endl;
            return 0;
        }
        ImageIO<RealImage> rio;
        ImageIO<LabelImage> lio;


        string maskInput;
        parser.GetStringTo("--mask", maskInput);
        LabelImage::Pointer mask;
        if (maskInput != "") {
            mask = lio.ReadImage(maskInput.c_str());
        }

        LabelPixel outMax = 10000;
        LabelPixel outMin = 0;
        if (args.size() > 2) {
            outMin = atoi(args[2].c_str());
        }
        if (args.size() > 3) {
            outMax = atoi(args[3].c_str());
        }
        RealImage::Pointer input = rio.ReadImage(args[0].c_str());
        ImageProcessing proc;
        LabelImage::Pointer output = proc.NormalizeToIntegralType(input, outMin, outMax, mask);
        lio.WriteImage(args[1].c_str(), output);
    } else if (parser.GetBool("--distancemap")) {
        if (args.size() < 2) {
            cout << "--distancemap requires [input-label-image] [output-vector-image]" << endl;
            return 0;
        }
        ImageIO<LabelImage> lio;
        ImageIO<VectorImage> vio;

        LabelImage::Pointer input = lio.ReadImage(args[0].c_str());
        ImageProcessing proc;
        VectorImage::Pointer output = proc.ComputeDistanceMap(input);
        vio.WriteImage(args[1].c_str(), output);
    } else if (parser.GetBool("--createmask")) {
        if (args.size() < 2) {
            cout << "--createmask requires [input-label-image] [input-label-image]" << endl;
            return 0;
        }
        ImageIO<LabelImage> lio;
        LabelImage::Pointer input = lio.ReadImage(args[0].c_str());
        ImageProcessing proc;
        LabelImage::Pointer output = proc.ThresholdToBinary(input);
        lio.WriteImage(args[1].c_str(), output);
    } else if (parser.GetBool("--transform")) {
//        string realInput, labelInput;
//        parser.GetStringTo("--inputimage", realInput);
//        parser.GetStringTo("--inputlabel", labelInput);
//        string transform;
//        ImageProcessing proc;
//        if (realInput != "") {
//            ImageIO<RealImage> io;
//            RealImage::Pointer output = proc.TransformImage<RealImage>(io.ReadImage(realInput.c_str()));
//            io.WriteImage(args[0].c_str(), output);
//        }
//        if (labelInput != "") {
//            ImageIO<LabelImage> io;
//            LabelImage::Pointer output = proc.TransformImage<LabelImage>(io.ReadImage(labelInput.c_str()));
//            io.WriteImage(args[0].c_str(), output);
//        }
    } else if (parser.GetBool("--align")) {
        // we align #1 to #0
        if (args.size() > -1) {
            cout << "--align requires [config.txt] [output.vtk]" << endl;
        }
        solver.LoadConfig(args[0].c_str());
        if (system.points() == 0) {
            cout << "no points loaded" << endl;
        }

        system[srcIdx].ComputeAlignment(system[dstIdx]);
        system[srcIdx].AlignmentTransformX2Y();
        cout << system[srcIdx].alignment << endl;
        cout << system[srcIdx].inverseAlignment << endl;
   } else if (parser.GetBool("--eval")) {
        if (args.size() < 2) {
            cout << "--eval requires [label1] [label2]" << endl;
            return 0;
        }
        AtlasSimilarityScore score;
        score.Compute(labelIO.ReadImage(args[0].c_str()), labelIO.ReadImage(args[1].c_str()));
        cout << score << endl;
        return 0; 
    } else if (parser.GetBool("--histo")) {
        computeHistogram(parser, args);
    } else if (parser.GetBool("--traceWarp")) {
        traceWarp(parser, args);
    } else if (parser.GetBool("--zerocrossing")) {
        zeroCrossing(parser, args);
    } else if (parser.GetBool("--removeborder")) {
        if (args.size() < 2) {
        	cout << "--removeborder requires [input-image] [output-image]" << endl;
            return 0;
        }

        LabelImage::Pointer img = labelIO.ReadImage(args[0].c_str());
        RemoveBoundary(img, parser.GetStringAsInt("--size", 0));
        labelIO.WriteImage(args[1].c_str(), img);

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

        cout << "Spreading Particles ..." << endl;
        solver.SpreadParticles();
        cout << "Start Running... " << endl;
        solver.Run();
        solver.SaveConfig(args[1].c_str());

        // final point location marking onto image
        StringVector& markingImages = solver.m_Options.GetStringVector("FinalMarking:");
        if (markingImages.size() > 0) {
            ImageIO<LabelImage> io;

            for (int i = 0; i < markingImages.size(); i++) {
                LabelImage::Pointer label = solver.m_System[i].GetLabel();
                LabelImage::Pointer canvas = io.NewImage(label);
                ParticleArray& data = system[i].m_Particles;
                MarkAtImage<ParticleArray>(data, data.size(), canvas, 1);
                io.WriteImage(markingImages[i].c_str(), canvas);
            }
        }

        StringVector& warpedLabels = solver.m_Options.GetStringVector("OutputLabelToFirst:");
        if (warpedLabels.size() > 0) {
            ImageIO<LabelImage> io;
            StringVector& labelImages = solver.m_Options.GetStringVector("LabelImages:");
            if (warpedLabels.size() != labelImages.size()) {
                cout << "OutputLabelToFirst: and LabeImages: size are different" << endl;
            } else {
                for (int i = 1; i < warpedLabels.size(); i++) {
                    LabelImage::Pointer label = io.ReadImage(labelImages[i].c_str());
                    // bspline resampling
                    LabelImage::Pointer output = warp_image<LabelImage>(system[0], system[i], label, label, true, false, false);
                    io.WriteImage(warpedLabels[i].c_str(), output);
                }
            }
        }
    }
    return 0;
}
