#include "piOptions.h"
#include "piParticleRunner.h"
#include "piImageProc.h"
#include "libconfig.h++"

/**
 
 @mainpage Particle-guided Image Registration Software

 This documentation provides the complete set of functions and usage for the use of Particle-guided Image Registration.

 */

using namespace pi;

// external function defined in piParticleRunner.cpp
namespace pi {
    void executeParticleRunner(Options& parser, StringVector& args);
    void executeDemonsRunner(Options& parser, StringVector& args);
    void executeQARunner(Options& parser, StringVector& args);
    void executeRxRunner(Options& parser, StringVector& args);
    void executeLabelFusionRunner(Options& parser, StringVector& args);

    /// @brief Measure the volume overlap ratio
    void executeVolumeOverlaps(Options& parser, StringVector& args);

    /// @brief Compute the entropy image from a list of intensity images.
    void executeEntropyImage(Options& parser, StringVector& args);

    /// @brief Compute the distance map
    void executeComputeDistanceMap(Options& parser, StringVector &args);
}


static void doMerge(Options& opts, StringVector& args);


static void end() {
    exit(EXIT_SUCCESS);
}

static void die() {
    exit(EXIT_FAILURE);
}

static void particle2mat(Options& opts, StringVector& args) {
    if (!opts.GetBool("--p2mat")) {
        return;
    }
    using namespace libconfig;
    Config config;
    config.readFile(args[0].c_str());

    ofstream of(args[1].c_str());
    Setting& list = config.lookup("particles");
    int n = list[0].getLength();
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < list.getLength(); i++) {
            double d = list[i][j];
            of << d << " ";
        }
        of << endl;
    }
    end();
}

static void doSlice(Options& opts, StringVector& args) {
    if (!opts.GetBool("--slice")) {
        return;
    }
    if (args.size() < 4) {
        cout << "--slice dim index imagefile outputfile" << endl;
        die();
    }

    int dim = atoi(args[0].c_str());
    int slice = atoi(args[1].c_str());

    ImageIO<RealImage3> io;
    ImageInfo info;
    RealImage3::Pointer image = io.ReadCastedImage(args[2], info);
    RealImage2Vector sliceImages = SliceVolume(image, dim);

    if (slice < sliceImages.size()) {
        ImageIO<RealImage2> wio;
        wio.WriteCastedImage(args[3], sliceImages[slice], info.componenttype);
    } else {
        cout << "slice index is out of range" << endl;
    }
    end();
}

static void doSeparate(Options& opts, StringVector& args) {
    if (!opts.GetBool("--separate")) {
        return;
    }

    ImageIO<DisplacementFieldType> io;
    DisplacementFieldType::Pointer input = io.ReadImage(args[0]);

    for (int i = 0; i < DisplacementFieldType::PixelType::Length; i++) {
        ImageIO<RealImage> realIO;
        RealImage::Pointer output = realIO.NewImageS<DisplacementFieldType>(input);

        RealImage::PixelType* oBuf = output->GetBufferPointer();
        itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType> iter(input, input->GetBufferedRegion());

        int j = 0;
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter, j++) {
            DisplacementFieldType::PixelType p = iter.Get();
            oBuf[j] = p[i];
        }

        realIO.WriteImage(args[i+1], output);
    }
    end();
}



int main(int argc, char* argv[]) {
    // show which dimension this executable is handling
    cout << argv[0] << " with dimension = " << __Dim << endl;

    
    Options opts;
    opts.addOption("-o", "Specify a filename for an output image", SO_REQ_SEP);
    opts.addOption("--fusion", "label fusion from a config", "`--fusion config-file output-file target-image`", SO_REQ_SEP);
    opts.addOption("--overlap", "Compute the overlap ratio (dice|jaccard). This option can take two or arguments. The first argument is a gold standard, and other arguments are multiple number of label images to be compared.", "--overlap dice output-text ref1 ref2-1 ref2-2 ... ref2-n", SO_REQ_SEP);
    opts.addOption("--p2mat", "point list to matrix", SO_NONE);
    opts.addOption("--slice", "extract a slice from 3d volume", "--slice dim index imagefile outputfile", SO_NONE);
    opts.addOption("--imageMerge", "merge 2D images into a 3d volume (--imageMerge output input1 input2 ...)", SO_REQ_SEP);
    opts.addOption("--qa", "extract a slice with a label map", SO_NONE);
    opts.addOption("--config", "[file] use this config file", SO_REQ_SEP);
    opts.addOption("--demons", "run Demons registration", SO_NONE);
    opts.addOption("--separate", "[input] [x] [y] [z] ... separate vector images into individual image files", SO_NONE);
    opts.addOption("--rx", "registration experiments ", SO_NONE);
    opts.addOption("--dots", "--rx --dots generate a series of gaussian dot images", SO_NONE);
    opts.addOption("--sigma", "sigma value [double]", "--sigma 0.8", SO_REQ_SEP);
    opts.addOption("--entropyImage", "Compute an entropy image from a set of given images", "`--entropyImage -o output.nrrd input1.nrrd input2.nrrd ...`", SO_NONE);
    opts.addOption("--test", "Run in a test mode. The test mode is context sensitive depending on the given argument. For example, if `--entropyImage --test` is given, it will automatically provide a set of input images and produce an output into a specific directory.", SO_NONE);
    opts.addOption("--distanceMap", "Compute the Danielsson's distance map. This will also generate distmag.nrrd, x.nrrd, y.nrrd, and z.nrrd that stores the component of the vector distance map for debugging purpose.", "--distanceMap input output-vector output-magnitude", SO_NONE);
    opts.addOption("--help", "print this message", SO_NONE);

    opts.ParseOptions(argc, argv, NULL);
    StringVector& args = opts.GetStringVector("args");

    if (opts.GetBool("--help") || opts.GetBool("-h")) {
        cout << "## ParticleRun Command Line Options" << endl;
        opts.PrintUsage();
        return 0;
    }


    particle2mat(opts, args);
    doSlice(opts, args);
    doSeparate(opts, args);


    if (opts.GetBool("--qa")) {
        executeQARunner(opts, args);
    } else if (opts.GetString("--imageMerge", "") != "" && args.size() > 0) {
        doMerge(opts, args);
    } else if (opts.GetBool("--demons")) {
        executeDemonsRunner(opts, args);
    } else if (opts.GetBool("--rx")) {
        executeRxRunner(opts, args);
    } else if (opts.GetString("--fusion", "") != "") {
        executeLabelFusionRunner(opts, args);
    } else if (opts.GetBool("--entropyImage")) {
        executeEntropyImage(opts, args);
    } else if (opts.GetString("--overlap") == "dice" || opts.GetString("--overlap") == "jaccard") {
        executeVolumeOverlaps(opts, args);
    } else if (opts.GetBool("--distanceMap")) {
        executeComputeDistanceMap(opts, args);
    } else {
        executeParticleRunner(opts, args);
    }
}


// merge a given list of arguments into the output-image that is given with the parameter --imageMerge
static void doMerge(Options& opts, StringVector& args) {
    string outputFile = opts.GetString("--imageMerge", "");
    if (outputFile == "") {
        return;
    }

    // first, create an empty image with the same x-y dimension with
    // the first argument image but has z-dimension with the number of arguments

    ImageIO<RealImage2> image2IO;
    ImageIO<RealImage3> image3IO;

    // read the first image and create the outputImage buffer
    RealImage2::Pointer inputImage = image2IO.ReadCastedImage(args[0]);
    RealImage3::Pointer outputImage = image3IO.NewImageT(
                inputImage->GetBufferedRegion().GetSize(0), inputImage->GetBufferedRegion().GetSize(1), args.size());

    // prepare the outputBuffer to write
    RealImage3::PixelType* outputBuffer = outputImage->GetBufferPointer();


    int i = 0;
    while (i < args.size()) {
        // prepare the inputImage buffer
        uint nSize = inputImage->GetPixelContainer()->Size();
        RealImage2::PixelType* inputBuffer = inputImage->GetBufferPointer();

        // copy the inputImage into the outputImage
        int nBytes = nSize * sizeof(RealImage3::PixelType);
        memcpy(outputBuffer, inputBuffer, nBytes);

        // forward the outputBuffer
        outputBuffer += nSize;

        // read the next image
        if (++i < args.size()) {
            inputImage = image2IO.ReadImage(args[i]);
        }
    }

    // output file preparation
    string outputFile = opts.GetString("--imageMerge");

    // write the 3d output image
    image3IO.WriteImage(outputFile, outputImage);

    // exit the program
    exit(0);
}
