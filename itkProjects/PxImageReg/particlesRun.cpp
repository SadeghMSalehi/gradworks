#include "piOptions.h"
#include "piParticleRunner.h"
#include "piImageProc.h"
#include "libconfig.h++"

using namespace pi;

// external function defined in piParticleRunner.cpp
namespace pi {
    void executeParticleRunner(Options& parser, StringVector& args);
    void executeDemonsRunner(Options& parser, StringVector& args);
    void executeQARunner(Options& parser, StringVector& args);
    void executeRxRunner(Options& parser, StringVector& args);
}


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

static void printHelp(Options& opts) {
    StringVector& specs = opts.GetOptionNames();
    for (int i = 0; i < specs.size(); i++) {
        string name = specs[i];
        string help = opts.GetOptionHelp(name);
        if (help == "") {
            continue;
        }
        cout << name << "\t" << help << endl;
    }
}

int main(int argc, char* argv[]) {
    Options opts;
    opts.addOption("--p2mat", "point list to matrix", SO_NONE);
    opts.addOption("--slice", "extract a slice", SO_NONE);
    opts.addOption("--qa", "extract a slice with a label map", SO_NONE);
    opts.addOption("--config", "[file] use this config file", SO_REQ_SEP);
    opts.addOption("--demons", "run Demons registration", SO_NONE);
    opts.addOption("--separate", "[input] [x] [y] [z] ... separate vector images into individual image files", SO_NONE);
    opts.addOption("--rx", "[fixed-image] [moving-image] [output-image] [output-transform]", SO_NONE);
    opts.addOption("--dots", "--rx --dots generate a series of gaussian dot images", SO_NONE);
    opts.addOption("--sigma", "sigma value [double]", SO_REQ_SEP);
    opts.addOption("--help", SO_NONE);

    opts.ParseOptions(argc, argv, NULL);
    StringVector& args = opts.GetStringVector("args");

    if (opts.GetBool("--help") || opts.GetBool("-h")) {
        printHelp(opts);
        return 0;
    }

    particle2mat(opts, args);
    doSlice(opts, args);
    doSeparate(opts, args);

    if (opts.GetBool("--qa")) {
        executeQARunner(opts, args);
    } else if (opts.GetBool("--demons")) {
        executeDemonsRunner(opts, args);
    } else if (opts.GetBool("--rx")) {
        executeRxRunner(opts, args);
    } else {
        executeParticleRunner(opts, args);
    }
}
