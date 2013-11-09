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

int main(int argc, char* argv[]) {
    Options opts;
    opts.addOption("--help", SO_NONE);
    opts.addOption("--p2mat", SO_NONE);
    opts.addOption("--slice", SO_NONE);
    opts.addOption("--qa", SO_NONE);
    opts.addOption("--config", SO_REQ_SEP);

	opts.ParseOptions(argc, argv, NULL);
	StringVector& args = opts.GetStringVector("args");

    particle2mat(opts, args);
    doSlice(opts, args);
    if (opts.GetBool("--qa")) {
        executeQARunner(opts, args);
    } else {
        executeParticleRunner(opts, args);
    }
	executeParticleRunner(opts, args);
}
