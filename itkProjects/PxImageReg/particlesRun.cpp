#include "piOptions.h"
#include "piParticleRunner.h"
#include "libconfig.h++"

using namespace pi;

// external function defined in piParticleRunner.cpp
namespace pi {
    void executeParticleRunner(Options& parser, StringVector& args);
    void executeDemonsRunner(Options& parser, StringVector& args);
}


static void end() {
    exit(EXIT_SUCCESS);
}

static void particle2mat(StringVector& args) {
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

int main(int argc, char* argv[]) {
    Options opts;
    opts.addOption("--help", SO_NONE);
    opts.addOption("--p2mat", SO_NONE);
	opts.ParseOptions(argc, argv, NULL);
	StringVector& args = opts.GetStringVector("args");

    if (opts.GetBool("--p2mat")) {
        particle2mat(args);
    }
	executeParticleRunner(opts, args);
}
