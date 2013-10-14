//
//  simul2main.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/23/13.
//
//

#include <iostream>
#include "QMainWindow"
#include "piOptions.h"
#include "piSimul2.h"
#include "piParticleCore.h"
#include "piImageIO.h"
#include "piPatchCompare.h"

using namespace std;

class MainApps: public QApplication {
public:
    MainApps(int &argc, char* argv[]): QApplication(argc, argv) {
    }

    virtual ~MainApps() {
    }

    virtual bool notify(QObject* rec, QEvent* ev) {
        try {
            return QApplication::notify(rec, ev);
        } catch (std::exception& ex) {
            cout << ex.what() << endl;
        }
        return true;
    }
};


void testUBlasVector() {
    using namespace pi;

    ParticleArray particles(10);
    for (int i = 0; i < particles.size(); i++) {
        particles[i].x[0] = i;
    }
    particles.erase_element(5);
    for (int i = 0; i < particles.size(); i++) {
        cout << particles[i].x[0] << endl;
    }
}


int main(int argc, char* argv[]) {
    using namespace pi;
    CSimpleOpt::SOption options[] ={
        { 1, "--gui", SO_NONE },
        { 2, "--test", SO_NONE },
        { 3, "--makePatch", SO_NONE },
        { 4, "--labelTransferWithPatch", SO_NONE },
        { 5, "-o", SO_REQ_SEP },
        SO_END_OF_OPTIONS
    };

    Options parser;
    StringVector& args = parser.ParseOptions(argc, argv, options);

    if (parser.GetBool("--test")) {
        testUBlasVector();
    } else if (parser.GetBool("--makePatch")) {
        ImageIO<RealImage> io;
        RealImage::Pointer image = io.ReadCastedImage(args[0]);
        PatchCompare patchMaker;
        PatchImage::Pointer patchImage = patchMaker.buildPatchImage(image);
        ImageIO<PatchImage> io2;
        io2.WriteImage(args[1], patchImage);
    } else if (parser.GetBool("--labelTransferWithPatch")) {
        transferLabelsWithPatch(args, parser.GetString("-o"));
    } else if (args.size() == 0 || parser.GetBool("--gui")) {
        MainApps apps(argc, argv);
        Simul2 w;
        w.show();
        return apps.exec();
    } else {

    }
}