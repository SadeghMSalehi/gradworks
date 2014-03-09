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
#include "piTestMain.h"
#include "vtkParticleHelper.h"

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
        { 6, "--searchRadius", SO_REQ_SEP },
        { 7, "--kNearest", SO_REQ_SEP },
        { 8, "--demons", SO_NONE },
        { 9, "--boost", SO_NONE },
        { 10, "--vtk", SO_REQ_SEP },
        SO_END_OF_OPTIONS
    };

    Options parser;
    StringVector& args = parser.ParseOptions(argc, argv, options);

    if (parser.GetBool("--test")) {
        testUBlasVector();
    }


    // delegate patch-related process to PatchCompare
    PatchCompare patchMaker;
    if (patchMaker.main(parser, args)) {
        return EXIT_SUCCESS;
    }

    TestMain testMain;
    if (testMain.main(parser, args)) {
        return EXIT_SUCCESS;
    }

    vtkParticleHelper helper;
    if (helper.main(parser, args)) {
        return EXIT_SUCCESS;
    }

    if (args.size() == 0 || parser.GetBool("--gui")) {
        MainApps apps(argc, argv);
        Simul2 w;
        w.show();
        return apps.exec();
    } else {

    }
}
