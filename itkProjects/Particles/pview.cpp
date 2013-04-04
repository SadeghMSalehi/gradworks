//
//  pview.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/14/13.
//
//

#include "pview.h"
#include <iostream>
#include "QMainWindow"
#include "QApplication"
#include "piOptions.h"
#include "pviewAIRWindow.h"
#include "piQTestModule.h"
#include "airCLI.h"

using namespace std;
using namespace pi;

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

int main(int argc, char* argv[]) {
    using namespace pi;
    CSimpleOpt::SOption options[] ={
        { 1, "--gui", SO_NONE },
        { 2, "--isoRG", SO_NONE },
        { 3, "-f", SO_REQ_SEP },
        { 4, "-b", SO_REQ_SEP },
        // image related command
        { 5, "--extractSlice", SO_NONE },
        { 6, "--dir", SO_REQ_SEP },
        { 7, "--range", SO_REQ_SEP },
        { 8, "--pasteSlice", SO_NONE },
        { 9, "--pasteLabel", SO_NONE },
        // non-image command
        { 11, "-test", SO_NONE },
        { 10, "--fitTest", SO_NONE },
        SO_END_OF_OPTIONS
    };

    Options parser;
    StringVector& args = parser.ParseOptions(argc, argv, options);

    if (parser.GetBool("-test")) {
        pi::QTestModule test;
        test.Run(&parser, args);
        return 0;
    }
    if (args.size() == 0 || parser.GetBool("--gui")) {
        MainApps apps(argc, argv);
        AIRWindow w;
        w.show();
        return apps.exec();
    }
    air::CommandLineTools cmdTools;
    cmdTools.Run(&parser, args);
    return 0;
}