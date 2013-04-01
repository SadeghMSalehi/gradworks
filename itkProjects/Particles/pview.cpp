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
#include "piImageSlice.h"
#include "airImageAlgorithm.h"

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
        SO_END_OF_OPTIONS
    };

    Options parser;
    StringVector& args = parser.ParseOptions(argc, argv, options);

    if (args.size() == 0 || parser.GetBool("--gui")) {
        MainApps apps(argc, argv);
        AIRWindow w;
        w.show();
        return apps.exec();
    } else if (parser.GetBool("--isoRG")) {
        int fg = parser.GetStringAsInt("-f", 2);
        int bg = parser.GetStringAsInt("-b", 3);
        if (args.size() < 3) {
            cout << "Isolated Connected Filter: -f [fgId=2] -b [bgId=3] [input-image] [input-label] [output-label]" << endl;
            return 0;
        }
        typedef air::ImageAlgorithm<AIRImage,AIRLabel> Algo;
        Algo algo;

        AIRImage::Pointer gray = __airImageIO.ReadCastedImage(args[0]);
        AIRLabel::Pointer label = __airLabelIO.ReadCastedImage(args[1]);

        AIRLabel::Pointer labelOut = algo.ExecuteIsolatedConnectedImageFilter(label, gray, fg, bg);
        __airLabelIO.WriteImage(args[2], labelOut);
    }
    return 0;
}