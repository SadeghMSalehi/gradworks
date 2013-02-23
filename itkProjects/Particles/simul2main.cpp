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

int main(int argc, char* argv[]) {
    using namespace pi;
    CSimpleOpt::SOption options[] ={
        { 1, "--gui", SO_NONE },
        SO_END_OF_OPTIONS
    };

    Options parser;
    StringVector& args = parser.ParseOptions(argc, argv, options);

    if (args.size() == 0 || parser.GetBool("--gui")) {
        MainApps apps(argc, argv);
        Simul2 w;
        w.show();
        return apps.exec();
    } else {
        
    }
}