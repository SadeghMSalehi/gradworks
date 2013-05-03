//
//  plutomain.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#include <QApplication>

#include "plutomain.h"
#include "piPlutoWindow.h"
#include "piOptions.h"

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
        pi::PlutoWindow w;
        w.show();
        return apps.exec();
    } else {
        
    }
}