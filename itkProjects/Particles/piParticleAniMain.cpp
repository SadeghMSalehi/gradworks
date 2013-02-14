//
//  pwParticleAniMain.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#include "piParticleAniMain.h"

#include "particleAni.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkCellTreeLocator.h"
#include "vtkTriangle.h"
#include "vtkGenericCell.h"

#include "QApplication"

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
    MainApps apps(argc, argv);
    AniWindow w;
    w.show();
    return apps.exec();
}
