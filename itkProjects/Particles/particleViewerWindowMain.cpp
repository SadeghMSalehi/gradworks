#include "particleViewerWindow.h"
#include "itkMacro.h"
#include "piImageDef.h"

#include "iostream"

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
        } catch (itk::ExceptionObject& i) {
            cout << i;
        }
        return true;
    }
};

int main(int argc, char* argv[]) {
    if (__Dim == 3) {
        cout << "running in 3d mode" << endl;
    }
    MainApps apps(argc, argv);
    ParticleViewerWindow w;
    if (argc > 1) {
        w.load(argv[1]);
    } else {
        w.load("/tmpfs/temp.txt");
    }
    w.show();
    return apps.exec();
}