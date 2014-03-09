//
//  bigview.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/31/13.
//
//

#include "bigview.h"
#include <iostream>
#include <QMainWindow>
#include <QApplication>
#include <QSplashScreen>
#include <QThread>
#include <QTimer>
#include "piOptions.h"
#include "bigViewWindow.h"

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

class Sleeper : public QThread
{
public:
    static void usleep(unsigned long usecs){QThread::usleep(usecs);}
    static void msleep(unsigned long msecs){QThread::msleep(msecs);}
    static void sleep(unsigned long secs){QThread::sleep(secs);}
};

int main(int argc, char* argv[]) {
    using namespace pi;
    CSimpleOpt::SOption options[] ={
        { 1, "--gui", SO_NONE },
        SO_END_OF_OPTIONS
    };

    Options parser;
    StringVector& args = parser.ParseOptions(argc, argv, options);
    MainApps apps(argc, argv);

    QPixmap pixmap(":/Icons/Images/LemonSplash.png");
    QSplashScreen splash(pixmap, Qt::WindowStaysOnTopHint);
    splash.show();

    QTimer timer;
    splash.connect(&timer, SIGNAL(timeout()), SLOT(close()));
    timer.setSingleShot(true);
    timer.start(2000);

    if (args.size() == 0 || parser.GetBool("--gui")) {
        BigViewWindow w;
        w.show();
        return apps.exec();
    } else {
        
    }
}