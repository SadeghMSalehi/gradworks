#include "mainwindow.h"
#include <QtGui>
#include <QApplication>
#include <iostream>
#include <exception>


using namespace std;


class MainApps: public QApplication {
public:
    MainApps(int &argc, char* argv[]):
    QApplication(argc, argv) {
    }

    virtual ~MainApps() {
    }

    virtual bool notify(QObject* rec, QEvent* ev) {
        try {
            return QApplication::notify(rec, ev);
        } catch (logic_error& e) {
            cout << e.what() << endl;
        }
        return true;
    }
};

int main(int argc, char* argv[]) {
    MainApps apps(argc, argv);
    MainWindow w;
    w.show();
    return apps.exec();
}
