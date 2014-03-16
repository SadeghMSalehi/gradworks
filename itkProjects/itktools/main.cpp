
#include "mainwindow.h"
#include <QtGui>
#include <QApplication>
#include <QStringList>
#include "itkMyCore.h"
#include "guiregTest.h"

using namespace std;
using namespace itk;

class MainApps: public QApplication {
public:
	MainApps(int &argc, char* argv[]) :
			QApplication(argc, argv) {

	}
	virtual ~MainApps() {
	}
	virtual bool notify(QObject *rec, QEvent *ev) {
		// cDebug() << "Called Application::notify()" << endl;
		try {
			return QApplication::notify(rec, ev);
		} catch (char const *str) {
			cerr << "EXCEPTION: " << str << endl;
			return false;
		} catch (itk::ExceptionObject e) {
			qCritical() << "itkException thrown:" << endl;
			e.Print(cout);
		} catch (std::exception& e) {
			qCritical() << "Exception thrown:" << e.what();
		} catch (...) {
			cerr << "Unknown exception!" << endl;
			abort();
		}
		return true;
	}
};

int main(int argc, char *argv[]) {
    if (argc == 1) {
        MainApps a(argc, argv);
		MainWindow w;
		w.loadDefaults();
		w.show();
		return a.exec();
    }

    itkMyCore core;
    core.ExecuteCommandLine(argc, argv);
	return 0;
}
