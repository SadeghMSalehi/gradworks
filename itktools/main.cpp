
#include "mainwindow.h"
#include <QtGui>
#include <QApplication>
#include <QStringList>
#include "itkMyCore.h"

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
	MainApps a(argc, argv);
	QStringList args = a.arguments();
	if (args.size() == 1) {
		MainWindow w;
		w.loadDefaults();
		w.show();
		return a.exec();
	}

	QString cmd = args.at(1);
	if (cmd == "registration") {
		itkMyCore core;
		core.LoadImage(args.at(2).toAscii().data());
		core.LoadTarget(args.at(3).toAscii().data());
		core.LoadLabelIfGrayImageLoaded(args.at(4).toAscii().data());
		core.PrepareRegistration();
		core.RunRegistration();
		core.ApplyLastTransform();
		if (core.InverseLabelSlice.IsNotNull()) {
			LabelType::Pointer label = core.InverseLabelSlice->GetViewImage();
			if (label.IsNotNull()) {
				itkcmds::itkImageIO<LabelType> io;
				io.WriteImageT(args.at(5).toAscii().data(), label);
			}
		}
		if (args.size() > 6) {
			core.WriteLastTransform(args.at(6).toAscii().data());
		}
	}
	return 0;
}
