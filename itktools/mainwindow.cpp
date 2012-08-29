#include "mainwindow.h"
#include "QMessageBox"
#include "QPainter"
#include "QPen"
#include "QFileDialog"
#include <QDebug>


using namespace std;

MainWindow::MainWindow(QWidget *parent) :
		QMainWindow(parent) {
	ui.setupUi(this);
	connect(ui.action_Draw, SIGNAL(triggered(bool)), ui.imageViewer,
			SLOT(toggleDraw(void)));
	connect(ui.action_Close, SIGNAL(triggered(bool)), this, SLOT(exit(void)));
}

MainWindow::~MainWindow() {

}

void MainWindow::on_action_Open_triggered(bool checked) {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
			"/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
	qDebug() << "open " << fileName;
}

void MainWindow::sayHello(void) {
	QMessageBox::information(this, tr("Apps"), tr("Hello~"), QMessageBox::Ok,
			QMessageBox::Ok);
}

void MainWindow::exit(void) {
	QApplication::exit();
}
