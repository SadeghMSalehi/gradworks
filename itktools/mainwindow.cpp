#include "mainwindow.h"
#include "QMessageBox"
#include "QPainter"
#include "QPen"
#include "QFileDialog"
#include "QTextStream"
#include <QDebug>

QTextStream qin(stdin, QIODevice::ReadOnly);
QTextStream qout(stdout, QIODevice::WriteOnly);
QTextStream qerr(stderr, QIODevice::WriteOnly);

using namespace std;

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
	ui.setupUi(this);
	// connect(ui.action_Draw, SIGNAL(triggered(bool)), ui.imageViewer, SLOT(toggleDraw(void)));
	connect(ui.action_Close, SIGNAL(triggered(bool)), this, SLOT(exit(void)));
    ui.graphicsView->setBackgroundBrush(QBrush(Qt::black, Qt::SolidPattern));
    ui.graphicsView->setScene(&_scene);
}

MainWindow::~MainWindow() {

}

void MainWindow::drawImage() {
    if (_core.CurrentSlice.IsNull()) {
        cout << "Current Slice is NULL!!" << endl;
        return;
    }
    int* buffer = _core.CurrentSlice->GetBufferPointer();
    for (int i = 0; i <_core.GrayImageSize[0] * _core.GrayImageSize[1]; i++) {
        buffer[i] = buffer[i] << 8 | buffer[i] << 16 | buffer[i] | 0xff << 24;
    }
    QImage img((unsigned char*) buffer, _core.GrayImageSize[0], _core.GrayImageSize[1], QImage::Format_RGB32);
    _scene.clear();
    _scene.addPixmap(QPixmap::fromImage(img));
    ui.graphicsView->update();
}

void MainWindow::on_action_Open_triggered(bool checked) {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    _core.LoadImage(fileName.toAscii().data());

    drawImage();
}

void MainWindow::sayHello(void) {
	QMessageBox::information(this, tr("Apps"), tr("Hello~"), QMessageBox::Ok, QMessageBox::Ok);
}

void MainWindow::exit(void) {
	QApplication::exit();
}
