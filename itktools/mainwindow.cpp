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

void MainWindow::moveSlice() {
    _currentSlice = ui.sliceSlider->sliderPosition();
    drawImage();
}

void MainWindow::drawImage(bool force) {
    if (!force && _currentSlice == _core.CurrentSliceIndex) {
        return;
    }

    _core.SetCurrentSlice(_currentSlice);
    _scene.clear();

    if (_core.SourceSlice.IsNotNull()) {
        if (ui.actionShowSource->isChecked()) {
            _core.SourceSlice->UpdateSlice(_currentSlice, 0xff);
            int* buffer = _core.SourceSlice->GetBitmapBuffer();
            QImage img((unsigned char*) buffer, _core.SourceSlice->GetSize()[0], _core.SourceSlice->GetSize()[1], QImage::Format_ARGB32);
            _scene.addPixmap(QPixmap::fromImage(img));
        }
    }

    if (_core.TargetSlice.IsNotNull()) {
        if (ui.actionShowTarget->isChecked()) {
            _core.TargetSlice->UpdateSlice(_currentSlice, 0xff);
            int* buffer = _core.TargetSlice->GetBitmapBuffer();
            QImage img((unsigned char*) buffer, _core.TargetSlice->GetSize()[0], _core.TargetSlice->GetSize()[1], QImage::Format_ARGB32);
            _scene.addPixmap(QPixmap::fromImage(img));
        }
    }

    if (_core.LabelSlice.IsNotNull()) {
        if (ui.actionShowLabel->isChecked()) {
            _core.LabelSlice->UpdateSlice(_currentSlice, ui.opacityDial->sliderPosition());
            int* buffer = _core.LabelSlice->GetBitmapBuffer();
            QImage img((unsigned char*) buffer, _core.LabelSlice->GetSize()[0], _core.LabelSlice->GetSize()[1], QImage::Format_ARGB32);
            _scene.addPixmap(QPixmap::fromImage(img));
        }
    }

    ui.graphicsView->update();
}

void MainWindow::on_loadTransformButton_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Transform Files (*.txt)"));
    if (fileName != NULL) {
        _core.LoadTransform(fileName.toAscii().data());
        if (ui.applyTransform->isChecked()) {
            _core.ApplyTransform();
            drawImage();
        }

    }
}

void MainWindow::on_actionOpenSource_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName != NULL) {
        _core.LoadImage(fileName.toAscii().data());
        _currentSlice = _core.CurrentSliceIndex;
        ui.sliceSlider->setMaximum(_core.GetMaxSliceIndex());
        ui.sliceSlider->setMinimum(_core.GetMinSliceIndex());
        ui.sliceSlider->setValue(_currentSlice);
        ui.actionOpenTarget->setEnabled(true);
        ui.actionOpenLabel->setEnabled(true);
        drawImage(true);
    }
}

void MainWindow::on_actionOpenTarget_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName != NULL) {
        _core.LoadTarget(fileName.toAscii().data());
        ui.actionShowTarget->setEnabled(true);
       drawImage();
    }
}

void MainWindow::on_actionOpenLabel_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName != NULL) {
        _core.LoadLabelIfGrayImageLoaded(fileName.toAscii().data());
        ui.actionShowLabel->setEnabled(true);
        ui.actionShowLabel->setChecked(true);
        ui.opacityDial->setEnabled(true);
        drawImage();
    }
}

void MainWindow::sayHello(void) {
	QMessageBox::information(this, tr("Apps"), tr("Hello~"), QMessageBox::Ok, QMessageBox::Ok);
}

void MainWindow::exit(void) {
	QApplication::exit();
}
