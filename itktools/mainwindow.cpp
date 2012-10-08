#include "mainwindow.h"
#include "QMessageBox"
#include "QPainter"
#include "QPen"
#include "QFileDialog"
#include "QTextStream"
#include "QtConcurrentRun"
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
    _registrationWatcher = NULL;
    _viewingDir = 2;
}

MainWindow::~MainWindow() {

}

void MainWindow::moveSlice() {
    _core.SetCurrentSlice(_viewingDir, ui.sliceSlider->value());
    drawImage();
}

void MainWindow::drawImage(bool force) {
    _scene.clear();
    if (ui.actionShowSource->isChecked()) {
        _scene.addPixmap(QPixmap::fromImage(_core.ConvertToQImage(Source, _viewingDir)));
    }
    if (ui.actionShowTarget->isChecked()) {
        _scene.addPixmap(QPixmap::fromImage(_core.ConvertToQImage(Target, _viewingDir)));
    }
    if (ui.actionShowLabel->isChecked()) {
        if (ui.applyTransformCheck->isChecked()) {
            _scene.addPixmap(QPixmap::fromImage(_core.ConvertToQImage(TransformedLabel, _viewingDir)));
        } else {
            _scene.addPixmap(QPixmap::fromImage(_core.ConvertToQImage(Label, _viewingDir)));
        }
    }
    ui.graphicsView->update();
}



void MainWindow::on_runRegistrationButton_clicked() {
    ui.runRegistrationButton->setEnabled(false);
    _core.PrepareRegistration();
    _registrationWatcher = new QFutureWatcher<RegistrationMethod::TransformHistoryType>();
    _registrationWatcher->setFuture(QtConcurrent::run(_core, &itkMyCore::RunRegistration));
    connect(_registrationWatcher, SIGNAL(finished()), this, SLOT(on_registrationFinished()));
}

void MainWindow::on_applyTransformCheck_stateChanged(int check) {
    if (check == Qt::Checked) {
        std::string transformParams = _core.ApplyTransform(ui.currentRegistrationStep->value());
        if (transformParams != "") {
            ui.transformParams->clear();
            ui.transformParams->appendPlainText(QString::fromUtf8(transformParams.c_str()));
        }
    } else {
        ui.transformParams->clear();
    }
    drawImage(true);
}

void MainWindow::loadDefaults() {
    const char* sourceFile = "/tmpfs/data/Atlas/p72_atlas_CLE_Parts.nrrd";
    const char* targetFile = "/tmpfs/data/Atlas/00.Parts.nrrd";
    const char* labelFile = "/tmpfs/data/Atlas/p72_atlas_CLE_Parts.nrrd";
//	const char* sourceFile =
//			"/biomed-resimg/work/joohwi/CLE2-Manual/MultiReg/00.T2.nrrd";
//	const char* targetFile =
//			"/biomed-resimg/work/joohwi/CLE2-Manual/MultiReg/19.T2.nrrd";
//	const char* labelFile =
//			"/biomed-resimg/work/joohwi/CLE2-Manual/MultiReg/00.ManualParts.nrrd";

    _core.LoadImage(sourceFile);

    changeSliceView(2);
    ui.actionOpenTarget->setEnabled(true);
    ui.actionOpenLabel->setEnabled(true);

    _core.LoadTarget(targetFile);
    ui.actionShowTarget->setEnabled(true);
    ui.registrationPanel->setEnabled(true);
    ui.optimizerType->setEnabled(true);
    ui.transformType->setEnabled(true);
    ui.runRegistrationButton->setEnabled(true);
    ui.regParams->setEnabled(true);

    _core.LoadLabel(labelFile);
    ui.actionShowLabel->setEnabled(true);
    ui.actionShowLabel->setChecked(true);
    ui.opacityDial->setEnabled(true);

    drawImage();
}

void MainWindow::on_actionOpenSource_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName != NULL) {
        _core.LoadImage(fileName.toAscii().data());
        ui.actionOpenTarget->setEnabled(true);
        ui.actionOpenLabel->setEnabled(true);
        changeSliceView(_viewingDir);
    }
}

void MainWindow::on_actionOpenTarget_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName != NULL) {
        _core.LoadTarget(fileName.toAscii().data());
        drawImage();

        ui.actionShowTarget->setEnabled(true);
        ui.registrationPanel->setEnabled(true);
        ui.optimizerType->setEnabled(true);
        ui.transformType->setEnabled(true);
        ui.runRegistrationButton->setEnabled(true);
        ui.regParams->setEnabled(true);
    }
}

void MainWindow::on_actionOpenLabel_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName != NULL) {
        _core.LoadLabel(fileName.toAscii().data());
        ui.actionShowLabel->setEnabled(true);
        ui.actionShowLabel->setChecked(true);
        ui.opacityDial->setEnabled(true);
        drawImage();
    }
}

void MainWindow::changeSliceView(int dir) {
//    ui.sliceSlider->setMaximum(_core.GetMaxSliceIndex(dir));
//    ui.sliceSlider->setValue(_core.GetCurrentSliceIndex(dir));
    _viewingDir = dir;
    ui.sliceSlider->setValue(_core.GetCurrentSliceIndex()[_viewingDir]);
    ui.sliceSlider->setMaximum(_core.GetMaxSliceIndex()[_viewingDir]);
    ui.sliceSlider->setMinimum(0);
    drawImage();
}

void MainWindow::sayHello(void) {
	QMessageBox::information(this, tr("Apps"), tr("Hello~"), QMessageBox::Ok, QMessageBox::Ok);
}

void MainWindow::exit(void) {
	QApplication::exit();
}

void MainWindow::on_loadTransformButton_clicked() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Transform Files (*.txt)"));
    if (fileName != NULL) {
        _core.LoadTransform(fileName.toAscii().data());
        on_transformAvailable(_core._registrationAlgorithm->GetNumberOfTransforms());

    }
}

void MainWindow::on_saveTransformButton_clicked() {
   	QString fileName = QFileDialog::getSaveFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Text Files (*.txt)"));
    if (fileName != NULL) {
        _core.SaveTransform(fileName.toAscii().data());
    }
}
