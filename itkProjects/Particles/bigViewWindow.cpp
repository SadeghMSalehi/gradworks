//
//  bigViewWindow.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/31/13.
//
//

#include "bigViewWindow.h"
#include "piImageDef.h"
#include "piImageSlice.h"
#include "piImageIO.h"
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QDesktopWidget>
#include <QGLWidget>
#include <QShortcut>
#include <QDoubleSpinBox>
#include "qtypedef.h"
#include "qutils.h"

using namespace pi;

BigViewWindow::BigViewWindow(QWidget* parent): QMainWindow(parent) {
    ui.setupUi(this);

    setupUi();
    setupShortcuts();
    connectSignalSlots();

    initialize();
}

BigViewWindow::~BigViewWindow() {
    delete _images;
}

#pragma mark -
#pragma mark Private Functions

void BigViewWindow::setupUi() {
    setWindowFlags( windowFlags() | Qt::WindowMaximizeButtonHint ); // doesn't remove maximise button
    ui.stackedWidget->hide();

    ui.toolBar->addWidget(ui.fileList);


    _lowIntensitySpinBox = new QDoubleSpinBox(ui.toolBar);
    _highIntensitySpinBox = new QDoubleSpinBox(ui.toolBar);

    ui.toolBar->addWidget(_lowIntensitySpinBox);
    ui.toolBar->addWidget(ui.intensitySlider);
    ui.toolBar->addWidget(_highIntensitySpinBox);

    QGLWidget* glWidget = new QGLWidget();
    ui.graphicsView->setViewport(glWidget);
    ui.graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);
    ui.graphicsView->setInteractive(true);
    
    centerToDesktop();
    raise();
}


void BigViewWindow::setupShortcuts() {
    QShortcut* zoomIn = new QShortcut(QKeySequence(tr("w")), this);
    zoomIn->setContext(Qt::ApplicationShortcut);
    connect(zoomIn, SIGNAL(activated()), this, SLOT(zoomIn()));

    QShortcut* zoomOut = new QShortcut(QKeySequence(tr("s")), this);
    zoomOut->setContext(Qt::ApplicationShortcut);
    connect(zoomOut, SIGNAL(activated()), this, SLOT(zoomOut()));
    
}

void BigViewWindow::centerToDesktop() {
    QDesktopWidget *desktop = QApplication::desktop();

    int screenWidth, width;
    int screenHeight, height;
    int x, y;
    QSize windowSize;

    screenWidth = desktop->width();
    screenHeight = desktop->height();

    windowSize = size();
    width = windowSize.width();
    height = windowSize.height();

    x = (screenWidth - width) / 2;
    y = (screenHeight - height) / 2;
    y -= 50;

    move ( x, y );
    setFixedSize(windowSize.width(), windowSize.height());
}

void BigViewWindow::connectSignalSlots() {
    connect(ui.actionShowMenu, SIGNAL(triggered(bool)), ui.stackedWidget, SLOT(setVisible(bool)));
    connect(this, SIGNAL(fileDropped(QString)), this, SLOT(openFile(QString)));
    connect(ui.graphicsView, SIGNAL(fileDropped(QString)), this, SLOT(openFile(QString)));
    connect(ui.fileList, SIGNAL(currentIndexChanged(int)), this, SLOT(volumeSelected(int)));
    connect(ui.intensitySlider, SIGNAL(sliderMoved(int)), this, SLOT(changeIntensity()));
    connect(ui.intensitySlider, SIGNAL(realLowValueChanged(double)), _lowIntensitySpinBox, SLOT(setValue(double)));
    connect(ui.intensitySlider, SIGNAL(realHighValueChanged(double)), _highIntensitySpinBox, SLOT(setValue(double)));
    connect(_lowIntensitySpinBox, SIGNAL(valueChanged(double)), ui.intensitySlider, SLOT(setRealLowValue(double)));
    connect(_highIntensitySpinBox, SIGNAL(valueChanged(double)), ui.intensitySlider, SLOT(setRealHighValue(double)));
    connect(this, SIGNAL(multipleFileDropeed(QList<QString>)), this, SLOT(openFiles(QList<QString>)));
    connect(ui.actionLoad, SIGNAL(triggered()), this, SLOT(openFile()));
}

void BigViewWindow::initialize() {
    _images = new AIRDisplayCollection();
    _sliceDirection = IJ;
    
    ui.graphicsView->setThumbsWidth(0);
    ui.graphicsView->setDisplayCollection(_images, false);
}

#pragma mark -
#pragma mark Public Slots

void BigViewWindow::openFile(QString fileName) {
    if (fileName.isEmpty()) {
        fileName = __fileManager.openFile(0, this, "Choose an image file");
        if (fileName.isEmpty()) {
            return;
        }
    }
    AIRImage::Pointer image = __airImageIO.ReadCastedImage(fileName.toStdString());
    if (image.IsNull()) {
        return;
    }
    _images->AddImage(image, fileName.toStdString());
    _images->SetReferenceId(0);
    if (_images->Count() == 1) {
        _images->SetReferenceSlice(_sliceDirection, 0);
    }
    ui.statusbar->showMessage(QString(tr("%1 loaded").arg(fileName)), 10000);
    ui.graphicsView->setVolumeToShow(_images->Count()-1);
    ui.graphicsView->updateDisplay();
    ui.graphicsView->fitToImage(_images->GetReferenceSize(_sliceDirection)/2);
    
    ui.fileList->addItem(fileName);
}

void BigViewWindow::openFiles(QList<QString> files) {
    for (int i = 0; i < files.size(); i++) {
        QString fileName = files[i];
        if (!fileName.isEmpty()) {
            AIRImage::Pointer image = __airImageIO.ReadCastedImage(fileName.toStdString());
            if (image.IsNull()) {
                return;
            }
            _images->AddImage(image, fileName.toStdString());
            _images->SetReferenceId(0);
            if (_images->Count() == 1) {
                _images->SetReferenceSlice(_sliceDirection, 0);
            }
            ui.statusbar->showMessage(QString(tr("%1 loaded").arg(fileName)), 10000);
            ui.graphicsView->setVolumeToShow(_images->Count()-1);
            ui.graphicsView->fitToImage(_images->GetReferenceSize(_sliceDirection)/2);
            ui.fileList->addItem(fileName);
        }
    }
    ui.graphicsView->updateDisplay();
}

void BigViewWindow::zoomIn() {
    ui.graphicsView->scale(1.1,1.1);
}

void BigViewWindow::zoomOut() {
    ui.graphicsView->scale(1/1.1,1/1.1);
}

void BigViewWindow::volumeSelected(int i) {
    AIRImageDisplay* disp = _images->at(i);
    _lowIntensitySpinBox->setMinimum(disp->histogram.dataMin);
    _lowIntensitySpinBox->setMaximum(disp->histogram.dataMax);
    _highIntensitySpinBox->setMinimum(disp->histogram.dataMin);
    _highIntensitySpinBox->setMaximum(disp->histogram.dataMax);
    _lowIntensitySpinBox->setValue(disp->histogram.rangeMin);
    _highIntensitySpinBox->setValue(disp->histogram.rangeMax);
    ui.intensitySlider->setRealMax(disp->histogram.dataMax);
    ui.intensitySlider->setRealMin(disp->histogram.dataMin);
    ui.intensitySlider->setRealLowValue(disp->histogram.rangeMin);
    ui.intensitySlider->setRealHighValue(disp->histogram.rangeMax);
    ui.graphicsView->moveToVolume(i);
}

void BigViewWindow::changeIntensity() {
    int idx = ui.fileList->currentIndex();
    AIRImageDisplay* disp = _images->at(idx);
    disp->histogram.rangeMin = ui.intensitySlider->realLowValue();
    disp->histogram.rangeMax = ui.intensitySlider->realHighValue();
    ui.graphicsView->updateDisplay();
}

#pragma mark -
#pragma mark Overrided Functions

void BigViewWindow::dragEnterEvent(QDragEnterEvent *event) {
    if (event->mimeData()->hasUrls()) {
        UrlList urls = event->mimeData()->urls();
        UrlList::ConstIterator iter;
        for (iter = urls.constBegin(); iter != urls.constEnd(); iter++) {
            const QUrl& url = (*iter);
            if (url.scheme() != "file") {
                event->setAccepted(false);
                return;
            }
        }
        event->setAccepted(true);
        event->acceptProposedAction();
    }
}

void BigViewWindow::dropEvent(QDropEvent *event) {
    QList<QString> fileNames;
    if (event->mimeData()->hasUrls()) {
        UrlList urls = event->mimeData()->urls();
        if (urls.size() > 1) {
            QList<QString> files;
            for (int i = 0; i < urls.size(); i++) {
                const QUrl& url = urls[i];
                if (url.scheme() != "file") {
                    return;
                }
                event->acceptProposedAction();
                QString filePath = url.path();
                files.append(filePath);
            }
            emit multipleFileDropeed(files);
        } else {
            const QUrl& url = urls[0];
            if (url.scheme() != "file") {
                return;
            }
            event->acceptProposedAction();
            QString filePath = url.path();
            emit fileDropped(filePath);
        }
    }
}