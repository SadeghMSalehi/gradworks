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
#include "qtypedef.h"

using namespace pi;

BigViewWindow::BigViewWindow(QWidget* parent): QMainWindow(parent) {
    ui.setupUi(this);

    setupUi();
    connectSignalSlots();

    initialize();
}

BigViewWindow::~BigViewWindow() {
    delete _images;
}

#pragma mark -
#pragma Private Functions

void BigViewWindow::setupUi() {
    ui.stackedWidget->hide();

    QGLWidget* glWidget = new QGLWidget();
    ui.graphicsView->setViewport(glWidget);
    
    centerToDesktop();
    raise();
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
}

void BigViewWindow::initialize() {
    _images = new AIRDisplayCollection();
    _sliceDirection = IJ;
    
    ui.graphicsView->setThumbsWidth(0);
    ui.graphicsView->setDisplayCollection(_images, false);
}

#pragma mark -
#pragma Public Slots

void BigViewWindow::openFile(QString fileName) {
    AIRImage::Pointer image = __airImageIO.ReadCastedImage(fileName.toStdString());
    if (image.IsNull()) {
        return;
    }
    _images->AddImage(image);
    _images->SetReferenceId(0);
    if (_images->Count() == 1) {
        _images->SetReferenceSlice(_sliceDirection, 0);
    }
    ui.statusbar->showMessage(QString(tr("%1 loaded").arg(fileName)), 10000);
    ui.graphicsView->setVolumeToShow(_images->Count()-1);
    ui.graphicsView->updateDisplay();
    ui.graphicsView->fitToImage(_images->GetReferenceSize(_sliceDirection)/2);
}

#pragma mark -
#pragma Overrided Functions

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
        const QUrl& url = urls[0];
        if (url.scheme() != "file") {
            return;
        }
        event->acceptProposedAction();
        QString filePath = url.path();
        emit fileDropped(filePath);
    }
}
