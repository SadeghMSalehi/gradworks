//
//  simpleSliceViewer.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/2/13.
//
//

#include "simpleSliceViewer.h"
#include "piImageIO.h"
#include "qutils.h"

#include <QActionGroup>
#include <QImage>
#include <QPixmap>

SimpleSliceViewer::SimpleSliceViewer(QWidget* parent): QMainWindow(parent) {
    ui.setupUi(this);
    ui.graphicsView->setScene(&_scene);

    QActionGroup* sliceViewGroup = new QActionGroup(this);
    sliceViewGroup->addAction(ui.actionIJ);
    sliceViewGroup->addAction(ui.actionJK);
    sliceViewGroup->addAction(ui.actionKI);
    sliceViewGroup->setExclusive(true);

}

SimpleSliceViewer::~SimpleSliceViewer() {
}

void SimpleSliceViewer::openImage(QString imageFile) {
    QString file = __fileManager.openFile(QFileManager::Image, this, "Open an Image");
    pi::ImageIO<ImageType> io;
    ImageType::Pointer image = io.ReadCastedImage(file.toUtf8().data());

    _sliceExtractor->SetInput(image);
    updateDisplay();
}


void SimpleSliceViewer::updateDisplay() {
    if (ui.actionIJ->isChecked()) {
        _sliceExtractor->SetSliceAxis(0, 1);
    } else if (ui.actionJK->isChecked()) {
        _sliceExtractor->SetSliceAxis(1, 2);
        _sliceExtractor->SetOutputFlip(itk::MajorFlip);
    } else if (ui.actionKI->isChecked()) {
        _sliceExtractor->SetSliceAxis(2, 0);
        _sliceExtractor->SetOutputFlip(itk::MajorFlip);
    }

    itk::ARGBSlice::Pointer slice = _sliceExtractor->GetOutput();
    const int w = slice->GetBufferedRegion().GetSize(0);
    const int h = slice->GetBufferedRegion().GetSize(1);
    _scene.addPixmap(QPixmap::fromImage(QImage((uchar*) slice->GetBufferPointer(), w, h, QImage::Format_ARGB32)));
}