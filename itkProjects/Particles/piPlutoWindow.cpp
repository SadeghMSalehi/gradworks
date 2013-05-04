//
//  piPlutoWindow.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#include "piPlutoWindow.h"


#include <QDesktopWidget>

#include "piImagePatch.h"

namespace pi {
    PlutoWindow::PlutoWindow(QWidget* parent): QMainWindow(parent) {
        _ui.setupUi(this);
        _ui.graphicsView->setScene(&_scene);

        _imageGroup = NULL;
    }
    
    PlutoWindow::~PlutoWindow() {
        
    }
    
    void PlutoWindow::centerToDesktop() {
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
        
        move(x, y);
        resize(windowSize.width(), windowSize.height());
    }
    
    void PlutoWindow::connectSignals() {
        
    }

    void PlutoWindow::flipImages() {
        for (int i = 0; i < _imageItems.size(); i++) {
            _imageItems[i]->setFlip(QGraphicsRealImageItem::UpDown);
        }
    }

    void PlutoWindow::on_actionOpen_triggered() {
        if (_imageGroup) {
            _scene.removeItem(_imageGroup);
        }

        // create new group as GraphicsScene removed it
        _imageGroup = new QGraphicsItemGroup();
        _scene.addItem(_imageGroup);

        RealImage3::Pointer volume = __real3IO.ReadCastedImage("/NIRAL/work/joohwi/nadia/E04Aligned2/C31_E04_T1.nii.gz");

        ImageHistogram<RealImage3> imageHisto;
        imageHisto.SetImage(volume);

        _images = __realImageTools.sliceVolume(volume, 1);
        _imageItems.resize(_images.size());

        int wPos = 0;
        for (int i = 0; i < _images.size(); i++) {
            _imageItems[i] = new QGraphicsRealImageItem();
            _imageItems[i]->setImage<RealImage2>(_images[i], false);
            _imageItems[i]->setRange(imageHisto.rangeMin, imageHisto.rangeMax);
            _imageItems[i]->refresh();
            _imageItems[i]->setPos(wPos, 0);
            _imageItems[i]->setInteraction(&_interaction);
            _imageItems[i]->setFlip(QGraphicsRealImageItem::UpDown);

            wPos += _imageItems[i]->boundingRect().width() + 3;

//            _imageGroup->addToGroup(_imageItems[i]);
            _scene.addItem(_imageItems[i]);

            QGraphicsTextItem* label = new QGraphicsTextItem(_imageItems[i]);
            label->setPlainText(QString("%1").arg(i));
            label->setPos(3,3);
            label->setFlag(QGraphicsItem::ItemIgnoresTransformations, true);
            label->setZValue(10);
            label->setDefaultTextColor(Qt::white);
        }

//        const int height = _images[0]->GetBufferedRegion().GetSize(1);
//        QTransform flipTransform;
//        flipTransform.setMatrix(1, 0, 0, 0, -1, 0, 0, height, 1);
//        _imageGroup->setTransform(flipTransform);

        _ui.graphicsView->fitInView(_imageItems[_images.size() / 2], Qt::KeepAspectRatio);
    }

    void PlutoWindow::on_actionReset_triggered() {
        _interaction.reset();
    }

    void PlutoWindow::on_actionStart_triggered() {
        ParticleVector& particles = _interaction.getParticles();
        if (particles.size() == 0) {
            return;
        }


        RealImage2::RegionType region;
        region.SetIndex(0, -3);
        region.SetIndex(1, -3);
        region.SetSize(0, 7);
        region.SetSize(1, 7);

        typedef itk::LinearInterpolateImageFunction<RealImage2> InterpolatorType;
        InterpolatorType::Pointer interp = InterpolatorType::New();
        interp->SetInputImage(_images[_images.size() / 2]);
        ImageSamples<RealImage2> samples(1, particles.size(), 49);
        samples.setSampleRegion(region);
        samples.addInterpolator(interp);
        samples.addParticles(&particles[0]);
        samples.sampleValues();


    }
}