//
//  piPlutoWindow.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#include "piPlutoWindow.h"

#include <fstream>
#include <QDesktopWidget>

#include "piImagePatch.h"

const char __contourFile[] = "/NIRAL/work/joohwi/nadia/contour.txt";

namespace pi {
    PlutoWindow::PlutoWindow(QWidget* parent): QMainWindow(parent) {
        _ui.setupUi(this);
        _ui.graphicsView->setScene(&_scene);
        _ui.graphicsMiniView->setScene(&_miniScene);
        _ui.graphicsMiniView->hide();

        _imageGroup = NULL;
        _patchItem = NULL;
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

        const int radius = 3;
        const int size = 7;
        RealImage2::RegionType region;
        region.SetIndex(0, -radius);
        region.SetIndex(1, -radius);
        region.SetSize(0, size);
        region.SetSize(1, size);

        typedef itk::LinearInterpolateImageFunction<RealImage2> InterpolatorType;
        InterpolatorType::Pointer interp = InterpolatorType::New();
        interp->SetInputImage(_images[_images.size() / 2]);
        ImageSamples<RealImage2> samples(1, particles.size(), size*size);
        samples.setSampleRegion(region);
        samples.addInterpolator(interp);
        samples.addParticles(&particles[0]);
        samples.sampleValues();
        RealImage2::Pointer image = samples.getValuesAsImage(size, size);


        if (_patchItem == NULL) {
            _patchItem = new QGraphicsRealImageItem();
            _miniScene.addItem(_patchItem);
            _ui.graphicsMiniView->show();
        }

        _patchItem->setImage<RealImage2>(image);
        _patchItem->setRange(_imageItems[0]->getRange(0), _imageItems[0]->getRange(1));
        _patchItem->refresh();

        _ui.graphicsMiniView->fitInView(0, 0, radius, radius, Qt::KeepAspectRatio);
    }

    void PlutoWindow::on_actionLoad_triggered() {
        ParticleVector contour;
        ifstream is(__contourFile);
        while (is.good()) {
            Particle p;
            is >> p;
            if (is.good()) {
                contour.push_back(p);
            }
        }

        cout << "Contour Length: " << contour.size() << endl;
        const int centerItem = _images.size() / 2.0;
        _interaction.setParticles(_imageItems[centerItem], contour);
    }

    void PlutoWindow::on_actionSave_triggered() {
        ParticleVector& particles = _interaction.getParticles();
        if (particles.size() == 0) {
            return;
        }

        ofstream os(__contourFile);
        for (int i = 0; i < particles.size(); i++) {
            os << particles[i] << endl;
        }
        os.close();
    }

    
}