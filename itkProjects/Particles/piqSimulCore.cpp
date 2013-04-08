//
//  piqSimulCore.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/7/13.
//
//

#include "piqSimulCore.h"
#include "QGraphicsImageItem.h"
#include "QParticlesGraphicsItem.h"

#include "qutils.h"
#include "piImageIO.h"
#include "piParticleCore.h"
#include "piParticleSystemSolver.h"

namespace piq {
    using namespace pi;

    static ImageIO<RealImage> __imageIO;
    static ImageIO<LabelImage> __labelIO;

    SimulCore::SimulCore(QWidget* parent): QObject(parent) {
        _parent = parent;
        _ui = NULL;
        
        _image1Item = NULL;
        _image2Item = NULL;

        _image1show = true;
        _image2show = false;

        _label1item = _label2item = NULL;
        _solver = new ParticleSystemSolver();
    }

    SimulCore::~SimulCore() {

    }

    void SimulCore::setUi(Ui_Simul2D *ui) {
        this->_ui = ui;
        connectSignals();

        _scene[0] = new QGraphicsScene(this);
        _scene[1] = new QGraphicsScene(this);

        _ui->graphicsView->setScene(_scene[0]);
        _ui->graphicsView2->setScene(_scene[1]);

        _particles[0] = new QParticlesGraphicsItem();
        _particles[1] = new QParticlesGraphicsItem();

        for (int i = 0; i < 2; i++) {
            _particles[i]->SetPen(Qt::NoPen);
            _particles[i]->SetBrush(QBrush(Qt::blue, Qt::SolidPattern));
            _particles[i]->setZValue(10);
        }
    }

    void SimulCore::setParticleSolver(ParticleSystemSolver* solver) {
        _solver = solver;

        ImageContext& context = _solver->GetImageContext();
        if (context.Count() > 0) {
            _image1Item = showImage(_scene[0], _image1Item, context.GetRealImage(0));
            _label1 = context.GetLabel(0);
            if (_label1.IsNotNull()) {
                _label1item = showLabel(_scene[0], _label1item, _label1, _image1Item);
            }
            _image2Item = showImage(_scene[1], _image2Item, context.GetRealImage(1));
            _label2 = context.GetLabel(1);
            if (_label2.IsNotNull()) {
                _label2item = showLabel(_scene[1], _label2item, _label2, _image2Item);
            }
        }
    }

    void SimulCore::connectSignals() {
        if (_ui != NULL) {
            connect(_ui->loadImage1, SIGNAL(fileDropped(QString)), SLOT(openImage1(QString)));
            connect(_ui->loadImage2, SIGNAL(fileDropped(QString)), SLOT(openImage2(QString)));
            connect(_ui->loadLabel1, SIGNAL(fileDropped(QString)), SLOT(openLabel1(QString)));
            connect(_ui->loadLabel2, SIGNAL(fileDropped(QString)), SLOT(openLabel2(QString)));
            connect(_ui->labelOpacity, SIGNAL(sliderMoved(int)), SLOT(labelOpacityChanged(int)));
            connect(_ui->actionRun, SIGNAL(triggered()), SLOT(run()));
        }
    }

    QGraphicsItem* SimulCore::getImageItem(int n) {
        return _image1Item;
    }

    QGraphicsScene* SimulCore::scene(int n) {
        return _scene[n];
    }


    void SimulCore::updateParticles() {
        _particles[0]->setParentItem(_image1Item);
        _particles[0]->SetParticles(&_solver->m_System[0][0], _solver->m_System[0].GetNumberOfPoints());
        _particles[1]->setParentItem(_image2Item);
        _particles[1]->SetParticles(&_solver->m_System[1][0], _solver->m_System[1].GetNumberOfPoints());

        for (int i = 0; i < 2; i++) {
            _particles[i]->update();
        }
    }

    SimulCore::QRealImageItem* SimulCore::showImage(QGraphicsScene* scene, QRealImageItem* item, RealImage::Pointer image) {
        if (item) {
            scene->removeItem(item);
        }
        item = new QGraphicsImageItem<RealImage>();
        item->setImage(image);
        item->setFlip(QRealImageItem::UpDown);
        item->refresh();

        scene->addItem(item);

        if (scene == _ui->graphicsView->scene()) {
            _ui->graphicsView->fitInView(item, Qt::KeepAspectRatio);
            _ui->graphicsView->centerOn(item);
        } else if (scene == _ui->graphicsView2->scene()) {
            _ui->graphicsView2->fitInView(item, Qt::KeepAspectRatio);
            _ui->graphicsView2->centerOn(item);
        }
        return item;
    }

    QGraphicsPixmapItem* SimulCore::showLabel(QGraphicsScene* scene, QGraphicsPixmapItem* item, LabelImage::Pointer labelImage, QGraphicsItem* parent) {
        if (item) {
            scene->removeItem(item);
        }
        if (parent == NULL) {
            return NULL;
        }

        QRectF rect = parent->boundingRect();
        QImage image((uchar*) labelImage->GetBufferPointer(),
                     rect.width(), rect.height(), QImage::Format_Indexed8);
        image.setColorCount(3);
        image.setColor(0, qRgba(0,0,0,0));
        image.setColor(1, qRgba(255,0,255,255));
        QPixmap pixmap = QPixmap::fromImage(image);

        item = new QGraphicsPixmapItem();
        item->setPixmap(pixmap);
        item->setOpacity(_ui->labelOpacity->value() / 255.0);
        item->setZValue(1);
        item->setParentItem(parent);
        return item;
    }


    void SimulCore::openImage1(QString filename) {
        _image1 = __imageIO.ReadCastedImage(filename.toUtf8().data());
        _image1Item = showImage(_scene[0], _image1Item, _image1);
    }

    void SimulCore::openImage2(QString filename) {
        if (_image2Item) {
            _scene[1]->removeItem(_image2Item);
        }
        _image2 = __imageIO.ReadCastedImage(filename.toUtf8().data());
        _image2Item = showImage(_scene[1], _image2Item, _image2);
    }

    void SimulCore::openLabel1(QString filename) {
        if (_image1Item == NULL) {
            return;
        }

        _label1 = __labelIO.ReadCastedImage(filename.toUtf8().data());
        showLabel(_scene[0], _label1item, _label1, _image1Item);
    }

    void SimulCore::openLabel2(QString filename) {
        _label2 = __labelIO.ReadCastedImage(filename.toUtf8().data());
    }

    void SimulCore::labelOpacityChanged(int value) {
        if (_label1item) {
            _label1item->setOpacity(value / 255.0);
        }
        if (_label2item) {
            _label2item->setOpacity(value / 255.0);
        }
    }
    
    void SimulCore::run() {
    }
}