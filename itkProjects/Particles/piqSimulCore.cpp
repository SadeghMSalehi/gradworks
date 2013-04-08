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

        for (int i = 0; i < 2; i++) {
            _imageItem[i] = NULL;
            _labelItem[i] = NULL;
        }
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
            _imageItem[0] = showImage(_scene[0], _imageItem[0], context.GetRealImage(0));
            _label[0] = context.GetLabel(0);
            if (_label[0].IsNotNull()) {
                _labelItem[0] = showLabel(_scene[0], _labelItem[0], _label[0], _imageItem[0]);
            }
            _imageItem[1] = showImage(_scene[1], _imageItem[1], context.GetRealImage(1));
            _label[1] = context.GetLabel(1);
            if (_label[1].IsNotNull()) {
                _labelItem[1] = showLabel(_scene[1], _labelItem[1], _label[1], _imageItem[1]);
            }
        }
    }

    void SimulCore::connectSignals() {
        if (_ui != NULL) {
            connect(_ui->labelOpacity, SIGNAL(sliderMoved(int)), SLOT(labelOpacityChanged(int)));
            connect(_ui->actionRun, SIGNAL(triggered()), SLOT(run()));
        }
    }

    QGraphicsItem* SimulCore::getImageItem(int n) {
        return _imageItem[n];
    }

    QGraphicsScene* SimulCore::scene(int n) {
        return _scene[n];
    }


    void SimulCore::updateParticles() {
        for (int i = 0; i < 2; i++) {
            _particles[i]->setParentItem(_imageItem[i]);
            _particles[i]->SetParticles(&_solver->m_System[i][0], _solver->m_System[i].GetNumberOfPoints());
            _particles[i]->update();
        }
    }

    void SimulCore::updateParticles(int n, pi::ParticleVector& particles) {
        _particles[n]->setParentItem(_imageItem[n]);
        _particles[n]->SetParticles(&particles[0], particles.size());
        _particles[n]->update();
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
        QTransform transform = _ui->graphicsView->transform();
        _ui->zoom->setValue(transform.m11() * 100);
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
        _image[0] = __imageIO.ReadCastedImage(filename.toUtf8().data());
        _imageItem[0] = showImage(_scene[0], _imageItem[0], _image[0]);
    }

    void SimulCore::openImage2(QString filename) {
        _image[1] = __imageIO.ReadCastedImage(filename.toUtf8().data());
        _imageItem[1] = showImage(_scene[1], _imageItem[1], _image[1]);
    }

    void SimulCore::openLabel1(QString filename) {
        if (_imageItem[0] == NULL) {
            return;
        }

        _label[0] = __labelIO.ReadCastedImage(filename.toUtf8().data());
        showLabel(_scene[0], _labelItem[0], _label[0], _imageItem[0]);
    }

    void SimulCore::openLabel2(QString filename) {
        _label[1] = __labelIO.ReadCastedImage(filename.toUtf8().data());
        showLabel(_scene[1], _labelItem[1], _label[1], _imageItem[1]);

    }

    void SimulCore::labelOpacityChanged(int value) {
        for (int i = 0; i < 2; i++) {
            if (_labelItem[i]) {
                _labelItem[i]->setOpacity(value / 255.0);
            }
        }
    }
    
    void SimulCore::run() {
    }
}