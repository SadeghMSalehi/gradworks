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
#include "QGraphicsRectWidget.h"

#include "qutils.h"
#include "piImageIO.h"
#include "piParticleCore.h"
#include "piParticleSystemSolver.h"
#include "piPatchTracking.h"
#include "piImageSlice.h"
#include "piImageEntropyComputer.h"

#include "itkARGBColorFunction.h"

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
            _auxImageItem[i] = NULL;
        }
        _solver = new ParticleSystemSolver();
        _miniScene = new QGraphicsScene();
    }

    SimulCore::~SimulCore() {

    }

    // assume setup is called only once
    void SimulCore::setUi(Ui_Simul2D *ui) {
        this->_ui = ui;

        _scene[0] = new QGraphicsScene(this);
        _scene[1] = new QGraphicsScene(this);

        _ui->graphicsView->setScene(_scene[0]);
        _ui->graphicsView2->setScene(_scene[1]);
        _ui->miniView->setScene(_miniScene);

        for (int i = 0; i < 2; i++) {
            _imageItem[i] = new QGraphicsRealImageItem();
            _labelItem[i] = new QGraphicsPixmapItem(_imageItem[i]);
            _scene[i]->addItem(_imageItem[i]);

            _auxImageItem[i] = new QGraphicsRealImageItem();
            _auxImageItem[i]->hide();
            _scene[i]->addItem(_auxImageItem[i]);

            if (i == 0) {
                _trackingItem[i] = new QGraphicsRectWidget();
                _trackingItem[i]->hide();
                _scene[i]->addItem(_trackingItem[i]);
                _trackingWidget = dynamic_cast<QGraphicsRectWidget*>(_trackingItem[i]);
                _trackingRect = new QGraphicsRectItem();
                _scene[1]->addItem(_trackingRect);

                QPen pen(Qt::yellow, 1, Qt::DashLine);
                pen.setCosmetic(true);
                _trackingRect->setPen(pen);
                _trackingRect->setZValue(11);
                _trackingRect->setOpacity(0.8);
                _trackingRect->hide();

            } else {
                _trackingItem[i] = new QGraphicsPolygonItem();

                QPen pen(Qt::green, 1);
                pen.setCosmetic(true);
                dynamic_cast<QGraphicsPolygonItem*>(_trackingItem[i])->setPen(pen);

                _scene[i]->addItem(_trackingItem[1]);
                _trackingItem[i]->hide();
            }

            _patchItem[i] = new QGraphicsRealImageItem();
            _patchItem[i]->hide();
            _miniScene->addItem(_patchItem[i]);
        }


        _ui->actionShowParticles->setChecked(true);
    }

    void SimulCore::setParticleSolver(ParticleSystemSolver* solver) {
        _solver = solver;

        ImageContext& context = _solver->GetImageContext();
        if (context.Count() > 0) {
            _image[0] = context.GetRealImage(0);
            _imageItem[0] = showImage(0, _image[0]);
            _label[0] = context.GetLabel(0);
            if (_label[0].IsNotNull()) {
                _labelItem[0] = showLabel(0, _label[0]);
            }
            _image[1] = context.GetRealImage(1);
            _imageItem[1] = showImage(1, _image[1]);
            _label[1] = context.GetLabel(1);
            if (_label[1].IsNotNull()) {
                _labelItem[1] = showLabel(1, _label[1]);
            }
        }

        for (int i = 0; i < 2; i++) {
            _auxImageItem[i]->hide();
            _particleItem[i].clear();
        }

        connectSignals();
    }

    void SimulCore::connectSignals() {
        if (_ui != NULL) {
            connect(_ui->labelOpacity, SIGNAL(sliderMoved(int)), SLOT(labelOpacityChanged(int)));
            connect(_ui->actionRun, SIGNAL(triggered()), SLOT(run()));
            connect(_ui->actionShowParticles, SIGNAL(toggled(bool)), SLOT(hideParticles(bool)));
            connect(_ui->placeButton, SIGNAL(clicked()), SLOT(startTrackingMode()));
            connect(_ui->trackButton, SIGNAL(clicked()), SLOT(trackPatch()));

            QGraphicsRectWidget* trackingWidget = dynamic_cast<QGraphicsRectWidget*>(_trackingItem[0]);
            connect(trackingWidget, SIGNAL(widgetMoved(QPointF)), SLOT(trackingWidgetMoved(QPointF)));
        }
    }

    QGraphicsItem* SimulCore::getImageItem(int n) {
        return _imageItem[n];
    }

    QGraphicsScene* SimulCore::scene(int n) {
        return _scene[n];
    }

    void SimulCore::hideParticles(bool show) {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < _particleItem[i].size(); j++) {
                if (!show) {
                    _particleItem[i][j]->show();
                } else {
                    _particleItem[i][j]->hide();
                }
            }
        }
    }

    void SimulCore::createParticleItems(int i, int n) {
        if (_particleItem[i].size() != n) {
            // check if need to create particle item
            for (int j = 0; j < _particleItem[i].size(); j++) {
                _scene[i]->removeItem(_particleItem[i][j]);
            }
            _particleItem[i].resize(n);
            // create particle item
            typedef itk::RGBAPixel<unsigned char> RGBA;
            typedef itk::Function::HSVColormapFunction<float, RGBA> HSVFunction;

            HSVFunction::Pointer hsvFunc = HSVFunction::New();
            hsvFunc->SetMinimumInputValue(0);
            hsvFunc->SetMaximumInputValue(n);

            double r = 2;
            for (int j = 0; j < n; j++) {
                RGBA color = hsvFunc->operator()(j);
                _particleItem[i][j] = new QGraphicsEllipseItem(_imageItem[i]);
                _particleItem[i][j]->setRect(-r/2.0, -r/2.0, r, r);
                _particleItem[i][j]->setZValue(10);
                _particleItem[i][j]->setOpacity(1);
                _particleItem[i][j]->setPen(Qt::NoPen);
                _particleItem[i][j]->setBrush(QBrush(qRgb(color[0], color[1], color[2]), Qt::SolidPattern));
                if (_ui->actionShowParticles->isChecked()) {
                    _particleItem[i][j]->hide();
                }
            }
        }
    }

    void SimulCore::updateParticles() {
        for (int i = 0; i < 2; i++) {
            const int n = _solver->m_System[i].size();
            createParticleItems(i, n);
            for (int j = 0; j < n; j++) {
                Particle& p = _solver->m_System[i][j];
                _particleItem[i][j]->setPos(p.x[0], p.x[1]);
            }
        }
    }

    void SimulCore::updateParticles(int i, pi::ParticleVector& particles) {
        const int n = particles.size();
        createParticleItems(i, n);
        for (int j = 0; j < n; j++) {
            Particle& p = particles[j];
            _particleItem[i][j]->setPos(p.x[0], p.x[1]);
        }
    }

    QGraphicsRealImageItem* SimulCore::showImage(int id, RealImage::Pointer image) {
        _imageItem[id]->setImage<RealImage>(image, true);
        _imageItem[id]->setFlip(QGraphicsRealImageItem::UpDown);
        _imageItem[id]->refresh();

        if (_scene[id] == _ui->graphicsView->scene()) {
            _ui->graphicsView->fitInView(_imageItem[id], Qt::KeepAspectRatio);
            _ui->graphicsView->centerOn(_imageItem[id]);
        } else if (_scene[id] == _ui->graphicsView2->scene()) {
            _ui->graphicsView2->fitInView(_imageItem[id], Qt::KeepAspectRatio);
            _ui->graphicsView2->centerOn(_imageItem[id]);
        }
        QTransform transform = _ui->graphicsView->transform();
        _ui->zoom->setValue(transform.m11() * 100);
        return _imageItem[id];
    }

    QGraphicsPixmapItem* SimulCore::showLabel(int id, LabelImage::Pointer labelImage) {
        QRectF rect = _imageItem[id]->boundingRect();
        QImage image((uchar*) labelImage->GetBufferPointer(),
                     rect.width(), rect.height(), QImage::Format_Indexed8);
        image.setColorCount(3);
        image.setColor(0, qRgba(0,0,0,0));
        image.setColor(1, qRgba(255,0,255,255));
        QPixmap pixmap = QPixmap::fromImage(image);

        _labelItem[id]->setPixmap(pixmap);
        _labelItem[id]->setOpacity(_ui->labelOpacity->value() / 255.0);
        _labelItem[id]->setZValue(1);
        _labelItem[id]->setParentItem(_imageItem[id]);
        return _labelItem[id];
    }


    void SimulCore::showAuxImage(int id, RealImage::Pointer image) {
        _auxImageItem[id]->setRange(_imageItem[id]->getRange(0), _imageItem[id]->getRange(1));
        _auxImageItem[id]->setImage<RealImage>(image);

        qreal xPos = _imageItem[id]->boundingRect().width();
        qreal yPos = 0;
        _auxImageItem[id]->setPos(xPos, yPos);
        _auxImageItem[id]->setTransform(_imageItem[id]->transform());
        _auxImageItem[id]->refresh();
        _auxImageItem[id]->show();
    }

    void SimulCore::openImage(int id, QString filename) {
        _image[id] = __imageIO.ReadCastedImage(filename.toUtf8().data());
        _imageItem[id] = showImage(id, _image[id]);
    }

    void SimulCore::openLabel(int id, QString filename) {
        if (_imageItem[id] == NULL) {
            return;
        }

        _label[id] = __labelIO.ReadCastedImage(filename.toUtf8().data());
        showLabel(id, _label[id]);
    }

    void SimulCore::labelOpacityChanged(int value) {
        for (int i = 0; i < 2; i++) {
            if (_labelItem[i]) {
                _labelItem[i]->setOpacity(value / 255.0);
            }
        }
    }

    void SimulCore::startTrackingMode() {
        _trackingWidget->setSize(_ui->rectSize->value());
        for (int i = 0; i < 2; i++) {
            _tracking.setImage(i, _image[i]);
            _trackingItem[i]->show();
        }

        QRectF imageRect = _imageItem[0]->mapRectFromItem(_trackingWidget, _trackingWidget->boundingRect());
        RealImage::RegionType region;
        rectToRegion(imageRect, region);
        
        _tracking.setPatchRegion(region);
        for (int i = 0; i < 2; i++) {
            _patchItem[i]->setTransform(_imageItem[i]->transform());
        }
        trackingWidgetMoved(QPointF());
    }

    void SimulCore::rectToRegion(QRectF& rect, RealImage::RegionType& region) {
        region.SetIndex(0, rect.x());
        region.SetIndex(1, rect.y());
        region.SetSize(0, rect.width());
        region.SetSize(1, rect.height());
    }

    
    void SimulCore::trackingWidgetMoved(QPointF pos) {
        QRectF imageRect = _imageItem[0]->mapRectFromItem(_trackingWidget, _trackingWidget->boundingRect());
        RealImage::RegionType region;
        rectToRegion(imageRect, region);

        _tracking.translatePatchRegion(region);
        _tracking.beginTracking();

        _trackingRect->setRect(_trackingWidget->boundingRect());
        _trackingRect->setPos(_trackingWidget->pos());
        _trackingRect->show();

        int h = 0;
        for (int i = 0; i < 2; i++) {
            // update mini view
            _patchItem[i]->setRange(_imageItem[i]->getRange(0), _imageItem[i]->getRange(1));
            _patchItem[i]->setImage<RealImage>(_tracking.getPatch(i));
            _patchItem[i]->setPos(QPointF(0, h));
            _patchItem[i]->show();

            h += (_patchItem[i]->boundingRect().height() + 1);
        }
        _ui->miniView->fitInView(_miniScene->sceneRect(), Qt::KeepAspectRatio);

        QGraphicsPolygonItem* item = dynamic_cast<QGraphicsPolygonItem*>(_trackingItem[1]);
        if (item != NULL) {
            item->setTransform(_imageItem[1]->transform());
            item->setPolygon(_tracking.getPatchPolygon(1));
            item->show();
        }
    }
    
    void SimulCore::showEntropyMeasurementImage() {
        QRectF imageRect = _imageItem[0]->mapRectFromItem(_trackingWidget, _trackingWidget->boundingRect());
        RealImage::RegionType region;
        rectToRegion(imageRect, region);
        
        ImageEntropyComputer comp;
        RealImage::Pointer image = comp.computeEntropy(_image[0], region, _image[1], _image[1]->GetBufferedRegion());
        
        showAuxImage(1, image);
    }

    void SimulCore::trackPatch() {
        trackingWidgetMoved(QPointF());
    }

    void SimulCore::run() {
    }
}