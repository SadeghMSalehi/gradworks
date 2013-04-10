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
            _imageItem[i] = new QGraphicsImageItem<RealImage>();
            _labelItem[i] = new QGraphicsPixmapItem(_imageItem[i]);
            _scene[i]->addItem(_imageItem[i]);

            _auxImageItem[i] = new QGraphicsImageItem<RealImage>();
            _auxImageItem[i]->hide();
            _scene[i]->addItem(_auxImageItem[i]);

            _rectItem[i] = new QGraphicsRectWidget();
            _rectItem[i]->hide();
            _scene[i]->addItem(_rectItem[i]);

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
            connect(_rectItem[0], SIGNAL(widgetMoved(QPointF)), SLOT(trackingWidgetMoved(QPointF)));
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

    SimulCore::QRealImageItem* SimulCore::showImage(int id, RealImage::Pointer image) {
        _imageItem[id]->setImage(image);
        _imageItem[id]->setFlip(QRealImageItem::UpDown);
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
        _auxImageItem[id]->setImage(image);
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
        for (int i = 0; i < 2; i++) {
            _rectItem[i]->setSize(_ui->rectSize->value());
            _rectItem[i]->show();

            samplePixels(i);
        }
        _ui->miniView->fitInView(_miniScene->sceneRect(), Qt::KeepAspectRatio);
    }

    void SimulCore::trackingWidgetMoved(QPointF pos) {
        _rectItem[1]->setTransform(_rectItem[0]->transform());
        _rectItem[1]->setPos(_rectItem[0]->pos());

        samplePixels(0);
        samplePixels(1);
    }

    void SimulCore::trackPatch() {
        PatchTracking tracking;
        for (int i = 0; i < 2; i++) {
            tracking.setImage(i, _image[i]);

            // rectangle in the image coordinate
            QRectF rect = _imageItem[i]->mapRectFromItem(_rectItem[i], _rectItem[i]->boundingRect());
            RealImage::RegionType region;
            region.SetIndex(0, rect.x());
            region.SetIndex(1, rect.y());
            region.SetSize(0, rect.width());
            region.SetSize(1, rect.height());
            tracking.setInitialRegion(i, region);
        }

        tracking.beginTracking();
        QRectF imagePos;
        imagePos.setCoords(tracking.getFinalIndex(0)[0], tracking.getFinalIndex(0)[1],
                           tracking.getFinalIndex(1)[0], tracking.getFinalIndex(1)[1]);
        QRectF scenePos = _imageItem[1]->mapRectToScene(imagePos);
        _rectItem[1]->setPos(scenePos.topLeft());

        _patchItem[1]->setImage(tracking.getPatch(1),
                    _imageItem[1]->histogram().dataMin, _imageItem[1]->histogram().dataMax);
        _patchItem[1]->refresh();
    }

    void SimulCore::samplePixels(int i) {
        QRectF rect = _imageItem[i]->mapRectFromItem(_rectItem[i], _rectItem[i]->boundingRect());

        RealImage::RegionType sampleRegion;
        sampleRegion.SetSize(0, rect.width());
        sampleRegion.SetSize(1, rect.height());
        sampleRegion.SetIndex(0, rect.x());
        sampleRegion.SetIndex(1, rect.y());

        if (_patch[i].IsNull() || _patch[i]->GetBufferedRegion().GetSize() != sampleRegion.GetSize()) {
            _patch[i] = __imageIO.NewImageT(sampleRegion.GetSize());
        }
        _patch[i]->SetBufferedRegion(sampleRegion);
        _patch[i]->SetRequestedRegion(sampleRegion);

        typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleFilter;
        static ResampleFilter::Pointer resample = ResampleFilter::New();
        resample->SetOutputStartIndex(sampleRegion.GetIndex());
        resample->SetSize(sampleRegion.GetSize());
        resample->SetInput(_image[i]);
        resample->GraftOutput(_patch[i]);
        resample->Update();

        int h = 0;
        for (int j = 0; j < i; j++) {
            h += _patchItem[j]->boundingRect().width();
            h += 1;
        }

        _patchItem[i]->setImage(_patch[i],
                                _imageItem[i]->histogram().dataMin, _imageItem[i]->histogram().dataMax);
        _patchItem[i]->setTransform(_imageItem[i]->transform());
        _patchItem[i]->refresh();
        _patchItem[i]->show();
        _patchItem[i]->setPos(0, h);
    }

    void SimulCore::run() {
    }
}