//
//  piqSimulCore.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/7/13.
//
//

#ifndef __ParticleGuidedRegistration__piqSimulCore__
#define __ParticleGuidedRegistration__piqSimulCore__

#include <iostream>

#include <QObject>
#include <QVector>

#include "ui_simul2d.h"
#include "piImageDef.h"
#include "piParticle.h"
#include "piPatchTracking.h"

template <class T> class QGraphicsImageItem;
typedef QGraphicsImageItem<pi::RealImage::PixelType> QGraphicsRealImageItem;

class QGraphicsScene;
class QGraphicsPixmapItem;
class QGraphicsEllipseItem;
class QGraphicsRectWidget;
class QAbstractGraphicsShapeItem;

namespace pi {
    class ParticleSystem;
    class ParticleSystemSolver;
}

/**
 * This class is responsible for running and executing registration algorithm.
 *
 */
namespace piq {
    class SimulCore: public QObject {
        Q_OBJECT

    public:
        SimulCore(QWidget* parent = NULL);
        virtual ~SimulCore();
        void setParticleSolver(pi::ParticleSystemSolver* solver);
        void showAuxImage(int, pi::RealImage::Pointer image);

    public slots:
        void setUi(Ui_Simul2D* ui);
        void openImage(int, QString);
        void openLabel(int, QString);
        void hideParticles(bool show);
        
        void startTrackingMode();
        void trackingWidgetMoved(QPointF);
        void trackPatch();
        
        void labelOpacityChanged(int value);
        void updateParticles();
        void updateParticles(int n, pi::ParticleVector& particles);

        void run();


        QGraphicsItem* getImageItem(int n);
        QGraphicsScene* scene(int n);

    private:
        void connectSignals();
        QGraphicsRealImageItem* showImage(int n, pi::RealImage::Pointer image);
        QGraphicsPixmapItem* showLabel(int n, pi::LabelImage::Pointer image);
        void createParticleItems(int id, int n);
        void rectToRegion(QRectF&, pi::RealImage::RegionType&);

    private:
        Ui_Simul2D* _ui;
        QWidget* _parent;

        pi::RealImage::Pointer _image[2];
        pi::LabelImage::Pointer _label[2];


        pi::PatchTracking _tracking;
        
        QGraphicsScene* _scene[2];
        QGraphicsScene* _miniScene;
        
        QGraphicsRealImageItem* _imageItem[2];
        QGraphicsRealImageItem* _auxImageItem[2];
        QGraphicsRealImageItem* _patchItem[2];

        
        QGraphicsPixmapItem* _labelItem[2];
        QVector<QGraphicsEllipseItem*> _particleItem[2];
        QAbstractGraphicsShapeItem* _trackingItem[2];
        QGraphicsRectItem* _trackingRect;
        QGraphicsRectWidget* _trackingWidget;

        pi::ParticleSystemSolver* _solver;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piqSimulCore__) */
